"""Python reimplementation of AraCore/run_rdt.sh for atom-to-atom mapping via RDT."""

import argparse
import base64
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from html import escape
from pathlib import Path


from unite_mappings import _resolve_reactions_dir


_N_FIELD = re.compile(
    r"""
    /N:       # atom-ordering field marker
    ([^/]+)   # capture: comma-separated 1-based indices
    """,
    re.VERBOSE,
)

_COMPARTMENT_RE = re.compile(r"\[([^\]]+)\]")


def extract_compartment(species_id: str) -> str:
    """Extract the compartment tag from a compartmented species ID.

    "M_GAP[h]" -> "h", "M_Glc[c]" -> "c".
    Returns "" for multi-species strings (containing spaces) or strings without a bracket tag.
    """
    if " " in species_id:
        return ""
    m = _COMPARTMENT_RE.search(species_id)
    return m.group(1) if m else ""


class SubprocessError(Exception):
    """Raised when a subprocess exits non-zero, capturing full output."""
    def __init__(self, cmd, returncode, stdout, stderr):
        self.cmd = cmd
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        cmd_str = " ".join(str(c) for c in cmd)
        stdout_str = stdout.decode(errors="replace") if isinstance(stdout, bytes) else (stdout or "")
        stderr_str = stderr.decode(errors="replace") if isinstance(stderr, bytes) else (stderr or "")
        super().__init__(
            f"Command '{cmd_str}' exited with code {returncode}\n"
            f"stdout:\n{stdout_str}\n"
            f"stderr:\n{stderr_str}"
        )


@dataclass(frozen=True)
class MappingEntry:
    """A single atom's role in the atom-to-atom mapping between reaction sides.

    RDT assigns each atom a global index across the whole reaction.
    InChI provides a canonical element-wise ordering (all C's, then all N's etc.).
    A MappingEntry connects one atom's RDT index to a species-specific label
    of the form metabolite_name[compartment]:element#inchi_element_wise_atom_position,
    e.g. `M_GAP[h]:C#1`

    Example:
        >>> MappingEntry(2, "from", "M_GAP[h]", "h", "C", 1)
        MappingEntry(rdt_atom_index=2, side='from', species_id='M_GAP[h]',
                     compartment='h', element='C', element_index=1)
        >>> MappingEntry(2, "from", "M_GAP[h]", "h", "C", 1).label
        'M_GAP[h]:C#1='

    Attributes:
        rdt_atom_index: Global atom index assigned by RDT across the reaction.
        side: "from" (reactant) or "to" (product). Determines the
            separator in `label`: `=` for from, `,` for to.
        species_id: Compartmented species identifier, e.g. "M_GAP[h]".
        compartment: Subcellular compartment tag, e.g. "h" for chloroplast,
            "c" for cytosol, "m" for mitochondria.
        element: Chemical element symbol, e.g. "C" or "O".
        element_index: 1-based counter within this element (C#1, C#2, N#1, ...).
    """
    rdt_atom_index: int
    side: str
    species_id: str
    compartment: str
    element: str
    element_index: int

    @property
    def label(self) -> str:
        sep = "=" if self.side == "from" else ","
        return f"{self.species_id}:{self.element}#{self.element_index}{sep}"


@dataclass
class MoleculeProcessingResult:
    """Per-molecule result of processing one molecule block.

    Attributes:
        mol_num: 1-based molecule number in the reaction.
        rdt_index: (element, rdt_atom_index) in MDL V2000 atom-table order.
        inchi_order: MDL line numbers in InChI canonical order.
        entries: MappingEntry objects in InChI order (same length as inchi_order).
        species_id: Identified species with compartment, e.g. "M_GAP[h]".
        compartment: Subcellular compartment tag.
        side: "from" (reactant) or "to" (product).
        inchi: InChI string (first line from obabel output, no AuxInfo).
    """
    mol_num: int
    rdt_index: list[tuple[str, int]]
    inchi_order: list[int]
    entries: list[MappingEntry]
    species_id: str
    compartment: str
    side: str
    inchi: str = ""


@dataclass
class ReactionProcessingResult:
    """Complete structured result of processing one reaction folder.

    Attributes:
        rxn_name: Folder name (e.g. "FBPA_h").
        smiles: Reaction SMILES string.
        from_num: Number of reactant molecules.
        to_num: Number of product molecules.
        molecules: Per-molecule processing results.
        mapping_lines_text: Full mapping_lines.txt content.
        mapping_text: Final mapping.txt content.
    """
    rxn_name: str
    smiles: str
    from_num: int
    to_num: int
    molecules: list[MoleculeProcessingResult]
    mapping_lines_text: str
    mapping_text: str

    @property
    def all_entries(self) -> list[MappingEntry]:
        """Flattened list of all MappingEntry objects across all molecules."""
        return [e for mol in self.molecules for e in mol.entries]


def parse_rxn_header(rxn_text: str) -> tuple[int, int]:
    """Extract the number of reactant and product molecules from an MDL RXN header.

    Replaces `grep -m1 -B1 '$MOL' ... | head -n 1` combined with
    `head -c3` / `tail -c+4` in the bash script.

    Args:
        rxn_text: Full contents of an `ECBLAST_smiles_AAM.rxn` file.

    Returns:
        `(from_num, to_num)`: the count of reactant and product molecules.

    Raises:
        ValueError: If no `$MOL` marker is found or the header line cannot be parsed.
    """
    mol_pos = rxn_text.find("$MOL")
    if mol_pos < 0:
        raise ValueError("No $MOL found in RXN text")
    pre_mol = rxn_text[:mol_pos]
    lines = pre_mol.split("\n")
    for line in reversed(lines):
        parts = line.split()
        if len(parts) == 2:
            try:
                return int(parts[0]), int(parts[1])
            except ValueError:
                continue
    raise ValueError("Could not parse RXN header counts")


def split_rxn_to_mols(rxn_text: str) -> list[str]:
    """Split an MDL RXN file into individual molecule blocks.

    Replaces `csplit -f MOL_ ... '/$MOL/' {*}`: each returned block
    is the text between consecutive `$MOL` markers (the header line
    `M00001`, atom table, bonds, `M  END`).

    Args:
        rxn_text: Full contents of an `ECBLAST_smiles_AAM.rxn` file.

    Returns:
        List of molecule block strings (one per molecule in the reaction).
    """
    parts = rxn_text.split("$MOL")
    blocks = []
    for part in parts[1:]:
        block = part
        if block.startswith("\n"):
            block = block[1:]
        blocks.append(block)
    return blocks


def parse_mdl_atom_table(mol_block: str) -> list[tuple[str, int]]:
    """Extract (element, rdt_atom_index) from V2000 atom lines in an MDL block.

    Replaces `awk '(NF==16){print $4"\\t"$14} (NF==15){print $4"\\t"(0+$13)}'`.
    Lines with 16 whitespace-separated fields carry the atom index in
    field 13 (0-based); lines with 15 fields use field 12 instead.

    Args:
        mol_block: MDL-format molecule text (as produced by `split_rxn_to_mols`).

    Returns:
        Ordered list of `(element_symbol, global_atom_index)` tuples,
        one per atom in the molecule.
    """
    lines = mol_block.split("\n")
    atoms = []
    for line in lines:
        fields = line.split()
        nf = len(fields)
        if nf == 16:
            element = fields[3]
            atom_index = int(fields[13])
            atoms.append((element, atom_index))
        elif nf == 15:
            element = fields[3]
            atom_index = int(fields[12])
            atoms.append((element, atom_index))
    return atoms


def parse_inchi_atom_order(inchi_text: str) -> list[int]:
    """Parse the InChI auxiliary `/N:` field to recover atom ordering.

    Replaces `grep 'AuxInfo' ... | sed 's/^.*\\/N://; s\\/.*$//; s/,/ /g'`.
    The `/N:` field lists 1-based line numbers in the MDL atom table
    that correspond to InChI's element-wise canonical ordering.

    Args:
        inchi_text: Contents of a `.inchi` file (standard InChI + AuxInfo lines).

    Returns:
        List of 1-based atom-table line numbers in InChI order.
        Returns `[1]` when no `/N:` field is present
        (correct for single-atom molecules like H+).

    Example:
        >>> parse_inchi_atom_order(
        ...     "InChI=1S/C3H7O6P/c4-1-3(5)2-9-10(6,7)8/h1,3,5H,2H2,(H2,6,7,8)/t3-/m0/s1\\n"
        ...     "AuxInfo=1/1/N:2,5,3,1,4,8,9,10,6,7/E:(6,7,8)/it:im/rA:10OCCOCOPOO-O-/rB:d1;s2;N3;s3;s5;s6;d7;s7;s7;/rC:-6.3212,3.8505,0;-5.0221,3.1005,0;-3.7231,3.8505,0;-3.7231,5.3505,0;-2.424,3.1005,0;-1.125,3.8505,0;.174,3.1005,0;-.576,1.8014,0;1.4731,2.3505,0;.924,4.3995,0;\\n"
        ... )
        [2, 5, 3, 1, 4, 8, 9, 10, 6, 7]
        >>> parse_inchi_atom_order("InChI=1S/H2O/h1H2\\n")
        [1]
    """
    for line in inchi_text.splitlines():
        if "AuxInfo" in line and "/N:" in line:
            match = _N_FIELD.search(line)
            if match:
                return [int(x) for x in match.group(1).split(",")]
    return [1]


def load_inchikey_table(text: str) -> list:
    """Parse `species_id_inchikey.txt` into a lookup table.

    Each line is `<species_id>\\t<inchikey>`.  The species ID may
    contain `_DASH_` placeholders that are replaced with literal `-`.

    Args:
        text: Full contents of a `species_id_inchikey.txt` file.

    Returns:
        List of `(inchikey, species_id_without_compartment)` tuples,
        preserving insertion order so that first/last semantics match the
        bash `grep` + `cut` pipeline.
    """
    table = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            species_id = parts[0].strip().replace("_DASH_", "-")
            inchikey = parts[1].strip()
            table.append((inchikey, species_id))
    return table


def lookup_species(inchikey: str, table: list) -> list[str]:
    """Look up species IDs by InChIKey with a 14-character prefix fallback.

    Replaces the bash `grep "$(cat ...inchikey)" species_id_inchikey.txt`
    exact match, followed by `grep "$(head -c14 ...inchikey)" ...` when
    the exact match fails.  Multiple matches can occur when different
    species share the same connectivity hash (first 14 chars).

    Args:
        inchikey: InChIKey string (typically 27 characters).
        table: Table from `load_inchikey_table`.

    Returns:
        List of matching species IDs (without compartment suffix).
        Empty list if no match is found at all.
    """
    matches = []
    for ref_key, species_id in table:
        if ref_key == inchikey:
            matches.append(species_id)
    if matches:
        return matches
    prefix = inchikey[:14]
    for ref_key, species_id in table:
        if ref_key[:14] == prefix:
            matches.append(species_id)
    return matches


def load_species_list(text: str) -> list[str]:
    """Load a newline-separated list of species IDs (with compartment tags).

    Used for `from_species_with_cmp` and `to_species_with_cmp` files.

    Args:
        text: File contents with one species ID per line.

    Returns:
        List of stripped, non-empty lines.
    """
    return [line.strip() for line in text.strip().split("\n") if line.strip()]


def find_species_with_cmp(species_no_cmp: str, cmp_list: list[str]) -> str | None:
    """Find compartmented species IDs whose base name matches exactly.

    Args:
        species_no_cmp: Species ID without compartment (e.g. "M_GAP").
        cmp_list: Species IDs with compartment tags (e.g. `["M_GAP[h]"]`).

    Returns:
        Space-joined string of all matching entries, or `None` if no match.
    """
    prefix = species_no_cmp + "["
    matches = [entry for entry in cmp_list if entry.startswith(prefix)]
    if matches:
        return " ".join(matches)
    return None


def find_species_with_cmp_multi(species_ids: list[str], cmp_list: list[str]) -> str | None:
    """Find compartmented species IDs matching any of several base IDs via
    exact prefix match.

    When `lookup_species` returns multiple candidates (e.g. both
    `M_Glc` and `M_starch1` share the same InChIKey), each candidate
    is matched against the compartmented list using `startswith`
    (e.g. "M_Glc" matches "M_Glc[c]" but not "M_Glc-SeA[c]").
    All unique matches are space-joined.

    Args:
        species_ids: Candidate species IDs without compartment.
        cmp_list: Species IDs with compartment tags.

    Returns:
        Space-joined string of all unique matching entries, or `None`.
    """
    all_matches = []
    for sid in species_ids:
        prefix = sid + "["
        for entry in cmp_list:
            if entry.startswith(prefix) and entry not in all_matches:
                all_matches.append(entry)
    if all_matches:
        return " ".join(all_matches)
    return None


def build_mapping_entries(
    rdt_index: list[tuple[str, int]],
    inchi_order: list[int],
    species_id: str,
    compartment: str,
    side: str,
) -> list[MappingEntry]:
    """Build structured mapping entries for one molecule.

    Converts RDT's global atom indices and InChI's canonical element ordering
    into species-labeled atom references suitable for assembly into the final
    mapping string.

    For each atom in InChI order, creates a `MappingEntry` with a per-element
    counter (C#1, C#2, ..., N#1, ...).

    Example:
        >>> rdt_index = [("C", 2), ("C", 5), ("O", 1)]
        >>> inchi_order = [1, 2, 3]  # InChI: C, C, O
        >>> entries = build_mapping_entries(rdt_index, inchi_order, "M_X[h]", "h", "from")
        >>> entries[0].label
        'M_X[h]:C#1='
        >>> entries[2].label
        'M_X[h]:O#1='

    Args:
        rdt_index: Per-atom `(element, global_atom_index)` from
            `parse_mdl_atom_table`.
        inchi_order: 1-based atom-table line numbers in InChI order,
            from `parse_inchi_atom_order`.
        species_id: Compartmented species identifier (e.g. `"M_GAP[h]"`).
        compartment: Subcellular compartment (e.g. `"h"`).  Empty string
            when species_id is multi-species.
        side: `"from"` for reactants, `"to"` for products.

    Returns:
        List of `MappingEntry` objects, one per atom.
    """
    counts: Counter[str] = Counter()
    entries: list[MappingEntry] = []
    for rdt_line_num in inchi_order:
        element, mapping_index = rdt_index[rdt_line_num - 1]
        counts[element] += 1
        entries.append(MappingEntry(
            rdt_atom_index=mapping_index,
            side=side,
            species_id=species_id,
            compartment=compartment,
            element=element,
            element_index=counts[element],
        ))
    return entries


def assemble_mapping(entries: list[MappingEntry]) -> str:
    """Assemble the final atom-to-atom mapping string from all entries.

    Takes mapping entries from all molecules in a reaction, sorts them by
    RDT's global atom index, drops hydrogen atoms (which are typically
    not matched between sides), and concatenates the remaining labels.

    The result is a single string encoding the full mapping:

        M_GAP[h]:O#1=M_FBP[h]:O#2,M_GAP[h]:C#1=M_FBP[h]:C#5,...

    Reactant atoms are terminated with `=`, product atoms with `,`
    (from `MappingEntry.label`).  These link across the `=` sign
    to show which reactant atom maps to which product atom.

    Args:
        entries: `MappingEntry` objects from all molecules, collected
            across multiple calls to `build_mapping_entries`.

    Returns:
        Single-line mapping string.  Empty string if all entries are hydrogen
        or the list is empty.
    """
    non_h = [e for e in entries if e.element != "H"]
    non_h.sort(key=lambda e: e.rdt_atom_index)
    return "".join(e.label for e in non_h).replace(" ", "").rstrip(",")


def run_rdt_java(smiles: str, rdt_jar: Path, cwd: Path) -> None:
    """Run the RDT Java tool to generate an atom-atom mapped `.rxn` file.

    Writes `ECBLAST_smiles_AAM.rxn` (and associated `.png`/`.txt`)
    into *cwd*.  Requires Java and the RDT JAR (v2.5.0).

    Args:
        smiles: Reaction SMILES string (educts>>products).
        rdt_jar: Path to the RDT JAR file.
        cwd: Working directory: RDT writes output here.

    Raises:
        SubprocessError: If the Java process exits non-zero.
    """
    cmd = [
        "java", "-jar", str(rdt_jar),
        "-Q", "SMI", "-q", smiles,
        "-g", "-c", "-b", "-j", "AAM", "-f", "TEXT",
    ]
    result = subprocess.run(cmd, cwd=str(cwd), capture_output=True)

    rxn_file = cwd / "ECBLAST_smiles_AAM.rxn"
    if rxn_file.exists() and rxn_file.stat().st_size > 0:
        return

    if result.returncode != 0:
        raise SubprocessError(cmd, result.returncode, result.stdout, result.stderr)
    raise SubprocessError(cmd, 0, b"", b"RDT exited 0 but produced no output files")


@dataclass
class RDTResult:
    """Holds the output of a successful RDT run for display in Jupyter."""

    rxn: str
    txt: str
    png: bytes

    def __repr__(self) -> str:
        parts = []
        if self.txt:
            parts.append("=== Text ===")
            parts.append(self.txt)
        if self.rxn:
            parts.append("=== RXN ===")
            parts.append(self.rxn)
        return "\n".join(parts)

    def _repr_html_(self) -> str:
        png_b64 = base64.b64encode(self.png).decode()
        return (
            '<div style="font-family: monospace;">'
            f'<img src="data:image/png;base64,{png_b64}" '
            'style="max-width:100%;" />'
            '<details open><summary>Text output</summary>'
            f'<pre style="font-size:0.85em;">{escape(self.txt)}</pre>'
            '</details>'
            '<details><summary>RXN file</summary>'
            f'<pre style="font-size:0.85em;">{escape(self.rxn)}</pre>'
            '</details>'
            '</div>'
        )


def run_rdt_jupyter(
    smiles: str,
    rdt_jar: Path,
    cwd: Path | None = None,
) -> RDTResult:
    """Run RDT and return an `RDTResult` for interactive / Jupyter use.

    Unlike `run_rdt_java` (which writes into a caller-specified
    directory and returns `None`), this function returns the contents
    of the generated files.  When *cwd* is `None` a temporary
    directory is created automatically.

    Args:
        smiles: Reaction SMILES string (educts>>products).
        rdt_jar: Path to the RDT JAR file.
        cwd: Optional working directory.  Defaults to a temp directory.

    Returns:
        `RDTResult` containing the `.rxn`, `.txt` and `.png`
        output produced by RDT.

    Raises:
        SubprocessError: If RDT fails to produce output files.
    """
    if cwd is None:
        cwd = Path(tempfile.mkdtemp())

    cmd = [
        "java", "-jar", str(rdt_jar),
        "-Q", "SMI", "-q", smiles,
        "-g", "-c", "-b", "-j", "AAM", "-f", "TEXT",
    ]
    result = subprocess.run(cmd, cwd=str(cwd), capture_output=True)

    rxn_path = cwd / "ECBLAST_smiles_AAM.rxn"
    txt_path = cwd / "ECBLAST_smiles_AAM.txt"
    png_path = cwd / "ECBLAST_smiles_AAM.png"

    if rxn_path.exists() and rxn_path.stat().st_size > 0:
        return RDTResult(
            rxn=rxn_path.read_text(),
            txt=txt_path.read_text() if txt_path.exists() else "",
            png=png_path.read_bytes() if png_path.exists() else b"",
        )

    raise SubprocessError(cmd, result.returncode, result.stdout, result.stderr)


def obabel_to_inchi(mol_block: str) -> str:
    """Convert a molecule block string to InChI with auxiliary info using OpenBabel.

    Replaces `obabel -i mdl ... -o inchi -xa -xT/nochg -O ...`.
    The `-xT/nochg` flag strips charge information from the InChI
    (for canonical atom ordering).  `-xa` requests auxiliary info
    containing the original atom positions needed for mapping.

    Args:
        mol_block: molecule block string, as produced by `split_rxn_to_mols`

    Raises:
        SubprocessError: If obabel exits non-zero.
    
    Returns:
        InChI string with auxiliary info
    """
    cmd = ["obabel", "-imdl", "-", "-oinchi", "-xa", "-xT/nochg"]
    result = subprocess.run(cmd, capture_output=True, input=mol_block.encode())
    if result.returncode != 0:
        raise SubprocessError(cmd, result.returncode, result.stdout, result.stderr)
    return result.stdout.decode()


def obabel_to_inchikey(mol_block: str) -> str:
    """Convert a molecule block string to an InChIKey using OpenBabel.

    Replaces `obabel -i mdl ... -oinchikey -O ...`.

    Args:
        mol_block: molecule block string, as produced by `split_rxn_to_mols`

    Raises:
        SubprocessError: If obabel exits non-zero.

    Returns:
        InChIKey string
    """
    cmd = ["obabel", "-imdl", "-", "-oinchikey"]
    result = subprocess.run(cmd, capture_output=True, input=mol_block.encode())
    if result.returncode != 0:
        raise SubprocessError(cmd, result.returncode, result.stdout, result.stderr)
    return result.stdout.decode().strip()

def process_reaction_data(rxn_dir: Path) -> ReactionProcessingResult:
    """Process a reaction folder and return structured data (no file writes).

    Reads the RXN file and species-lookup files, runs obabel for InChI
    generation, identifies species, and builds the full mapping.  Returns
    a structured `ReactionProcessingResult` that can be used for display
    or written to disk.

    Args:
        rxn_dir: Path to a reaction subfolder (must contain
            ``ECBLAST_smiles_AAM.rxn`` and species-lookup files).

    Returns:
        ``ReactionProcessingResult`` on success.

    Raises:
        FileNotFoundError: If the RXN file is missing or empty.
        ValueError: If the RXN file has no ``$MOL`` marker.
        SubprocessError: If obabel fails.
    """
    rxn_file = rxn_dir / "ECBLAST_smiles_AAM.rxn"
    if not rxn_file.exists() or rxn_file.stat().st_size == 0:
        raise FileNotFoundError(
            f"RXN file missing or empty in {rxn_dir}"
        )

    rxn_text = rxn_file.read_text()
    if "$MOL" not in rxn_text:
        raise ValueError(
            f"No $MOL marker in RXN file in {rxn_dir}"
        )

    smiles_path = rxn_dir / "rxn.smiles"
    smiles = smiles_path.read_text().strip() if smiles_path.exists() else ""

    from_num, to_num = parse_rxn_header(rxn_text)
    mol_blocks = split_rxn_to_mols(rxn_text)

    species_inchikey_text = (rxn_dir / "species_id_inchikey.txt").read_text()
    inchikey_table = load_inchikey_table(species_inchikey_text)

    from_species_text = (rxn_dir / "from_species_with_cmp").read_text()
    to_species_text = (rxn_dir / "to_species_with_cmp").read_text()
    from_species = load_species_list(from_species_text)
    to_species = load_species_list(to_species_text)

    molecules: list[MoleculeProcessingResult] = []
    counter = 1

    for i, mol_block in enumerate(mol_blocks):
        rdt_index = parse_mdl_atom_table(mol_block)
        inchi_text = obabel_to_inchi(mol_block)
        inchi_string = inchi_text.split("\n")[0] if inchi_text else ""
        inchikey = obabel_to_inchikey(mol_block)

        if not inchikey:
            counter += 1
            continue

        matches = lookup_species(inchikey, inchikey_table)
        if not matches:
            counter += 1
            continue

        if counter <= from_num:
            species_id = find_species_with_cmp_multi(matches, from_species)
            side = "from"
        else:
            species_id = find_species_with_cmp_multi(matches, to_species)
            side = "to"

        if species_id is None:
            counter += 1
            continue

        compartment = extract_compartment(species_id)
        inchi_order = parse_inchi_atom_order(inchi_text)
        entries = build_mapping_entries(rdt_index, inchi_order, species_id, compartment, side)

        molecules.append(MoleculeProcessingResult(
            mol_num=i + 1,
            rdt_index=rdt_index,
            inchi_order=inchi_order,
            entries=entries,
            species_id=species_id,
            compartment=compartment,
            side=side,
            inchi=inchi_string,
        ))
        counter += 1

    all_entries = [e for mol in molecules for e in mol.entries]
    mapping_lines_text = (
        "\n".join(f"{e.rdt_atom_index}\t{e.side}\t{e.label}" for e in all_entries)
        + ("\n" if all_entries else "")
    )
    mapping_text = assemble_mapping(all_entries)

    return ReactionProcessingResult(
        rxn_name=rxn_dir.name,
        smiles=smiles,
        from_num=from_num,
        to_num=to_num,
        molecules=molecules,
        mapping_lines_text=mapping_lines_text,
        mapping_text=mapping_text,
    )


def postprocess_reaction(rxn_dir: Path) -> tuple[bool, str, str]:
    """Post-process one reaction folder: split RXN -> identify species -> build mapping.

    Thin wrapper around ``process_reaction_data`` that writes the output
    files (``mapping.txt`` and ``mapping_lines.txt``) to disk and returns
    flat strings for backward compatibility.

    On failure, prints a message to stderr and returns ``(False, "", "")``
    without writing any output files.

    Args:
        rxn_dir: Path to a reaction subfolder inside ``reaction_intermediates/``.

    Returns:
        ``(success, mapping_lines_text, mapping_text)``.
    """
    try:
        result = process_reaction_data(rxn_dir)
    except Exception as e:
        print(f"Error processing {rxn_dir.name}: {e}", file=sys.stderr)
        return False, "", ""
    (rxn_dir / "mapping.txt").write_text(result.mapping_text)
    (rxn_dir / "mapping_lines.txt").write_text(result.mapping_lines_text)
    return True, result.mapping_lines_text, result.mapping_text


def process_reaction(rxn_dir: Path, rdt_jar: Path) -> tuple[bool, str, str]:
    """Run the full pipeline for a single reaction: RDT + postprocessing.

    Calls `run_rdt_java` to generate the `.rxn` file, then
    delegates to `postprocess_reaction` for splitting, species
    identification, and mapping assembly.

    Args:
        rxn_dir: Path to a reaction subfolder containing `rxn.smiles`.
        rdt_jar: Path to the RDT JAR file.

    Returns:
        Same as `postprocess_reaction`: `tuple[bool, str, str]`
        `(success, mapping_lines_text, mapping_text)`.
        Returns `(False, "", "")` if `rxn.smiles` is missing
        or RDT fails.
    """
    rxn_file = rxn_dir / "ECBLAST_smiles_AAM.rxn"
    rxn_smiles_path = rxn_dir / "rxn.smiles"
    if not rxn_smiles_path.exists():
        return False, "", ""

    smiles = rxn_smiles_path.read_text().strip()
    try:
        run_rdt_java(smiles, rdt_jar, rxn_dir)
    except SubprocessError as e:
        print(f"Error processing {rxn_dir.name}: {e}", file=sys.stderr)
        return False, "", ""

    return postprocess_reaction(rxn_dir)


def main():
    """CLI entry point: iterate over all reaction folders and run the pipeline."""
    _script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Run RDT atom mapping pipeline for AraCore reactions"
    )
    parser.add_argument(
        "--rdt-jar",
        type=Path,
        default=Path(os.environ.get(
            "RDT_JAR",
            str(_script_dir.parent / "rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar"),
        )),
        help="Path to RDT JAR file",
    )
    parser.add_argument(
        "--reactions-dir",
        type=Path,
        default=_script_dir / "reaction_intermediates.zip",
        help="Directory or .zip archive containing reaction subfolders",
    )
    parser.add_argument(
        "--postprocess-only",
        action="store_true",
        help="Only postprocess existing RDT output (skip RDT Java step)",
    )
    args = parser.parse_args()

    rdt_jar = args.rdt_jar.resolve()
    try:
        reactions_dir = _resolve_reactions_dir(args.reactions_dir)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    rxn_folders = sorted(reactions_dir.iterdir())
    total = len(rxn_folders)
    success = 0

    for i, rxn_folder in enumerate(rxn_folders, 1):
        if not rxn_folder.is_dir():
            continue
        print(rxn_folder.name)

        if args.postprocess_only:
            ok, _, _ = postprocess_reaction(rxn_folder)
            if ok:
                success += 1
        else:
            ok, _, _ = process_reaction(rxn_folder, rdt_jar)
            if ok:
                success += 1

    print(f"\nProcessed {success}/{total} reactions successfully", file=sys.stderr)


if __name__ == "__main__":
    main()
