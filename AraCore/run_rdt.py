"""Python reimplementation of AraCore/run_rdt.sh for atom-to-atom mapping via RDT."""

import argparse
import base64
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import List, Optional, Tuple

from unite_mappings import _resolve_reactions_dir


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


def parse_rxn_header(rxn_text: str) -> Tuple[int, int]:
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


def split_rxn_to_mols(rxn_text: str) -> List[str]:
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


def mol_block_to_mdl(mol_block: str) -> str:
    """Convert a raw molecule block (after `$MOL` split) into MDL text.

    In AraCore mode the `M CHG` lines are preserved (unlike the
    MetaCyc variant which strips them).  The function is essentially a
    pass-through because `split_rxn_to_mols` already strips the
    `$MOL` marker line.

    Args:
        mol_block: Text of one molecule block from `split_rxn_to_mols`.

    Returns:
        MDL-format molecule text ready for obabel consumption.
    """
    return mol_block


def parse_mdl_atom_table(mdl_text: str) -> List[Tuple[str, int]]:
    """Extract (element, rdt_atom_index) from V2000 atom lines in an MDL block.

    Replaces `awk '(NF==16){print $4"\\t"$14} (NF==15){print $4"\\t"(0+$13)}'`.
    Lines with 16 whitespace-separated fields carry the atom index in
    field 13 (0-based); lines with 15 fields use field 12 instead.

    Args:
        mdl_text: MDL-format molecule text (as produced by `mol_block_to_mdl`).

    Returns:
        Ordered list of `(element_symbol, global_atom_index)` tuples,
        one per atom in the molecule.
    """
    lines = mdl_text.split("\n")
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


def parse_inchi_atom_order(inchi_text: str) -> List[int]:
    """Parse the InChI auxiliary `/N:` field to recover atom ordering.

    Replaces `grep 'AuxInfo' ... | sed 's/^.*\\/N://; s\\/.*$//; s/,/ /g'`.
    The `/N:` field lists 1-based line numbers in the MDL atom table
    that correspond to InChI's element-wise canonical ordering.

    Args:
        inchi_text: Contents of a `.inchi` file (standard InChI + AuxInfo lines).

    Returns:
        List of 1-based atom-table line numbers in InChI order.
        Defaults to `[1]` when no `/N:` field is present (single atom).
    """
    for line in inchi_text.split("\n"):
        if "AuxInfo" in line and "/N:" in line:
            n_part = re.sub(r"^.*?/N:", "", line)
            n_part = re.sub(r"/.*$", "", n_part)
            indices = [int(x) for x in n_part.split(",")]
            return indices
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


def lookup_species(inchikey: str, table: list) -> List[str]:
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


def load_species_list(text: str) -> List[str]:
    """Load a newline-separated list of species IDs (with compartment tags).

    Used for `from_species_with_cmp` and `to_species_with_cmp` files.

    Args:
        text: File contents with one species ID per line.

    Returns:
        List of stripped, non-empty lines.
    """
    return [line.strip() for line in text.strip().split("\n") if line.strip()]


def find_species_with_cmp(species_no_cmp: str, cmp_list: List[str]) -> Optional[str]:
    """Find all compartmented species IDs that contain the given base ID.

    Replaces `grep "$species_id_without_cmp" from_species_with_cmp`.
    Multiple matches are space-joined, mirroring the bash behaviour where
    `echo $(grep ...)` collapses newlines into spaces.

    Args:
        species_no_cmp: Species ID without compartment (e.g. "M_GAP").
        cmp_list: Species IDs with compartment tags (e.g. `["M_GAP[h]"]`).

    Returns:
        Space-joined string of all matching entries, or `None` if no match.
    """
    matches = [entry for entry in cmp_list if species_no_cmp in entry]
    if matches:
        return " ".join(matches)
    return None


def find_species_with_cmp_multi(species_ids: List[str], cmp_list: List[str]) -> Optional[str]:
    """Find compartmented species IDs matching any of several base IDs.

    When `lookup_species` returns multiple candidates (e.g. both
    `M_Glc` and `M_starch1` share the same InChIKey prefix), the
    bash `grep` searches for all of them against the species list at
    once.  This function replicates that: it tries every candidate
    against the compartmented list and returns all unique hits.

    Args:
        species_ids: Candidate species IDs without compartment.
        cmp_list: Species IDs with compartment tags.

    Returns:
        Space-joined string of all unique matching entries, or `None`.
    """
    all_matches = []
    for sid in species_ids:
        for entry in cmp_list:
            if sid in entry and entry not in all_matches:
                all_matches.append(entry)
    if all_matches:
        return " ".join(all_matches)
    return None


def build_mapping_lines(
    rdt_index: List[Tuple[str, int]],
    inchi_order: List[int],
    species_id: str,
    side: str,
) -> List[str]:
    """Build individual mapping-line entries for one molecule.

    Walks the InChI atom order, tracks an element-wise counter (C#1,
    C#2, ..., N#1, ...), and emits one tab-separated line per atom in the
    format::

        <rdt_atom_index>\\t<from|to>\\t<species>:<element>#<counter><separator>

    where the separator is `=` for reactants (from-side) and `,`
    for products (to-side).

    Replaces the inner `for rdt_line in $inchi_index` loop in the
    bash script.

    Args:
        rdt_index: Per-atom `(element, global_atom_index)` from
            `parse_mdl_atom_table`.
        inchi_order: 1-based atom-table line numbers in InChI order,
            from `parse_inchi_atom_order`.
        species_id: Compartmented species identifier (e.g. "M_GAP[h]").
        side: "from" for reactants, "to" for products.

    Returns:
        List of formatted mapping-line strings.
    """
    mapping_end = "=" if side == "from" else ","
    lines = []
    last_element = None
    atom_counter = 0
    for rdt_line_num in inchi_order:
        element, mapping_index = rdt_index[rdt_line_num - 1]
        if last_element == element:
            atom_counter += 1
        else:
            last_element = element
            atom_counter = 1
        line = f"{mapping_index}\t{side}\t{species_id}:{element}#{atom_counter}{mapping_end}"
        lines.append(line)
    return lines


def assemble_mapping(mapping_lines_text: str) -> str:
    """Assemble the final `mapping.txt` content from individual lines.

    Replaces `sort -n mapping_lines.txt | grep -v ':H#' | cut -f3 |
    tr '\\n' ' ' | sed 's/ //g; s/,$//'`.  Lines are sorted
    numerically by the first column (RDT atom index), hydrogen atoms
    (`:H#`) are dropped, and the third column entries are concatenated
    into a single comma-separated string.

    Args:
        mapping_lines_text: Newline-separated mapping lines as produced
            by `build_mapping_lines`.

    Returns:
        Single-line mapping string without trailing newline, e.g.
        "M_GAP[h]:O#1=M_FBP[h]:O#2,M_GAP[h]:C#1=...".
    """
    lines = mapping_lines_text.strip().split("\n")
    lines = [l for l in lines if l.strip()]
    lines.sort(key=lambda l: int(l.split("\t")[0]))
    parts = []
    for line in lines:
        if ":H#" in line:
            continue
        fields = line.split("\t")
        if len(fields) >= 3:
            parts.append(fields[2])
    return "".join(parts).replace(" ", "").rstrip(",")


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
    cwd: Optional[Path] = None,
) -> RDTResult:
    """Run RDT and return an ``RDTResult`` for interactive / Jupyter use.

    Unlike :func:`run_rdt_java` (which writes into a caller-specified
    directory and returns ``None``), this function returns the contents
    of the generated files.  When *cwd* is ``None`` a temporary
    directory is created automatically.

    Args:
        smiles: Reaction SMILES string (educts>>products).
        rdt_jar: Path to the RDT JAR file.
        cwd: Optional working directory.  Defaults to a temp directory.

    Returns:
        :class:`RDTResult` containing the ``.rxn``, ``.txt`` and ``.png``
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


def obabel_to_inchi(mdl_path: Path, out_path: Path) -> None:
    """Convert an MDL file to InChI with auxiliary info using OpenBabel.

    Replaces `obabel -i mdl ... -o inchi -xa -xT/nochg -O ...`.
    The `-xT/nochg` flag strips charge information from the InChI
    (for canonical atom ordering).  `-xa` requests auxiliary info
    containing the original atom positions needed for mapping.

    Args:
        mdl_path: Path to the input MDL file.
        out_path: Path for the output `.inchi` file.

    Raises:
        SubprocessError: If obabel exits non-zero.
    """
    cmd = ["obabel", "-imdl", str(mdl_path), "-oinchi", "-xa", "-xT/nochg",
           "-O", str(out_path)]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        raise SubprocessError(cmd, result.returncode, result.stdout, result.stderr)


def obabel_to_inchikey(mdl_path: Path, out_path: Path) -> None:
    """Convert an MDL file to an InChIKey using OpenBabel.

    Replaces `obabel -i mdl ... -oinchikey -O ...`.

    Args:
        mdl_path: Path to the input MDL file.
        out_path: Path for the output `.inchikey` file.

    Raises:
        SubprocessError: If obabel exits non-zero.
    """
    cmd = ["obabel", "-imdl", str(mdl_path), "-oinchikey", "-O", str(out_path)]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        raise SubprocessError(cmd, result.returncode, result.stdout, result.stderr)


def postprocess_reaction(rxn_dir: Path) -> bool:
    """Post-process one reaction folder: split RXN → identify species → build mapping.

    Reads the existing `ECBLAST_smiles_AAM.rxn`, splits it into
    individual molecules, runs obabel for InChI/InChIKey generation,
    identifies each molecule's species ID, and assembles the final
    `mapping.txt` and `mapping_lines.txt`.

    Handles edge cases gracefully: empty or missing RXN files produce
    empty mapping files, matching the bash script's behaviour.

    Args:
        rxn_dir: Path to a reaction subfolder inside `reaction_intermediates/`.

    Returns:
        `True` on success, `False` on error (with a message to stderr).
    """
    try:
        rxn_file = rxn_dir / "ECBLAST_smiles_AAM.rxn"
        if not rxn_file.exists() or rxn_file.stat().st_size == 0:
            (rxn_dir / "mapping.txt").write_text("")
            (rxn_dir / "mapping_lines.txt").write_text("")
            return True

        rxn_text = rxn_file.read_text()
        if "$MOL" not in rxn_text:
            (rxn_dir / "mapping.txt").write_text("")
            (rxn_dir / "mapping_lines.txt").write_text("")
            return True

        from_num, to_num = parse_rxn_header(rxn_text)

        mol_blocks = split_rxn_to_mols(rxn_text)

        species_inchikey_text = (rxn_dir / "species_id_inchikey.txt").read_text()
        inchikey_table = load_inchikey_table(species_inchikey_text)

        from_species_text = (rxn_dir / "from_species_with_cmp").read_text()
        to_species_text = (rxn_dir / "to_species_with_cmp").read_text()
        from_species = load_species_list(from_species_text)
        to_species = load_species_list(to_species_text)

        all_mapping_lines = []
        counter = 1

        for i, mol_block in enumerate(mol_blocks):
            mol_num = i + 1
            mol_prefix = f"MOL_{mol_num:02d}"

            mdl_text = mol_block_to_mdl(mol_block)
            mdl_path = rxn_dir / f"{mol_prefix}.mdl"
            mdl_path.write_text(mdl_text)

            rdt_index = parse_mdl_atom_table(mdl_text)
            rdt_index_path = rxn_dir / f"{mol_prefix}.rdt_index"
            rdt_index_path.write_text(
                "\n".join(f"{elem}\t{idx}" for elem, idx in rdt_index) + "\n"
            )

            inchi_path = rxn_dir / f"{mol_prefix}.inchi"
            inchikey_path = rxn_dir / f"{mol_prefix}.inchikey"

            obabel_to_inchi(mdl_path, inchi_path)
            obabel_to_inchikey(mdl_path, inchikey_path)

            if not inchikey_path.exists() or inchikey_path.stat().st_size == 0:
                counter += 1
                continue

            inchikey = inchikey_path.read_text().strip()
            if not inchikey:
                counter += 1
                continue

            matches = lookup_species(inchikey, inchikey_table)

            if not matches:
                counter += 1
                continue

            if counter <= from_num:
                species_id = find_species_with_cmp_multi(matches, from_species)
                mapping_side = "from"
            else:
                species_id = find_species_with_cmp_multi(matches, to_species)
                mapping_side = "to"

            if species_id is None:
                counter += 1
                continue

            species_id_path = rxn_dir / f"{mol_prefix}.species_id"
            species_id_path.write_text(species_id + "\n")

            inchi_text = inchi_path.read_text()
            inchi_order = parse_inchi_atom_order(inchi_text)

            lines = build_mapping_lines(rdt_index, inchi_order, species_id, mapping_side)
            all_mapping_lines.extend(lines)

            counter += 1

        mapping_lines_path = rxn_dir / "mapping_lines.txt"
        mapping_lines_path.write_text("\n".join(all_mapping_lines) + "\n")

        mapping_text = assemble_mapping("\n".join(all_mapping_lines))
        (rxn_dir / "mapping.txt").write_text(mapping_text)

        return True
    except SubprocessError as e:
        print(f"Error processing {rxn_dir.name}: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error processing {rxn_dir.name}: {e}", file=sys.stderr)
        return False


def process_reaction(rxn_dir: Path, rdt_jar: Path) -> bool:
    """Run the full pipeline for a single reaction: RDT + postprocessing.

    Calls `run_rdt_java` to generate the `.rxn` file, then
    delegates to `postprocess_reaction` for splitting, species
    identification, and mapping assembly.

    Args:
        rxn_dir: Path to a reaction subfolder containing `rxn.smiles`.
        rdt_jar: Path to the RDT JAR file.

    Returns:
        `True` on success, `False` if `rxn.smiles` is missing,
        RDT fails, or postprocessing fails.
    """
    rxn_file = rxn_dir / "ECBLAST_smiles_AAM.rxn"
    rxn_smiles_path = rxn_dir / "rxn.smiles"
    if not rxn_smiles_path.exists():
        return False

    smiles = rxn_smiles_path.read_text().strip()
    try:
        run_rdt_java(smiles, rdt_jar, rxn_dir)
    except SubprocessError as e:
        print(f"Error processing {rxn_dir.name}: {e}", file=sys.stderr)
        return False

    return postprocess_reaction(rxn_dir)


def main():
    """CLI entry point: iterate over all reaction folders and run the pipeline.

    Supports two modes via flags:

    * **Default** (`python run_rdt.py`): run RDT Java on each reaction
      then postprocess.
    * `--postprocess-only`: skip the RDT step and re-derive mapping
      files from existing `.rxn` output.
    """
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
            if postprocess_reaction(rxn_folder):
                success += 1
        else:
            if process_reaction(rxn_folder, rdt_jar):
                success += 1

    print(f"\nProcessed {success}/{total} reactions successfully", file=sys.stderr)


if __name__ == "__main__":
    main()
