import hashlib
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

import run_rdt

ARACORE_DIR = Path(__file__).resolve().parent.parent
REACTIONS_ZIP = ARACORE_DIR / "reaction_intermediates.zip"
REACTIONS_DIR = ARACORE_DIR / "reaction_intermediates"
GOLDEN_DIR = Path(__file__).resolve().parent / "golden"
RDT_JAR = Path(os.environ.get(
    "RDT_JAR", ARACORE_DIR.parent / "rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar"
))

KNOWN_PYTHON_BASH_DIFFS = {"DPE12_h", "DPE2_c", "OrnAT_m"}


def _extract_rxn_from_zip(rxn_name, dest_dir):
    prefix = f"reaction_intermediates/{rxn_name}/"
    with zipfile.ZipFile(REACTIONS_ZIP, "r") as z:
        for name in z.namelist():
            if name.startswith(prefix):
                z.extract(name, dest_dir)
    return dest_dir / "reaction_intermediates" / rxn_name


def _strip_generated_files(rxn_dir):
    for pattern in ["mapping*", "MOL_*.mdl", "MOL_*.inchi",
                    "MOL_*.inchikey", "MOL_*.rdt_index", "MOL_*.species_id"]:
        for f in rxn_dir.glob(pattern):
            f.unlink()


def _all_rxn_names_from_zip():
    with zipfile.ZipFile(REACTIONS_ZIP, "r") as z:
        return sorted({
            n.split("/")[1] for n in z.namelist()
            if n.startswith("reaction_intermediates/")
            and n.count("/") == 2
            and n.endswith("/")
        })


def _run_rdt(smiles, rdt_jar, cwd):
    cmd = [
        "java", "-jar", str(rdt_jar),
        "-Q", "SMI", "-q", smiles,
        "-g", "-c", "-b", "-j", "AAM", "-f", "TEXT",
    ]
    result = subprocess.run(cmd, cwd=str(cwd), capture_output=True)
    if result.returncode != 0:
        raise run_rdt.SubprocessError(cmd, result.returncode, result.stdout, result.stderr)


@pytest.mark.integration
def test_postprocess_single_reaction(tmp_path):
    rxn_dir = _extract_rxn_from_zip("FBPA_h", tmp_path)
    _strip_generated_files(rxn_dir)

    success, _, mapping_text = run_rdt.postprocess_reaction(rxn_dir)
    assert success is True

    golden = GOLDEN_DIR / "FBPA_h.mapping.txt"
    assert golden.exists()
    expected = golden.read_text()
    actual = (rxn_dir / "mapping.txt").read_text()
    assert actual == expected


RXN_NAMES = _all_rxn_names_from_zip()


@pytest.mark.integration
@pytest.mark.parametrize("rxn_name", RXN_NAMES)
def test_python_matches_bash_for_one_reaction(tmp_path, rxn_name):
    rxn_dir = _extract_rxn_from_zip(rxn_name, tmp_path / rxn_name)
    bash_mapping = (rxn_dir / "mapping.txt").read_text()
    bash_mapping_lines_path = rxn_dir / "mapping_lines.txt"
    if bash_mapping_lines_path.exists():
        bash_mapping_lines = bash_mapping_lines_path.read_text()
        _strip_generated_files(rxn_dir)
    else:
        bash_mapping_lines = None
        _strip_generated_files(rxn_dir)

    success, mapping_lines_text, mapping_text = run_rdt.postprocess_reaction(rxn_dir)
    py_mapping = (rxn_dir / "mapping.txt").read_text()
    assert mapping_text == py_mapping

    if rxn_name in KNOWN_PYTHON_BASH_DIFFS:
        assert py_mapping != bash_mapping, (
            f"Expected known Python-vs-bash difference for {rxn_name}, "
            f"but mappings now match (bug may have been fixed)"
        )
    else:
        assert py_mapping == bash_mapping, (
            f"Unexpected Python-vs-bash mismatch for {rxn_name}"
        )

    py_mapping_lines = (rxn_dir / "mapping_lines.txt").read_text()
    assert mapping_lines_text == py_mapping_lines

    if bash_mapping_lines is not None and rxn_name not in KNOWN_PYTHON_BASH_DIFFS:
        assert py_mapping_lines == bash_mapping_lines, (
            f"Unexpected Python-vs-bash mapping_lines mismatch for {rxn_name}"
        )


@pytest.mark.integration
def test_postprocess_reaction_is_deterministic(tmp_path):
    rxn_name = "AspAT_h"
    outputs = set()
    for run_idx in range(3):
        rxn_dir = _extract_rxn_from_zip(rxn_name, tmp_path / f"run{run_idx}")
        _strip_generated_files(rxn_dir)

        run_rdt.postprocess_reaction(rxn_dir)
        outputs.add((rxn_dir / "mapping.txt").read_text())
        shutil.rmtree(rxn_dir.parent.parent, ignore_errors=True)

    assert len(outputs) == 1


@pytest.mark.integration
@pytest.mark.parametrize("rxn_name", RXN_NAMES)
def test_mapping_txt_consistent_with_mapping_lines(tmp_path, rxn_name):
    rxn_dir = _extract_rxn_from_zip(rxn_name, tmp_path / rxn_name)
    _strip_generated_files(rxn_dir)

    success, mapping_lines_text, mapping_text = run_rdt.postprocess_reaction(rxn_dir)
    assert success is True

    reconstructed = run_rdt.assemble_mapping(mapping_lines_text)
    assert reconstructed == mapping_text, (
        f"mapping.txt is not derivable from mapping_lines.txt for {rxn_name}"
    )


def _bash_pipeline_assemble(mapping_lines_text):
    """Reproduce bash: sort -n mapping_lines.txt | grep -v ':H#' | cut -f3 | tr '\\n' ' ' | sed 's/ //g; s/,$//'"""
    lines = [l for l in mapping_lines_text.strip().split("\n") if l.strip()]
    lines.sort(key=lambda l: int(l.split("\t")[0]))
    parts = []
    for line in lines:
        if ":H#" in line:
            continue
        fields = line.split("\t")
        if len(fields) >= 3:
            parts.append(fields[2])
    return "".join(parts).replace(" ", "").rstrip(",")


class TestBashScriptBugs:
    """Evidence that run_rdt.sh produces incorrect output, extracted from
    the bash-produced files stored in reaction_intermediates.zip."""

    def test_dpe12_mol02_empty_species_in_bash_output(self, tmp_path):
        """Bash left DPE12_h MOL_02 (from-side) with an empty species_id.

        MOL_02's InChIKey matches both M_Glc and M_starch1 in
        species_id_inchikey.txt.  The multi-line lookup result caused the
        subsequent grep against from_species_with_cmp to fail on the system
        that produced the zip, resulting in an empty species_id.
        """
        rxn_dir = _extract_rxn_from_zip("DPE12_h", tmp_path)
        mol02_species = (rxn_dir / "MOL_02.species_id").read_text().strip()
        assert mol02_species == ""

    def test_dpe12_mol04_filename_leaked_as_species(self, tmp_path):
        """Bash embedded the literal grep target filename in DPE12_h MOL_04's species_id.

        The species_id file contains 'to_species_with_cmp:M_Glc[h]' where
        'to_species_with_cmp' is the FILENAME of the grep target, not a
        species identifier.  A valid species has the form M_<name>[<cmp>].
        """
        rxn_dir = _extract_rxn_from_zip("DPE12_h", tmp_path)
        mol04_species = (rxn_dir / "MOL_04.species_id").read_text().strip()
        assert "to_species_with_cmp" in mol04_species, (
            f"Expected filename leak, got: {mol04_species!r}"
        )

    def test_dpe12_bash_mapping_has_empty_species_and_filename(self, tmp_path):
        """The bash-produced mapping.txt for DPE12_h contains empty species
        names and a literal filename as species — both self-evidently broken."""
        rxn_dir = _extract_rxn_from_zip("DPE12_h", tmp_path)
        mapping = (rxn_dir / "mapping.txt").read_text()
        assert ":O#1=" in mapping or ":C#1=" in mapping, (
            "Bash mapping should have entries with empty species (e.g. ':C#1=')"
        )
        assert "to_species_with_cmp" in mapping, (
            "Bash mapping should contain the filename 'to_species_with_cmp'"
        )

    def test_dpe12_correct_species_exists_for_mol02(self, tmp_path):
        """MOL_02 SHOULD map to M_starch1[h] — it is present in
        from_species_with_cmp and M_starch1 is a valid InChIKey match."""
        rxn_dir = _extract_rxn_from_zip("DPE12_h", tmp_path)
        inchikey = (rxn_dir / "MOL_02.inchikey").read_text().strip()
        species_table = run_rdt.load_inchikey_table(
            (rxn_dir / "species_id_inchikey.txt").read_text()
        )
        from_species = (rxn_dir / "from_species_with_cmp").read_text().strip().splitlines()

        matches = run_rdt.lookup_species(inchikey, species_table)
        assert "M_starch1" in matches
        assert any("M_starch1" in f for f in from_species)

    def test_dpe12_correct_species_exists_for_mol04(self, tmp_path):
        """MOL_04 SHOULD map to M_Glc[h] — it is present in
        to_species_with_cmp and M_Glc is a valid InChIKey match."""
        rxn_dir = _extract_rxn_from_zip("DPE12_h", tmp_path)
        inchikey = (rxn_dir / "MOL_04.inchikey").read_text().strip()
        species_table = run_rdt.load_inchikey_table(
            (rxn_dir / "species_id_inchikey.txt").read_text()
        )
        to_species = (rxn_dir / "to_species_with_cmp").read_text().strip().splitlines()

        matches = run_rdt.lookup_species(inchikey, species_table)
        assert "M_Glc" in matches
        assert any("M_Glc" in f for f in to_species)

    def test_dpe12_python_produces_valid_species(self, tmp_path):
        """Python correctly identifies M_starch1[h] and M_Glc[h] for DPE12_h,
        producing valid species names without empty entries or filename leaks.

        Every entry in the mapping must start with a species name (M_...),
        never with a bare ':element' or a filename like 'to_species_with_cmp'.
        """
        import re
        rxn_dir = _extract_rxn_from_zip("DPE12_h", tmp_path)
        _strip_generated_files(rxn_dir)
        run_rdt.postprocess_reaction(rxn_dir)
        mapping = (rxn_dir / "mapping.txt").read_text()

        entries = re.split(r'[,]', mapping)
        for entry in entries:
            from_to = entry.split("=")
            for side in from_to:
                if not side:
                    continue
                assert side.startswith("M_"), (
                    f"Expected species to start with 'M_', got: {side!r} "
                    f"in entry {entry!r}"
                )

        assert "to_species_with_cmp" not in mapping
        assert "M_starch1[h]" in mapping
        assert "M_Glc[h]" in mapping

    def test_ornata_mapping_lines_has_substring_false_positive(self, tmp_path):
        """Bash grep 'M_Glu' against to_species_with_cmp matched BOTH
        M_Glu[m] and M_Glu-SeA[m] — a substring false positive.

        The mapping_lines.txt records 'M_Glu[m] M_Glu-SeA[m]' as the
        species for MOL_03, but M_Glu and M_Glu-SeA are different species.
        """
        rxn_dir = _extract_rxn_from_zip("OrnAT_m", tmp_path)
        mapping_lines = (rxn_dir / "mapping_lines.txt").read_text()
        assert "M_Glu[m] M_Glu-SeA[m]" in mapping_lines

    def test_ornata_mol03_inchikey_matches_only_m_glu(self, tmp_path):
        """MOL_03's InChIKey matches ONLY M_Glu (not M_Glu-SeA) in the
        species table — the bash grep substring match was a false positive."""
        rxn_dir = _extract_rxn_from_zip("OrnAT_m", tmp_path)
        inchikey = (rxn_dir / "MOL_03.inchikey").read_text().strip()
        species_table = run_rdt.load_inchikey_table(
            (rxn_dir / "species_id_inchikey.txt").read_text()
        )
        matches = run_rdt.lookup_species(inchikey, species_table)
        assert matches == ["M_Glu"], (
            f"InChIKey {inchikey} should match only M_Glu, got: {matches}"
        )

    def test_ornata_mapping_txt_inconsistent_with_mapping_lines(self, tmp_path):
        """OrnAT_m's mapping.txt CANNOT be derived from its mapping_lines.txt
        via the bash pipeline. The zip contains inconsistent data from mixed
        sources — the mapping_lines.txt has the multi-species grep result,
        but mapping.txt was apparently regenerated separately."""
        rxn_dir = _extract_rxn_from_zip("OrnAT_m", tmp_path)
        mapping_txt = (rxn_dir / "mapping.txt").read_text()
        mapping_lines = (rxn_dir / "mapping_lines.txt").read_text()
        reconstructed = _bash_pipeline_assemble(mapping_lines)
        assert reconstructed != mapping_txt

    def test_ornata_bash_pipeline_produces_concatenated_species(self, tmp_path):
        """When the bash assembly pipeline is applied to OrnAT_m's
        mapping_lines.txt, the space-stripping step concatenates the two
        species names into the invalid string 'M_Glu[m]M_Glu-SeA[m]'."""
        rxn_dir = _extract_rxn_from_zip("OrnAT_m", tmp_path)
        mapping_lines = (rxn_dir / "mapping_lines.txt").read_text()
        reconstructed = _bash_pipeline_assemble(mapping_lines)
        assert "M_Glu[m]M_Glu-SeA[m]" in reconstructed, (
            "Bash pipeline concatenates multi-match species after space removal"
        )

    def test_grep_substring_false_positive(self, tmp_path):
        """Running actual grep 'M_Glu' against a file containing both
        M_Glu[m] and M_Glu-SeA[m] matches BOTH lines.

        This is the root cause of the OrnAT_m bug: grep uses substring
        matching, so M_Glu matches the unrelated species M_Glu-SeA.
        """
        species_file = tmp_path / "to_species_with_cmp"
        species_file.write_text("M_Glu[m]\nM_Glu-SeA[m]\n")
        result = subprocess.run(
            ["grep", "M_Glu", str(species_file)],
            capture_output=True, text=True,
        )
        matches = result.stdout.strip().splitlines()
        assert "M_Glu[m]" in matches
        assert "M_Glu-SeA[m]" in matches, (
            f"grep 'M_Glu' falsely matches M_Glu-SeA[m] — "
            f"all matches: {matches}"
        )

    def test_exact_prefix_match_avoids_false_positive(self, tmp_path):
        """Exact prefix matching (what the script SHOULD do) matches only
        M_Glu[m], not M_Glu-SeA[m]."""
        species_list = ["M_Glu[m]", "M_Glu-SeA[m]"]
        exact_matches = [
            s for s in species_list if s.split("[")[0] == "M_Glu"
        ]
        assert exact_matches == ["M_Glu[m]"]


@pytest.mark.integration
def test_rdt_non_determinism_diagnostic(tmp_path):
    if not RDT_JAR.exists():
        pytest.skip(f"RDT JAR not found at {RDT_JAR}")

    rxn_name = "AspAT_h"
    rxn_dir = _extract_rxn_from_zip(rxn_name, tmp_path)

    smiles = (rxn_dir / "rxn.smiles").read_text().strip()

    NUM_RUNS = 5
    rxn_hashes = set()

    for run_idx in range(NUM_RUNS):
        run_dir = tmp_path / f"{rxn_name}_run{run_idx}"
        run_dir.mkdir()

        try:
            _run_rdt(smiles, RDT_JAR, run_dir)
        except run_rdt.SubprocessError as e:
            print(f"Run {run_idx} failed: {e}", file=sys.stderr)

        rxn_file = run_dir / "ECBLAST_smiles_AAM.rxn"
        if rxn_file.exists():
            rxn_text = rxn_file.read_text()
            mol_blocks = run_rdt.split_rxn_to_mols(rxn_text)
            atom_lines = []
            for mol_block in mol_blocks:
                atoms = run_rdt.parse_mdl_atom_table(mol_block)
                for elem, idx in atoms:
                    atom_lines.append(f"{elem}\t{idx}")
            rxn_hashes.add(
                hashlib.md5("\n".join(atom_lines).encode()).hexdigest()
            )

        shutil.rmtree(run_dir, ignore_errors=True)

    if len(rxn_hashes) > 1:
        print(
            f"WARNING: RDT produced {len(rxn_hashes)} distinct atom mappings "
            f"for {rxn_name} across {NUM_RUNS} runs. This confirms RDT is "
            f"non-deterministic for reactions with symmetric molecules."
        )
    else:
        print(
            f"NOTE: RDT produced identical output for {rxn_name} across "
            f"{NUM_RUNS} runs on this machine. RDT's non-determinism may be "
            f"environment-dependent or intermittent."
        )
