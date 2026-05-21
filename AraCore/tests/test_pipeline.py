import hashlib
import os
import shutil
import subprocess
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
    subprocess.run(cmd, cwd=str(cwd),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


@pytest.mark.integration
def test_postprocess_single_reaction(tmp_path):
    rxn_dir = _extract_rxn_from_zip("FBPA_h", tmp_path)
    _strip_generated_files(rxn_dir)

    result = run_rdt.postprocess_reaction(rxn_dir)
    assert result is True

    golden = GOLDEN_DIR / "FBPA_h.mapping.txt"
    assert golden.exists()
    expected = golden.read_text()
    actual = (rxn_dir / "mapping.txt").read_text()
    assert actual == expected


@pytest.mark.integration
def test_python_matches_bash_for_all_reactions(tmp_path):
    failures = []
    for rxn_name in _all_rxn_names_from_zip():
        rxn_dir = _extract_rxn_from_zip(rxn_name, tmp_path / rxn_name)
        bash_mapping = (rxn_dir / "mapping.txt").read_text()
        _strip_generated_files(rxn_dir)

        run_rdt.postprocess_reaction(rxn_dir)
        py_mapping = (rxn_dir / "mapping.txt").read_text()

        if py_mapping != bash_mapping:
            failures.append(rxn_name)

    unexpected = [f for f in failures if f not in KNOWN_PYTHON_BASH_DIFFS]
    assert unexpected == [], (
        f"{len(unexpected)} unexpected Python-vs-bash mismatches:\n"
        + "\n".join(unexpected[:20])
    )
    assert set(failures) == KNOWN_PYTHON_BASH_DIFFS, (
        f"Expected known diffs {sorted(KNOWN_PYTHON_BASH_DIFFS)}, "
        f"got {sorted(failures)}"
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
def test_rdt_non_determinism_diagnostic(tmp_path):
    if not RDT_JAR.exists():
        pytest.skip(f"RDT JAR not found at {RDT_JAR}")

    rxn_name = "AspAT_h"
    rxn_dir = REACTIONS_DIR / rxn_name
    if not rxn_dir.is_dir():
        pytest.skip(f"{rxn_name} not in reaction_intermediates/ (run make prepare-aracore)")

    smiles = (rxn_dir / "rxn.smiles").read_text().strip()

    NUM_RUNS = 5
    rxn_hashes = set()

    for run_idx in range(NUM_RUNS):
        run_dir = tmp_path / f"{rxn_name}_run{run_idx}"
        run_dir.mkdir()

        _run_rdt(smiles, RDT_JAR, run_dir)

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
