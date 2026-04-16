import shutil
from pathlib import Path

import pytest

import run_rdt

ARAcore_DIR = Path(__file__).resolve().parent.parent
REACTIONS_DIR = ARAcore_DIR / "reaction_intermediates"
GOLDEN_DIR = Path(__file__).resolve().parent / "golden"


@pytest.mark.integration
def test_postprocess_single_reaction(sample_rxn_dir, tmp_path):
    rxn_copy = tmp_path / "FBPA_h"
    shutil.copytree(sample_rxn_dir, rxn_copy)

    for f in rxn_copy.glob("mapping*"):
        f.unlink()
    for f in rxn_copy.glob("MOL_*.mdl"):
        f.unlink()
    for f in rxn_copy.glob("MOL_*.inchi"):
        f.unlink()
    for f in rxn_copy.glob("MOL_*.inchikey"):
        f.unlink()
    for f in rxn_copy.glob("MOL_*.rdt_index"):
        f.unlink()
    for f in rxn_copy.glob("MOL_*.species_id"):
        f.unlink()

    result = run_rdt.postprocess_reaction(rxn_copy)
    assert result is True

    golden = GOLDEN_DIR / "FBPA_h.mapping.txt"
    assert golden.exists()
    expected = golden.read_text()
    actual = (rxn_copy / "mapping.txt").read_text()
    assert actual == expected


@pytest.mark.integration
def test_all_mappings_match_golden_reference(tmp_path):
    failures = []
    rxn_folders = sorted(REACTIONS_DIR.iterdir())
    for rxn_folder in rxn_folders:
        if not rxn_folder.is_dir():
            continue
        golden = GOLDEN_DIR / f"{rxn_folder.name}.mapping.txt"
        if not golden.exists():
            continue
        golden_text = golden.read_text()
        actual_path = rxn_folder / "mapping.txt"
        if not actual_path.exists():
            failures.append(f"{rxn_folder.name}: mapping.txt missing")
            continue
        actual_text = actual_path.read_text()
        if actual_text != golden_text:
            failures.append(f"{rxn_folder.name}: mapping mismatch")
    assert failures == [], (
        f"{len(failures)} mapping mismatches:\n" + "\n".join(failures[:20])
    )


@pytest.mark.integration
def test_postprocess_all_reactions_produce_same_output(tmp_path):
    rxn_folders = sorted(REACTIONS_DIR.iterdir())
    failures = []
    for rxn_folder in rxn_folders:
        if not rxn_folder.is_dir():
            continue
        golden = GOLDEN_DIR / f"{rxn_folder.name}.mapping.txt"
        if not golden.exists():
            continue

        rxn_copy = tmp_path / rxn_folder.name
        shutil.copytree(rxn_folder, rxn_copy)

        for f in rxn_copy.glob("mapping*"):
            f.unlink()
        for f in rxn_copy.glob("MOL_*.mdl"):
            f.unlink()
        for f in rxn_copy.glob("MOL_*.inchi"):
            f.unlink()
        for f in rxn_copy.glob("MOL_*.inchikey"):
            f.unlink()
        for f in rxn_copy.glob("MOL_*.rdt_index"):
            f.unlink()
        for f in rxn_copy.glob("MOL_*.species_id"):
            f.unlink()

        try:
            run_rdt.postprocess_reaction(rxn_copy)
            expected = golden.read_text()
            actual = (rxn_copy / "mapping.txt").read_text()
            if actual != expected:
                failures.append(f"{rxn_folder.name}: mismatch")
        except Exception as e:
            failures.append(f"{rxn_folder.name}: error: {e}")

        shutil.rmtree(rxn_copy)

    assert failures == [], (
        f"{len(failures)} failures:\n" + "\n".join(failures[:20])
    )
