import shutil
import zipfile
from pathlib import Path

import pytest

ARACORE_DIR = Path(__file__).resolve().parent.parent
REACTIONS_ZIP = ARACORE_DIR / "reaction_intermediates.zip"
REACTIONS_DIR = ARACORE_DIR / "reaction_intermediates"
GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

SAMPLE_RXN = "FBPA_h"


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


@pytest.fixture
def aracore_dir():
    return ARACORE_DIR


@pytest.fixture
def reactions_dir():
    return REACTIONS_DIR


@pytest.fixture
def golden_dir():
    return GOLDEN_DIR


@pytest.fixture
def sample_rxn_dir(tmp_path):
    return _extract_rxn_from_zip(SAMPLE_RXN, tmp_path)


@pytest.fixture
def sample_rxn_text(sample_rxn_dir):
    return (sample_rxn_dir / "ECBLAST_smiles_AAM.rxn").read_text()


@pytest.fixture
def sample_mol_01_text(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_01").read_text()


@pytest.fixture
def sample_mol_01_mdl(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_01.mdl").read_text()


@pytest.fixture
def sample_mol_01_inchi(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_01.inchi").read_text()


@pytest.fixture
def sample_mol_01_inchikey(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_01.inchikey").read_text()


@pytest.fixture
def sample_mol_01_rdt_index(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_01.rdt_index").read_text()


@pytest.fixture
def sample_mol_02_inchi(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_02.inchi").read_text()


@pytest.fixture
def sample_mol_02_rdt_index(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_02.rdt_index").read_text()


@pytest.fixture
def sample_mol_03_inchi(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_03.inchi").read_text()


@pytest.fixture
def sample_mol_03_rdt_index(sample_rxn_dir):
    return (sample_rxn_dir / "MOL_03.rdt_index").read_text()


@pytest.fixture
def sample_species_inchikey(sample_rxn_dir):
    return (sample_rxn_dir / "species_id_inchikey.txt").read_text()


@pytest.fixture
def sample_from_species(sample_rxn_dir):
    return (sample_rxn_dir / "from_species_with_cmp").read_text()


@pytest.fixture
def sample_to_species(sample_rxn_dir):
    return (sample_rxn_dir / "to_species_with_cmp").read_text()


@pytest.fixture
def sample_mapping_lines(sample_rxn_dir):
    return (sample_rxn_dir / "mapping_lines.txt").read_text()


@pytest.fixture
def sample_mapping_txt(sample_rxn_dir):
    return (sample_rxn_dir / "mapping.txt").read_text()
