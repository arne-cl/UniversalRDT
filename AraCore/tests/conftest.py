from pathlib import Path

import pytest

ARAcore_DIR = Path(__file__).resolve().parent.parent
REACTIONS_DIR = ARAcore_DIR / "reaction_intermediates"
GOLDEN_DIR = Path(__file__).resolve().parent / "golden"

SAMPLE_RXN = "FBPA_h"


@pytest.fixture
def aracore_dir():
    return ARAcore_DIR


@pytest.fixture
def reactions_dir():
    return REACTIONS_DIR


@pytest.fixture
def golden_dir():
    return GOLDEN_DIR


@pytest.fixture
def sample_rxn_dir():
    return REACTIONS_DIR / SAMPLE_RXN


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
