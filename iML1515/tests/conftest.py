import sys
import json
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

IML1515_DIR = Path(__file__).resolve().parent.parent


@pytest.fixture
def small_model():
    """Minimal model with a few reactions for unit testing."""
    return {
        "metabolites": [
            {
                "id": "dhap_c",
                "name": "Dihydroxyacetone phosphate",
                "compartment": "c",
                "formula": "C3H5O6P",
                "annotation": {
                    "inchi_key": ["GNGACRATGGDKBX-UHFFFAOYSA-L"],
                    "seed.compound": ["cpd00095"],
                    "chebi": ["CHEBI:16108"],
                },
            },
            {
                "id": "gtp_c",
                "name": "GTP",
                "compartment": "c",
                "formula": "C10H16N5O14P3",
                "annotation": {
                    "inchi_key": ["ZKHQWZAMYRWXGA-KQYNXXCUSA-K"],
                    "seed.compound": ["cpd00037"],
                    "chebi": ["CHEBI:15996"],
                },
            },
            {
                "id": "h2o_c",
                "name": "H2O",
                "compartment": "c",
                "formula": "H2O",
                "annotation": {
                    "inchi_key": ["XLYOFNOQVPJJNP-UHFFFAOYSA-N"],
                    "seed.compound": ["cpd00001"],
                    "chebi": ["CHEBI:15377"],
                },
            },
            {
                "id": "h_c",
                "name": "H+",
                "compartment": "c",
                "formula": "H",
                "annotation": {
                    "inchi_key": ["GPRLSGONYQIRFK-UHFFFAOYSA-N"],
                    "seed.compound": ["cpd00067"],
                    "chebi": ["CHEBI:15378"],
                },
            },
            {
                "id": "pi_c",
                "name": "Phosphate",
                "compartment": "c",
                "formula": "HO4P",
                "annotation": {
                    "inchi_key": ["NBIIXXVUZAFLBC-UHFFFAOYSA-K"],
                    "seed.compound": ["cpd00009"],
                    "chebi": ["CHEBI:18367"],
                },
            },
            {
                "id": "ppi_c",
                "name": "Diphosphate",
                "compartment": "c",
                "formula": "H2P2O7",
                "annotation": {
                    "inchi_key": ["MNUQJQIHOXMAMT-UHFFFAOYSA-J"],
                    "seed.compound": ["cpd00008"],
                    "chebi": ["CHEBI:18307"],
                },
            },
            {
                "id": "glc__D_e",
                "name": "D-Glucose",
                "compartment": "e",
                "formula": "C6H12O6",
                "annotation": {
                    "inchi_key": ["WQZGKKKJIJFFOK-PHYPRBDUSA-N"],
                    "seed.compound": ["cpd00027"],
                    "chebi": ["CHEBI:17634"],
                },
            },
            {
                "id": "stuck_met_c",
                "name": "Metabolite without SMILES",
                "compartment": "c",
                "formula": "X",
                "annotation": {},
            },
        ],
        "reactions": [
            {
                "id": "PPA",
                "name": "Pyrophosphatase",
                "metabolites": {
                    "h2o_c": -1.0,
                    "h_c": 1.0,
                    "pi_c": 2.0,
                    "ppi_c": -1.0,
                },
                "lower_bound": -1000.0,
                "upper_bound": 1000.0,
                "subsystem": "test",
            },
            {
                "id": "EX_glc_e",
                "name": "D-Glucose exchange",
                "metabolites": {"glc__D_e": -1.0},
                "lower_bound": 0.0,
                "upper_bound": 1000.0,
                "subsystem": "Exchange",
            },
            {
                "id": "BIOMASS_Ec_iML1515_core_75p37M",
                "name": "Biomass",
                "metabolites": {"dhap_c": -1.0, "h_c": 70.0},
                "lower_bound": 0.0,
                "upper_bound": 1000.0,
                "subsystem": "Biomass",
            },
            {
                "id": "R_STUCK",
                "name": "Reaction with missing SMILES",
                "metabolites": {"stuck_met_c": -1.0, "dhap_c": 1.0},
                "lower_bound": 0.0,
                "upper_bound": 1000.0,
                "subsystem": "test",
            },
            {
                "id": "OK_RXN",
                "name": "Normal reaction",
                "metabolites": {"dhap_c": -1.0, "gtp_c": -1.0, "h2o_c": 1.0, "h_c": 1.0},
                "lower_bound": -1000.0,
                "upper_bound": 1000.0,
                "subsystem": "test",
            },
            {
                "id": "SINGLE_MET_BOUNDARY",
                "name": "Single metabolite boundary",
                "metabolites": {"dhap_c": 1.0},
                "lower_bound": 0.0,
                "upper_bound": 1000.0,
                "subsystem": "test",
            },
        ],
        "id": "test_model",
        "compartments": {"c": "cytosol", "e": "extracellular"},
        "version": 1,
    }


@pytest.fixture
def tmp_model_path(small_model, tmp_path):
    """Write the small model as a gzipped JSON file."""
    import gzip

    path = tmp_path / "model.json.gz"
    with gzip.open(path, "wt") as f:
        json.dump(small_model, f)
    return path


@pytest.fixture
def small_model_smiles_map():
    return {
        "dhap_c": ("O=C[C@H](O)COP(=O)([O-])[O-]", "GNGACRATGGDKBX-UHFFFAOYSA-L"),
        "gtp_c": ("C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N=C(NC2=O)N", "ZKHQWZAMYRWXGA-KQYNXXCUSA-K"),
        "h2o_c": ("O", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"),
        "h_c": ("[H+]", "GPRLSGONYQIRFK-UHFFFAOYSA-N"),
        "pi_c": ("OP(=O)([O-])[O-]", "NBIIXXVUZAFLBC-UHFFFAOYSA-K"),
        "ppi_c": ("OP(=O)([O-])OP(=O)([O-])[O-]", "MNUQJQIHOXMAMT-UHFFFAOYSA-J"),
        "glc__D_e": ("C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", "WQZGKKKJIJFFOK-PHYPRBDUSA-N"),
    }


@pytest.fixture
def tmp_output_dir(tmp_path):
    return tmp_path / "reaction_intermediates"
