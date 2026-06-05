import pytest
import preprocess_iml1515 as pp


class TestIsExchange:
    def test_ex_prefix(self, small_model):
        rxn = small_model["reactions"][1]
        assert rxn["id"] == "EX_glc_e"
        assert pp.is_exchange_reaction(rxn) is True

    def test_single_met_boundary(self, small_model):
        rxn = small_model["reactions"][5]
        assert rxn["id"] == "SINGLE_MET_BOUNDARY"
        assert pp.is_exchange_reaction(rxn) is True

    def test_normal_reaction_not_exchange(self, small_model):
        rxn = small_model["reactions"][0]
        assert rxn["id"] == "PPA"
        assert pp.is_exchange_reaction(rxn) is False

    def test_ok_rxn_not_exchange(self, small_model):
        rxn = small_model["reactions"][4]
        assert rxn["id"] == "OK_RXN"
        assert pp.is_exchange_reaction(rxn) is False


class TestIsBiomass:
    def test_biomass_detected(self, small_model):
        rxn = small_model["reactions"][2]
        assert rxn["id"] == "BIOMASS_Ec_iML1515_core_75p37M"
        assert pp.is_biomass_reaction(rxn) is True

    def test_non_biomass(self, small_model):
        rxn = small_model["reactions"][0]
        assert pp.is_biomass_reaction(rxn) is False

    def test_exchange_not_biomass(self, small_model):
        rxn = small_model["reactions"][1]
        assert pp.is_biomass_reaction(rxn) is False


class TestHasAllSmiles:
    def test_all_smiles_present(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][0]
        assert pp.has_all_smiles(rxn, small_model_smiles_map) is True

    def test_missing_smiles_detected(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][3]
        assert rxn["id"] == "R_STUCK"
        assert pp.has_all_smiles(rxn, small_model_smiles_map) is False


class TestFilterReactions:
    def test_ppa_included(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][0]
        assert pp.filter_reaction(rxn, small_model_smiles_map) is True

    def test_exchange_excluded(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][1]
        assert pp.filter_reaction(rxn, small_model_smiles_map) is False

    def test_biomass_excluded(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][2]
        assert pp.filter_reaction(rxn, small_model_smiles_map) is False

    def test_no_smiles_excluded(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][3]
        assert pp.filter_reaction(rxn, small_model_smiles_map) is False

    def test_single_met_boundary_excluded(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][5]
        assert pp.filter_reaction(rxn, small_model_smiles_map) is False
