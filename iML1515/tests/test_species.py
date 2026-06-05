import pytest
import preprocess_iml1515 as pp


class TestBiggToBracket:
    def test_simple_cytosol(self):
        base, bracket = pp.bigg_to_bracket("dhap_c")
        assert base == "dhap"
        assert bracket == "dhap[c]"

    def test_extracellular(self):
        base, bracket = pp.bigg_to_bracket("glc__D_e")
        assert base == "glc__D"
        assert bracket == "glc__D[e]"

    def test_periplasm(self):
        base, bracket = pp.bigg_to_bracket("some_met_p")
        assert base == "some_met"
        assert bracket == "some_met[p]"

    def test_double_underscore_preserved(self):
        base, bracket = pp.bigg_to_bracket("cysi__L_c")
        assert base == "cysi__L"
        assert bracket == "cysi__L[c]"

    def test_ala__D_c(self):
        base, bracket = pp.bigg_to_bracket("ala__D_c")
        assert base == "ala__D"
        assert bracket == "ala__D[c]"


class TestSpeciesIdUtils:
    def test_bracket_to_bigg_roundtrip(self):
        assert pp.bigg_to_bracket("dhap_c") == ("dhap", "dhap[c]")

    def test_base_name_used_for_inchikey(self):
        _, bracket = pp.bigg_to_bracket("gtp_c")
        assert bracket == "gtp[c]"
