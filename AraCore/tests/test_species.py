import run_rdt


def test_load_inchikey_table(sample_species_inchikey):
    table = run_rdt.load_inchikey_table(sample_species_inchikey)
    table_dict = {k: v for k, v in table}
    assert table_dict["LXJXRIRHZLFYRP-VKHMYHEASA-L"] == "M_GAP"
    assert table_dict["GNGACRATGGDKBX-UHFFFAOYSA-L"] == "M_DHAP"
    assert table_dict["RNBGYGVWRKECFJ-VRPWFDPXSA-J"] == "M_FBP"


def test_lookup_species_exact(sample_species_inchikey):
    table = run_rdt.load_inchikey_table(sample_species_inchikey)
    result = run_rdt.lookup_species("LXJXRIRHZLFYRP-VKHMYHEASA-L", table)
    assert result == ["M_GAP"]


def test_lookup_species_prefix_fallback(sample_species_inchikey):
    table = run_rdt.load_inchikey_table(sample_species_inchikey)
    result = run_rdt.lookup_species("LXJXRIRHZLFYRP-XXXXXXXXXX-N", table)
    assert result == ["M_GAP"]


def test_lookup_species_miss(sample_species_inchikey):
    table = run_rdt.load_inchikey_table(sample_species_inchikey)
    result = run_rdt.lookup_species("AAAAAAAAAAAAA-XXXXXXXXXX-N", table)
    assert result == []


def test_lookup_species_multiple_prefix_matches():
    table = [
        ("GMKMEZVLHJARHF-WHFBIAKZSA-N", "M_DAP"),
        ("GMKMEZVLHJARHF-SYDPRGILSA-N", "M_mDAP"),
    ]
    result = run_rdt.lookup_species("GMKMEZVLHJARHF-UHFFFAOYSA-N", table)
    assert result == ["M_DAP", "M_mDAP"]


def test_load_species_list(sample_from_species):
    species = run_rdt.load_species_list(sample_from_species)
    assert species == ["M_DHAP[h]", "M_GAP[h]"]


def test_find_species_with_cmp():
    result = run_rdt.find_species_with_cmp("M_GAP", ["M_GAP[h]", "M_DHAP[h]"])
    assert result == "M_GAP[h]"


def test_find_species_with_cmp_multiple():
    result = run_rdt.find_species_with_cmp("M_Pi", ["M_ADP[c]", "M_Pi[c]", "M_Pi[h]"])
    assert result == "M_Pi[c] M_Pi[h]"


def test_find_species_with_cmp_miss():
    result = run_rdt.find_species_with_cmp("M_XXX", ["M_GAP[h]", "M_DHAP[h]"])
    assert result is None


def test_find_species_with_cmp_multi_substring_false_positive():
    """find_species_with_cmp_multi MUST NOT use substring matching.

    When ``species_ids`` contains ``"M_Glu"``, it should match only
    ``"M_Glu[m]"``, not ``"M_Glu-SeA[m]"`` — these are different
    metabolites and substring matching falsely conflates them.
    """
    result = run_rdt.find_species_with_cmp_multi(
        species_ids=["M_Glu"],
        cmp_list=["M_Glu[m]", "M_Glu-SeA[m]"],
    )
    assert result == "M_Glu[m]"


def test_find_species_with_cmp_multi_multiple_candidates_different_ichikey():
    """When two candidates share the same InChIKey (e.g. M_Glc and M_starch1),
    both are searched against the compartmented list.  Only the one with a
    matching compartment entry is returned.

    This is the DPE reaction case: M_Glc and M_starch1 have identical
    InChIKeys because starch is a glucose polymer.  The compartment filter
    (from_species_with_cmp vs to_species_with_cmp) resolves the ambiguity.
    """
    result = run_rdt.find_species_with_cmp_multi(
        species_ids=["M_Glc", "M_starch1"],
        cmp_list=["M_starch1[h]"],
    )
    assert result == "M_starch1[h]"

    result = run_rdt.find_species_with_cmp_multi(
        species_ids=["M_Glc", "M_starch1"],
        cmp_list=["M_Glc[h]"],
    )
    assert result == "M_Glc[h]"
