import pytest
import preprocess_iml1515 as pp


class TestSplitBySign:
    def test_ppa_reaction(self, small_model):
        rxn = small_model["reactions"][0]
        mets = rxn["metabolites"]
        reactants, products = pp.split_by_sign(mets)
        assert set(reactants.keys()) == {"h2o_c", "ppi_c"}
        assert set(products.keys()) == {"h_c", "pi_c"}

    def test_ok_rxn(self, small_model):
        rxn = small_model["reactions"][4]
        mets = rxn["metabolites"]
        reactants, products = pp.split_by_sign(mets)
        assert set(reactants.keys()) == {"dhap_c", "gtp_c"}
        assert set(products.keys()) == {"h2o_c", "h_c"}


class TestFlattenCoefficients:
    def test_unit_coefficients(self):
        mets = {"dhap_c": -1.0, "gtp_c": -1.0}
        result = pp.flatten_coefficients(mets)
        assert result == ["dhap_c", "gtp_c"]

    def test_non_unit_coefficients(self):
        mets = {"h2o_c": -1.0, "pi_c": 2.0, "h_c": 1.0, "ppi_c": -1.0}
        result = pp.flatten_coefficients(mets)
        assert result == ["h2o_c", "pi_c", "pi_c", "h_c", "ppi_c"]

    def test_large_coefficient(self):
        mets = {"h_c": 75.37723}
        result = pp.flatten_coefficients(mets)
        assert len(result) == 75
        assert all(x == "h_c" for x in result)

    def test_zero_coefficient_skipped(self):
        mets = {"a_c": 0.0, "b_c": -1.0}
        result = pp.flatten_coefficients(mets)
        assert result == ["b_c"]

    def test_negative_retains_sign_for_reactants(self):
        mets = {"x_c": -2.0}
        result = pp.flatten_coefficients(mets)
        assert result == ["x_c", "x_c"]


class TestBuildSmilesString:
    def test_basic(self):
        result = pp.build_smiles_string(["O"], ["[H+]"])
        assert result == "O>>[H+]"

    def test_multiple_reactants_products(self):
        result = pp.build_smiles_string(
            ["O", "CC(=O)[O-]"],
            ["[H+]", "OP(=O)([O-])[O-]"],
        )
        assert result == "O.CC(=O)[O-]>>[H+].OP(=O)([O-])[O-]"

    def test_no_reactants(self):
        result = pp.build_smiles_string([], ["C"])
        assert result == ">>C"

    def test_no_products(self):
        result = pp.build_smiles_string(["C"], [])
        assert result == "C>>"

    def test_flattened_ppi(self):
        result = pp.build_smiles_string(
            ["O", "OP(=O)([O-])OP(=O)([O-])[O-]", "OP(=O)([O-])OP(=O)([O-])[O-]"],
            ["[H+]", "OP(=O)([O-])[O-]", "OP(=O)([O-])[O-]"],
        )
        parts = result.split(">>")
        left = parts[0].split(".")
        right = parts[1].split(".")
        assert len(left) == 3
        assert len(right) == 3


class TestBuildSmilesFromReaction:
    def test_ppa(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][0]
        result = pp.build_smiles_from_reaction(rxn, small_model_smiles_map)
        assert ">>" in result
        left, right = result.split(">>")
        left_mols = left.split(".")
        right_mols = right.split(".")
        assert len(left_mols) == 2
        assert len(right_mols) == 3


class TestGetFromToSpecies:
    def test_ppa(self, small_model):
        rxn = small_model["reactions"][0]
        from_species, to_species = pp.get_from_to_species(rxn)
        assert sorted(from_species) == sorted(["h2o[c]", "ppi[c]"])
        assert sorted(to_species) == sorted(["h[c]", "pi[c]", "pi[c]"])

    def test_ppa_deduplicated(self, small_model):
        from_species, to_species = pp.get_from_to_species(
            small_model["reactions"][0]
        )
        from_unique = sorted(set(from_species))
        to_unique = sorted(set(to_species))
        assert from_unique == ["h2o[c]", "ppi[c]"]
        assert to_unique == ["h[c]", "pi[c]"]


class TestBuildInchiKeyMap:
    def test_ppa(self, small_model, small_model_smiles_map):
        rxn = small_model["reactions"][0]
        result = pp.build_inchikey_map(rxn, small_model_smiles_map)
        assert result["h2o"] == "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        assert result["h"] == "GPRLSGONYQIRFK-UHFFFAOYSA-N"
        assert result["pi"] == "NBIIXXVUZAFLBC-UHFFFAOYSA-K"
        assert result["ppi"] == "MNUQJQIHOXMAMT-UHFFFAOYSA-J"

    def test_base_name_no_compartment(self, small_model, small_model_smiles_map):
        result = pp.build_inchikey_map(
            small_model["reactions"][0], small_model_smiles_map
        )
        for key in result:
            assert "[" not in key
            assert "_" not in key or key.count("_") <= 1


class TestWriteReactionOutput:
    def test_creates_four_files(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][0]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        rxn_dir = tmp_output_dir / "PPA"
        assert rxn_dir.is_dir()
        assert (rxn_dir / "rxn.smiles").is_file()
        assert (rxn_dir / "from_species_with_cmp").is_file()
        assert (rxn_dir / "to_species_with_cmp").is_file()
        assert (rxn_dir / "species_id_inchikey.txt").is_file()

    def test_rxn_smiles_content(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][4]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        content = (tmp_output_dir / "OK_RXN" / "rxn.smiles").read_text().strip()
        assert ">>" in content

    def test_inchikey_file_tsv_format(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][0]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        lines = (tmp_output_dir / "PPA" / "species_id_inchikey.txt").read_text().strip().split("\n")
        for line in lines:
            parts = line.split("\t")
            assert len(parts) == 2
            assert parts[0]
            assert parts[1]

    def test_from_species_content(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][0]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        lines = (tmp_output_dir / "PPA" / "from_species_with_cmp").read_text().strip().split("\n")
        assert "h2o[c]" in lines
        assert "ppi[c]" in lines
        assert "h[c]" not in lines
        assert "pi[c]" not in lines

    def test_to_species_content(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][0]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        lines = (tmp_output_dir / "PPA" / "to_species_with_cmp").read_text().strip().split("\n")
        assert "h[c]" in lines
        assert "pi[c]" in lines
        assert "h2o[c]" not in lines
        assert "ppi[c]" not in lines

    def test_skipped_reaction_not_written(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][1]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        assert not (tmp_output_dir / "EX_glc_e").exists()


class TestWriteReactionOutputIntegration:
    def test_ppa_all_files_match(self, small_model, small_model_smiles_map, tmp_output_dir):
        rxn = small_model["reactions"][0]
        pp.write_reaction_output(rxn, small_model_smiles_map, tmp_output_dir)
        rxn_dir = tmp_output_dir / "PPA"

        smiles = (rxn_dir / "rxn.smiles").read_text().strip()
        assert smiles.count(">>") == 1
        left, right = smiles.split(">>")
        assert left
        assert right

        from_lines = (rxn_dir / "from_species_with_cmp").read_text().strip().split("\n")
        to_lines = (rxn_dir / "to_species_with_cmp").read_text().strip().split("\n")
        ik_lines = (rxn_dir / "species_id_inchikey.txt").read_text().strip().split("\n")

        assert sorted(from_lines) == ["h2o[c]", "ppi[c]"]
        assert sorted(to_lines) == ["h[c]", "pi[c]"]
        assert len(ik_lines) == 4
