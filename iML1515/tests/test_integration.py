import pytest

import preprocess_iml1515 as pp


class TestFullPipeline:
    """Integration tests running the full preprocessing pipeline."""

    def test_process_all_reactions(self, small_model, small_model_smiles_map, tmp_output_dir):
        model = small_model
        stats = pp.process_model(model, small_model_smiles_map, tmp_output_dir)
        assert stats["total"] == 6
        assert stats["excluded_exchange"] == 2
        assert stats["excluded_biomass"] == 1
        assert stats["excluded_no_smiles"] == 1
        assert stats["included"] == 2

    def test_output_dirs_created(self, small_model, small_model_smiles_map, tmp_output_dir):
        pp.process_model(model=small_model, smiles_map=small_model_smiles_map, output_dir=tmp_output_dir)
        assert (tmp_output_dir / "PPA").is_dir()
        assert (tmp_output_dir / "OK_RXN").is_dir()
        assert not (tmp_output_dir / "EX_glc_e").exists()
        assert not (tmp_output_dir / "BIOMASS_Ec_iML1515_core_75p37M").exists()
        assert not (tmp_output_dir / "R_STUCK").exists()

    def test_excluded_count(self, small_model, tmp_model_path, small_model_smiles_map, tmp_output_dir):
        stats = pp.process_model(small_model, small_model_smiles_map, tmp_output_dir)
        assert stats["total"] == 6
        assert stats["included"] == 2
        assert stats["excluded_exchange"] == 2
        assert stats["excluded_biomass"] == 1
        assert stats["excluded_no_smiles"] == 1
