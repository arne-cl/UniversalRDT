import json
import pytest
import preprocess_iml1515 as pp


class TestLoadModel:
    def test_loads_json_gz(self, tmp_model_path):
        model = pp.load_model(tmp_model_path)
        assert model["id"] == "test_model"
        assert len(model["metabolites"]) == 8
        assert len(model["reactions"]) == 6

    def test_nonexistent_path_raises(self):
        with pytest.raises(FileNotFoundError):
            pp.load_model("/nonexistent/path.json.gz")


class TestCache:
    def test_empty_cache(self, tmp_path):
        path = tmp_path / "cache.json"
        cache = pp.load_cache(path)
        assert cache == {}

    def test_nonexistent_cache_returns_empty(self, tmp_path):
        cache = pp.load_cache(tmp_path / "nonexistent.json")
        assert cache == {}

    def test_save_and_load(self, tmp_path):
        path = tmp_path / "cache.json"
        data = {"cpd00001": {"smiles": "O", "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N"}}
        pp.save_cache(data, path)
        loaded = pp.load_cache(path)
        assert loaded == data


class TestCLI:
    def test_defaults(self):
        args = pp.parse_args([])
        assert args.model_path is not None
        assert args.output_dir is not None
        assert args.cache_file is not None
        assert args.verbose == 0

    def test_custom_paths(self):
        args = pp.parse_args([
            "--model-path", "/custom/model.json.gz",
            "--output-dir", "/out",
            "--cache-file", "/cache.json",
            "-v",
        ])
        assert str(args.model_path) == "/custom/model.json.gz"
        assert str(args.output_dir) == "/out"
        assert str(args.cache_file) == "/cache.json"
        assert args.verbose == 1

    def test_verbose_repeatable(self):
        args = pp.parse_args(["-vv"])
        assert args.verbose == 2
