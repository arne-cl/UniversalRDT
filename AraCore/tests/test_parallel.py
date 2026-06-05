"""Tests for parallel reaction processing in main()."""
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

import run_rdt


@pytest.fixture
def reaction_dirs(tmp_path):
    """Create 6 fake reaction dirs with rxn.smiles files."""
    dirs = []
    for i in range(6):
        d = tmp_path / f"RXN_{i:03d}"
        d.mkdir()
        (d / "rxn.smiles").write_text(f"C>>O  # reaction {i}")
        dirs.append(d)
    return dirs


class TestRunWorker:
    """_run_worker dispatches to the correct function based on mode."""

    def test_full_pipeline_calls_process_reaction(self, reaction_dirs):
        jar = Path("/fake.jar")
        with patch("run_rdt.process_reaction", return_value=(True, "a", "b")) as mock:
            result = run_rdt._run_worker(reaction_dirs[0], jar, False)
            mock.assert_called_once_with(reaction_dirs[0], jar)
            assert result == (True, "a", "b")

    def test_postprocess_only_calls_postprocess(self, reaction_dirs):
        jar = Path("/fake.jar")
        with patch("run_rdt.postprocess_reaction", return_value=(False, "", "")) as mock:
            result = run_rdt._run_worker(reaction_dirs[0], jar, True)
            mock.assert_called_once_with(reaction_dirs[0])
            assert result == (False, "", "")


class TestProcessReactionsParallel:
    """process_reactions_parallel dispatches work to a process pool."""

    def test_uses_process_pool_executor_with_workers(self, reaction_dirs):
        with patch("run_rdt._run_worker", return_value=(True, "", "")):
            with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
                executor = MagicMock()
                mock_pool_cls.return_value.__enter__ = MagicMock(return_value=executor)
                mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
                run_rdt.process_reactions_parallel(
                    reaction_dirs, rdt_jar=Path("/fake.jar"), workers=4,
                )
                mock_pool_cls.assert_called_once_with(max_workers=4)

    def test_submits_all_folders(self, reaction_dirs):
        with patch("run_rdt._run_worker", return_value=(True, "", "")):
            with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
                executor = MagicMock()
                mock_pool_cls.return_value.__enter__ = MagicMock(return_value=executor)
                mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
                run_rdt.process_reactions_parallel(
                    reaction_dirs, rdt_jar=Path("/fake.jar"), workers=4,
                )
                assert executor.submit.call_count == len(reaction_dirs)

    def test_submits_run_worker_with_correct_args(self, reaction_dirs):
        jar = Path("/path/to/rdt.jar")
        with patch("run_rdt._run_worker", return_value=(True, "", "")):
            with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
                executor = MagicMock()
                mock_pool_cls.return_value.__enter__ = MagicMock(return_value=executor)
                mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
                run_rdt.process_reactions_parallel(
                    reaction_dirs, rdt_jar=jar, workers=4,
                )
                for i, c in enumerate(executor.submit.call_args_list):
                    assert c[0][0] == run_rdt._run_worker
                    assert c[0][1] == reaction_dirs[i]
                    assert c[0][2] == jar
                    assert c[0][3] is False

    def test_counts_successes_from_futures(self, reaction_dirs):
        results = [(True, "", ""), (True, "", ""), (False, "", ""),
                   (True, "", ""), (False, "", ""), (True, "", "")]
        with patch("run_rdt._run_worker", side_effect=results):
            with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
                mock_pool_cls.return_value.__enter__ = MagicMock()
                mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
                mock_pool_cls.return_value.__enter__.return_value.__iter__ = MagicMock(
                    return_value=iter([])
                )

                def fake_submit(fn, *args, **kwargs):
                    future = MagicMock()
                    result = fn(*args, **kwargs)
                    future.result.return_value = result
                    return future

                mock_pool_cls.return_value.__enter__.return_value.submit = fake_submit
                success, total = run_rdt.process_reactions_parallel(
                    reaction_dirs, rdt_jar=Path("/fake.jar"), workers=2,
                )
                assert success == 4
                assert total == 6

    def test_workers_default_is_none(self, reaction_dirs):
        with patch("run_rdt._run_worker", return_value=(True, "", "")):
            with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
                executor = MagicMock()
                mock_pool_cls.return_value.__enter__ = MagicMock(return_value=executor)
                mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
                run_rdt.process_reactions_parallel(
                    reaction_dirs, rdt_jar=Path("/fake.jar"),
                )
                mock_pool_cls.assert_called_once_with(max_workers=None)

    def test_empty_input(self):
        with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
            executor = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=executor)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
            success, total = run_rdt.process_reactions_parallel(
                [], rdt_jar=Path("/fake.jar"), workers=4,
            )
            assert success == 0
            assert total == 0

    def test_postprocess_only_passes_true_to_worker(self, reaction_dirs):
        with patch("run_rdt._run_worker", return_value=(True, "", "")):
            with patch("run_rdt.ProcessPoolExecutor") as mock_pool_cls:
                executor = MagicMock()
                mock_pool_cls.return_value.__enter__ = MagicMock(return_value=executor)
                mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
                run_rdt.process_reactions_parallel(
                    reaction_dirs, rdt_jar=Path("/fake.jar"),
                    postprocess_only=True,
                )
                for c in executor.submit.call_args_list:
                    assert c[0][3] is True


class TestWorkersCliArg:
    """main() accepts --workers and passes it through."""

    def test_default_workers_is_none(self, tmp_path):
        fake_dir = tmp_path / "rxns"
        fake_dir.mkdir()
        with patch("run_rdt.process_reactions_parallel", return_value=(0, 0)) as mock:
            with patch("run_rdt._resolve_reactions_dir", return_value=fake_dir):
                run_rdt.main(["--reactions-dir", str(fake_dir)])
                _, kwargs = mock.call_args
                assert kwargs["workers"] is None

    def test_workers_flag_passed_through(self, tmp_path):
        fake_dir = tmp_path / "rxns"
        fake_dir.mkdir()
        with patch("run_rdt.process_reactions_parallel", return_value=(0, 0)) as mock:
            with patch("run_rdt._resolve_reactions_dir", return_value=fake_dir):
                run_rdt.main(["--reactions-dir", str(fake_dir), "--workers", "8"])
                _, kwargs = mock.call_args
                assert kwargs["workers"] == 8
