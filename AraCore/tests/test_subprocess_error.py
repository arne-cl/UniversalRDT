import subprocess
import zipfile
from pathlib import Path

import pytest

import run_rdt

ARACORE_DIR = Path(__file__).resolve().parent.parent
REACTIONS_ZIP = ARACORE_DIR / "reaction_intermediates.zip"
RDT_JAR = Path(
    __import__("os").environ.get(
        "RDT_JAR",
        ARACORE_DIR.parent / "rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar",
    )
)


class TestSubprocessError:
    def test_attributes(self):
        err = run_rdt.SubprocessError(
            cmd=["java", "-jar", "test.jar"],
            returncode=1,
            stdout=b"out",
            stderr=b"err",
        )
        assert err.cmd == ["java", "-jar", "test.jar"]
        assert err.returncode == 1
        assert err.stdout == b"out"
        assert err.stderr == b"err"

    def test_str_contains_cmd_and_returncode(self):
        err = run_rdt.SubprocessError(
            cmd=["java", "-jar", "test.jar"],
            returncode=42,
            stdout=b"",
            stderr=b"bad input",
        )
        s = str(err)
        assert "java -jar test.jar" in s
        assert "42" in s
        assert "bad input" in s

    def test_str_decodes_bytes(self):
        err = run_rdt.SubprocessError(
            cmd=["test"], returncode=1, stdout=b"hello", stderr=b"world",
        )
        s = str(err)
        assert "hello" in s
        assert "world" in s

    def test_str_handles_none(self):
        err = run_rdt.SubprocessError(
            cmd=["test"], returncode=1, stdout=None, stderr=None,
        )
        s = str(err)
        assert "test" in s

    def test_str_handles_string_output(self):
        err = run_rdt.SubprocessError(
            cmd=["test"], returncode=1, stdout="text out", stderr="text err",
        )
        s = str(err)
        assert "text out" in s
        assert "text err" in s


class TestRunRdtJavaCapturesOutput:
    @pytest.mark.integration
    def test_invalid_smiles_raises_subprocess_error(self, tmp_path):
        if not RDT_JAR.exists():
            pytest.skip(f"RDT JAR not found at {RDT_JAR}")
        with pytest.raises(run_rdt.SubprocessError) as exc_info:
            run_rdt.run_rdt_java("NOT_A_VALID_SMILES", RDT_JAR, tmp_path)
        err = exc_info.value
        assert err.returncode != 0
        assert len(err.stdout) > 0 or len(err.stderr) > 0

    @pytest.mark.integration
    def test_missing_jar_raises_subprocess_error(self, tmp_path):
        fake_jar = tmp_path / "nonexistent.jar"
        with pytest.raises(run_rdt.SubprocessError) as exc_info:
            run_rdt.run_rdt_java("C>>C", fake_jar, tmp_path)
        err = exc_info.value
        assert err.returncode != 0
        assert err.stderr is not None


class TestObabelCapturesOutput:
    def test_obabel_to_inchi_empty_input_empty_output(self):
        assert run_rdt.obabel_to_inchi('') == ''

    def test_obabel_to_inchikey_missing_file_returns_none(self):
        assert run_rdt.obabel_to_inchikey('') == ''


class TestProcessReactionErrorHandling:
    @pytest.mark.integration
    def test_process_reaction_returns_false_on_rdt_failure(self, tmp_path, capsys):
        if not RDT_JAR.exists():
            pytest.skip(f"RDT JAR not found at {RDT_JAR}")
        rxn_dir = tmp_path / "bad_reaction"
        rxn_dir.mkdir()
        (rxn_dir / "rxn.smiles").write_text("NOT_VALID_SMILES>>ALSO_BAD")
        success, mapping_lines, mapping_text = run_rdt.process_reaction(rxn_dir, RDT_JAR)
        assert success is False
        assert mapping_lines == ""
        assert mapping_text == ""
        captured = capsys.readouterr()
        assert "Error processing bad_reaction" in captured.err

    def test_process_reaction_returns_false_missing_smiles(self, tmp_path):
        rxn_dir = tmp_path / "no_smiles"
        rxn_dir.mkdir()
        success, mapping_lines, mapping_text = run_rdt.process_reaction(rxn_dir, Path("/fake/jar.jar"))
        assert success is False
        assert mapping_lines == ""
        assert mapping_text == ""

    def test_process_reaction_returns_false_missing_jar(self, tmp_path, capsys):
        rxn_dir = tmp_path / "missing_jar"
        rxn_dir.mkdir()
        (rxn_dir / "rxn.smiles").write_text("C>>C")
        fake_jar = tmp_path / "nonexistent.jar"
        success, mapping_lines, mapping_text = run_rdt.process_reaction(rxn_dir, fake_jar)
        assert success is False
        assert mapping_lines == ""
        assert mapping_text == ""
        captured = capsys.readouterr()
        assert "Error processing missing_jar" in captured.err
