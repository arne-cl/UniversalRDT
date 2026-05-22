import pytest

import unite_mappings

GOLDEN_UNITE_DIR = (
    __import__("pathlib").Path(__file__).resolve().parent / "golden_unite"
)
GOLDEN_DIR = __import__("pathlib").Path(__file__).resolve().parent / "golden"
REACTIONS_DIR = (
    __import__("pathlib").Path(__file__).resolve().parent.parent
    / "reaction_intermediates.zip"
)


# ── Unit tests ───────────────────────────────────────────────────────


class TestSplitAndSortMapping:
    def test_basic_split(self):
        text = "B=C,A=B,C=A"
        result = unite_mappings.split_and_sort_mapping(text)
        assert result == ["A=B", "B=C", "C=A"]

    def test_single_pair(self):
        assert unite_mappings.split_and_sort_mapping("A=B") == ["A=B"]

    def test_empty(self):
        assert unite_mappings.split_and_sort_mapping("") == []
        assert unite_mappings.split_and_sort_mapping("  \n") == []

    def test_preserves_dashes(self):
        text = "M_A-CoA[h]:N#1=M_CoA[h]:N#1,M_A-CoA[h]:N#2=M_CoA[h]:N#2"
        result = unite_mappings.split_and_sort_mapping(text)
        assert len(result) == 2
        assert all("M_A-CoA" in pair for pair in result)


class TestFilterNitrogenMappings:
    def test_filters_nitrogen(self):
        lines = [
            "rxn1 M_A[h]:C#1=M_B[h]:C#1",
            "rxn1 M_A[h]:N#1=M_B[h]:N#1",
            "rxn2 M_C[h]:O#1=M_D[h]:O#1",
            "rxn2 M_C[h]:N#1=M_D[h]:N#1",
        ]
        result = unite_mappings.filter_nitrogen_mappings(lines)
        assert len(result) == 2
        assert all(":N#" in line for line in result)

    def test_no_nitrogen(self):
        lines = ["rxn M_A[h]:C#1=M_B[h]:C#1"]
        assert unite_mappings.filter_nitrogen_mappings(lines) == []

    def test_all_nitrogen(self):
        lines = ["rxn M_A[h]:N#1=M_B[h]:N#1"]
        assert unite_mappings.filter_nitrogen_mappings(lines) == lines


class TestCountPerReaction:
    def test_basic_count(self):
        n_lines = [
            "rxn1 M_A:N#1=M_B:N#1",
            "rxn1 M_A:N#2=M_B:N#2",
            "rxn2 M_C:N#1=M_D:N#1",
        ]
        result = unite_mappings.count_per_reaction(n_lines)
        lines = result.strip().split("\n")
        assert len(lines) == 2
        assert lines[0].endswith("rxn2")
        assert "1" in lines[0]
        assert lines[1].endswith("rxn1")
        assert "2" in lines[1]

    def test_uniq_c_format(self):
        n_lines = ["rxn A:N#1=B:N#1"] * 5
        result = unite_mappings.count_per_reaction(n_lines)
        assert result.startswith("      5 rxn")

    def test_empty(self):
        assert unite_mappings.count_per_reaction([]) == ""


class TestMakeHistogram:
    def test_basic_histogram(self):
        count_text = "      1 A\n      1 B\n      2 C\n      2 D\n      3 E\n"
        result = unite_mappings.make_histogram(count_text)
        lines = result.rstrip("\n").split("\n")
        assert len(lines) == 3
        assert lines[0] == "      2 1"
        assert lines[1] == "      2 2"
        assert lines[2] == "      1 3"

    def test_uniq_c_format(self):
        count_text = "      1 A\n      2 B\n"
        result = unite_mappings.make_histogram(count_text)
        lines = result.strip().split("\n")
        for line in lines:
            count_val = int(line.strip().split()[0])
            assert 1 <= count_val <= 10

    def test_empty(self):
        assert unite_mappings.make_histogram("") == ""


class TestExtractNitrogenAtoms:
    def test_basic_extraction(self):
        n_lines = [
            "rxn1 M_A[h]:N#1=M_B[h]:N#1",
            "rxn1 M_A[h]:N#2=M_B[h]:N#2",
        ]
        atoms = unite_mappings.extract_nitrogen_atoms(n_lines)
        assert len(atoms) == 4
        assert "M_A[h]:N#1" in atoms
        assert "M_B[h]:N#1" in atoms

    def test_empty(self):
        assert unite_mappings.extract_nitrogen_atoms([]) == []


class TestCountAtoms:
    def test_unique_sort(self):
        atoms = ["B:N#2", "A:N#1", "B:N#2", "A:N#1"]
        result = unite_mappings.count_atoms(atoms, unique=True)
        lines = result.strip().split("\n")
        assert lines == ["A:N#1", "B:N#2"]

    def test_count_sort(self):
        atoms = ["A:N#1", "A:N#1", "B:N#2", "C:N#1", "C:N#1", "C:N#1"]
        result = unite_mappings.count_atoms(atoms, unique=False)
        lines = result.rstrip("\n").split("\n")
        assert len(lines) == 3
        assert lines[0].endswith("B:N#2")
        assert lines[0].strip().startswith("1")
        assert lines[1].endswith("A:N#1")
        assert lines[1].strip().startswith("2")
        assert lines[2].endswith("C:N#1")
        assert lines[2].strip().startswith("3")

    def test_uniq_c_width(self):
        atoms = ["A:N#1"] * 5
        result = unite_mappings.count_atoms(atoms, unique=False)
        assert result.startswith("      5 A:N#1")

    def test_empty(self):
        assert unite_mappings.count_atoms([], unique=True) == ""
        assert unite_mappings.count_atoms([], unique=False) == ""


# ── Integration tests ────────────────────────────────────────────────


@pytest.mark.integration
class TestShellEquivalence:
    """Prove the Python port produces identical output to unite_mappings.sh."""

    def test_python_matches_shell_all_files(self, tmp_path):
        unite_mappings.unite_mappings(REACTIONS_DIR, tmp_path)
        failures = []
        for golden in GOLDEN_UNITE_DIR.iterdir():
            if not golden.is_file():
                continue
            actual = tmp_path / golden.name
            if not actual.exists():
                failures.append(f"{golden.name}: missing from Python output")
                continue
            expected_text = golden.read_text()
            actual_text = actual.read_text()
            if expected_text != actual_text:
                failures.append(f"{golden.name}: content mismatch")
        assert failures == [], (
            f"{len(failures)} mismatches:\n" + "\n".join(failures)
        )

    def test_python_matches_shell_per_reaction_sorted(self, tmp_path):
        unite_mappings.unite_mappings(REACTIONS_DIR, tmp_path)
        golden = GOLDEN_UNITE_DIR / "all_mapping.N.sorted.txt"
        actual = tmp_path / "all_mapping.N.sorted.txt"
        assert golden.read_text() == actual.read_text()

    def test_python_matches_shell_rxn_counts(self, tmp_path):
        unite_mappings.unite_mappings(REACTIONS_DIR, tmp_path)
        for name in [
            "all_rxn_N_count.txt",
            "all_rxn_N_count.histo",
        ]:
            golden = GOLDEN_UNITE_DIR / name
            actual = tmp_path / name
            assert golden.read_text() == actual.read_text(), f"{name} mismatch"

    def test_python_matches_shell_atom_stats(self, tmp_path):
        unite_mappings.unite_mappings(REACTIONS_DIR, tmp_path)
        for name in [
            "all_atoms.N.sorted.txt",
            "all_atoms.N.count.txt",
            "all_atoms.N.count.histo",
        ]:
            golden = GOLDEN_UNITE_DIR / name
            actual = tmp_path / name
            assert golden.read_text() == actual.read_text(), f"{name} mismatch"

    def test_python_matches_shell_global_mapping(self, tmp_path):
        unite_mappings.unite_mappings(REACTIONS_DIR, tmp_path)
        golden = GOLDEN_UNITE_DIR / "all_mapping.sorted.txt"
        actual = tmp_path / "all_mapping.sorted.txt"
        assert golden.read_text() == actual.read_text()
