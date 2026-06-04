import pytest

import run_rdt


class TestMappingEntry:
    def test_label_from_side_ends_with_equals(self):
        entry = run_rdt.MappingEntry(1, "from", "M_X[h]", "C", 1)
        assert entry.label == "M_X[h]:C#1="

    def test_label_to_side_ends_with_comma(self):
        entry = run_rdt.MappingEntry(1, "to", "M_X[h]", "C", 1)
        assert entry.label == "M_X[h]:C#1,"

    def test_is_frozen(self):
        entry = run_rdt.MappingEntry(1, "from", "M_X[h]", "C", 1)
        with pytest.raises(AttributeError):
            entry.element = "N"


class TestBuildMappingEntries:
    def test_from_side_gap(self):
        rdt_index = [
            ("O", 1), ("C", 2), ("C", 3), ("O", 4), ("C", 5),
            ("O", 6), ("P", 7), ("O", 8), ("O", 9), ("O", 10),
        ]
        inchi_order = [2, 5, 3, 1, 4, 8, 9, 10, 6, 7]
        entries = run_rdt.build_mapping_entries(
            rdt_index, inchi_order, "M_GAP[h]", "from"
        )
        assert len(entries) == 10
        assert entries[0] == run_rdt.MappingEntry(2, "from", "M_GAP[h]", "C", 1)
        assert entries[1] == run_rdt.MappingEntry(5, "from", "M_GAP[h]", "C", 2)
        assert entries[4] == run_rdt.MappingEntry(4, "from", "M_GAP[h]", "O", 2)
        assert entries[9] == run_rdt.MappingEntry(7, "from", "M_GAP[h]", "P", 1)
        assert entries[0].label == "M_GAP[h]:C#1="
        assert entries[9].label == "M_GAP[h]:P#1="

    def test_to_side_fbp(self):
        rdt_index = [
            ("O", 1), ("C", 2), ("O", 14), ("C", 13), ("O", 10),
            ("O", 9), ("O", 6), ("C", 5), ("O", 20), ("O", 19),
            ("O", 16), ("C", 15), ("C", 12), ("O", 11), ("C", 3),
            ("O", 4), ("P", 17), ("O", 18), ("P", 7), ("O", 8),
        ]
        inchi_order = [12, 8, 13, 4, 2, 15, 3, 1, 16, 9, 10, 18, 5, 6, 20, 11, 7, 14, 17, 19]
        entries = run_rdt.build_mapping_entries(
            rdt_index, inchi_order, "M_FBP[h]", "to"
        )
        assert len(entries) == 20
        assert entries[0] == run_rdt.MappingEntry(15, "to", "M_FBP[h]", "C", 1)
        assert entries[-1] == run_rdt.MappingEntry(7, "to", "M_FBP[h]", "P", 2)
        assert entries[0].label == "M_FBP[h]:C#1,"
        assert entries[-1].label == "M_FBP[h]:P#2,"

    def test_single_atom(self):
        rdt_index = [("O", 42)]
        inchi_order = [1]
        entries = run_rdt.build_mapping_entries(
            rdt_index, inchi_order, "M_H2O[c]", "from"
        )
        assert len(entries) == 1
        assert entries[0] == run_rdt.MappingEntry(42, "from", "M_H2O[c]", "O", 1)

    def test_counter_independent_per_element(self):
        rdt_index = [("C", 1), ("N", 2), ("C", 3)]
        inchi_order = [1, 2, 3]
        entries = run_rdt.build_mapping_entries(
            rdt_index, inchi_order, "M_X[h]", "from"
        )
        assert entries[0].element_index == 1  # C#1
        assert entries[1].element_index == 1  # N#1
        assert entries[2].element_index == 2  # C#2


class TestAssembleMapping:
    def test_basic(self):
        entries = [
            run_rdt.MappingEntry(2, "from", "M_X[h]", "C", 1),
            run_rdt.MappingEntry(5, "to", "M_Y[h]", "C", 1),
        ]
        result = run_rdt.assemble_mapping(entries)
        assert result == "M_X[h]:C#1=M_Y[h]:C#1"

    def test_filters_hydrogen_by_element(self):
        entries = [
            run_rdt.MappingEntry(1, "from", "M_X[h]", "C", 1),
            run_rdt.MappingEntry(2, "from", "M_X[h]", "H", 1),
            run_rdt.MappingEntry(3, "to", "M_Y[h]", "C", 1),
        ]
        result = run_rdt.assemble_mapping(entries)
        assert "H" not in result
        assert result == "M_X[h]:C#1=M_Y[h]:C#1"

    def test_sorts_by_rdt_atom_index(self):
        entries = [
            run_rdt.MappingEntry(5, "from", "M_X[h]", "C", 1),
            run_rdt.MappingEntry(2, "from", "M_X[h]", "O", 1),
        ]
        result = run_rdt.assemble_mapping(entries)
        assert result == "M_X[h]:O#1=M_X[h]:C#1="

    def test_empty_entries(self):
        assert run_rdt.assemble_mapping([]) == ""
