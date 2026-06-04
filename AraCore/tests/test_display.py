import pytest
import run_rdt


# ── pair_entries unit tests ───────────────────────────────────────────


class TestPairEntries:
    """Group MappingEntry objects by rdt_atom_index into from/to pairs."""

    def test_basic_pair(self):
        from_entry = run_rdt.MappingEntry(1, "from", "M_X[h]", "h", "C", 1)
        to_entry = run_rdt.MappingEntry(1, "to", "M_Y[h]", "h", "C", 1)
        mapping = run_rdt.pair_entries([from_entry, to_entry])
        assert mapping.pairs == [(from_entry, to_entry)]
        assert mapping.unmapped_from == []
        assert mapping.unmapped_to == []

    def test_multiple_pairs(self):
        entries = [
            run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1),
            run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "C", 2),
            run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "C", 1),
            run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "C", 2),
        ]
        mapping = run_rdt.pair_entries(entries)
        assert len(mapping.pairs) == 2

    def test_sorted_by_rdt_index(self):
        entries = [
            run_rdt.MappingEntry(5, "from", "M_A[h]", "h", "C", 1),
            run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "O", 1),
            run_rdt.MappingEntry(5, "to", "M_B[h]", "h", "C", 1),
            run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "O", 1),
        ]
        mapping = run_rdt.pair_entries(entries)
        assert mapping.pairs == [
            (entries[1], entries[3]),  # rdt_idx=2 first
            (entries[0], entries[2]),  # rdt_idx=5 second
        ]

    def test_filters_hydrogen(self):
        entries = [
            run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "H", 1),
            run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "H", 1),
            run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "C", 1),
            run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "C", 1),
        ]
        mapping = run_rdt.pair_entries(entries)
        assert len(mapping.pairs) == 1

    def test_empty_entries(self):
        mapping = run_rdt.pair_entries([])
        assert mapping.pairs == []
        assert mapping.unmapped_from == []
        assert mapping.unmapped_to == []

    def test_unmapped_from(self):
        from_entry = run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1)
        to_entry = run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "O", 1)
        mapping = run_rdt.pair_entries([from_entry, to_entry])
        assert len(mapping.pairs) == 0
        assert mapping.unmapped_from == [from_entry]
        assert mapping.unmapped_to == [to_entry]

    def test_mixed_elements(self):
        c_from = run_rdt.MappingEntry(2, "from", "M_GAP[h]", "h", "C", 1)
        o_from = run_rdt.MappingEntry(1, "from", "M_GAP[h]", "h", "O", 1)
        c_to = run_rdt.MappingEntry(2, "to", "M_FBP[h]", "h", "C", 5)
        o_to = run_rdt.MappingEntry(1, "to", "M_FBP[h]", "h", "O", 2)
        mapping = run_rdt.pair_entries([c_from, o_from, c_to, o_to])
        assert len(mapping.pairs) == 2
        assert mapping.pairs[0] == (o_from, o_to)  # idx 1 first
        assert mapping.pairs[1] == (c_from, c_to)  # idx 2 second


# ── AtomMapping display unit tests ────────────────────────────────────


class TestAtomMappingStr:
    """__str__ produces aligned text table."""

    def test_basic_table(self):
        pairs = [
            (run_rdt.MappingEntry(2, "from", "M_GAP[h]", "h", "C", 1),
             run_rdt.MappingEntry(2, "to", "M_FBP[h]", "h", "C", 5)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[])
        s = str(mapping)
        assert "RDT" in s
        assert "M_GAP[h]" in s
        assert "M_FBP[h]" in s
        assert "C#1" in s
        assert "C#5" in s

    def test_mixed_elements(self):
        pairs = [
            (run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "O", 1),
             run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "O", 2)),
            (run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "C", 1),
             run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "C", 1)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[])
        s = str(mapping)
        assert "O#1" in s
        assert "O#2" in s
        assert "C#1" in s

    def test_empty_pairs(self):
        mapping = run_rdt.AtomMapping(pairs=[], unmapped_from=[], unmapped_to=[])
        s = str(mapping)
        assert "No atom pairs" in s or len(s) < 50


class TestAtomMappingSummary:
    """Summary provides quick overview."""

    def test_single_pair(self):
        pairs = [
            (run_rdt.MappingEntry(2, "from", "M_GAP[h]", "h", "C", 1),
             run_rdt.MappingEntry(2, "to", "M_FBP[h]", "h", "C", 5)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[])
        s = mapping.summary()
        assert "1" in s
        assert "C" in s

    def test_with_rxn_name(self):
        pairs = [
            (run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1),
             run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "C", 1)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[],
                                      rxn_name="FBPA_h")
        s = mapping.summary()
        assert "FBPA_h" in s

    def test_element_breakdown(self):
        pairs = [
            (run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1),
             run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "C", 1)),
            (run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "O", 1),
             run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "O", 1)),
            (run_rdt.MappingEntry(3, "from", "M_A[h]", "h", "N", 1),
             run_rdt.MappingEntry(3, "to", "M_B[h]", "h", "N", 1)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[])
        s = mapping.summary()
        assert "3" in s
        assert "C" in s
        assert "O" in s
        assert "N" in s

    def test_with_unmapped(self):
        pairs = [
            (run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1),
             run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "C", 1)),
        ]
        u_from = [run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "O", 1)]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=u_from, unmapped_to=[])
        s = mapping.summary()
        assert "1 unmapped" in s


class TestAtomMappingReprHtml:
    """_repr_html_ produces element-colored HTML table."""

    def test_contains_table(self):
        pairs = [
            (run_rdt.MappingEntry(2, "from", "M_GAP[h]", "h", "C", 1),
             run_rdt.MappingEntry(2, "to", "M_FBP[h]", "h", "C", 5)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[])
        html = mapping._repr_html_()
        assert "<table" in html
        assert "M_GAP[h]" in html
        assert "M_FBP[h]" in html

    def test_element_colors_present(self):
        pairs = [
            (run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1),
             run_rdt.MappingEntry(1, "to", "M_B[h]", "h", "C", 1)),
            (run_rdt.MappingEntry(2, "from", "M_A[h]", "h", "O", 1),
             run_rdt.MappingEntry(2, "to", "M_B[h]", "h", "O", 1)),
        ]
        mapping = run_rdt.AtomMapping(pairs=pairs, unmapped_from=[], unmapped_to=[])
        html = mapping._repr_html_()
        # Should contain color/style attributes
        assert "style=" in html or "background" in html.lower()

    def test_unmapped_section(self):
        u_from = [run_rdt.MappingEntry(1, "from", "M_A[h]", "h", "C", 1)]
        mapping = run_rdt.AtomMapping(pairs=[], unmapped_from=u_from, unmapped_to=[])
        html = mapping._repr_html_()
        assert "unmapped" in html.lower()

    def test_empty_mapping(self):
        mapping = run_rdt.AtomMapping(pairs=[], unmapped_from=[], unmapped_to=[])
        html = mapping._repr_html_()
        assert html  # Should produce some non-empty string


# ── Integration: run_rdt_and_map_from_folder ──────────────────────────


class TestRunRdtAndMapFromFolder:
    """Integration tests using extracted reaction data (no Java RDT needed)."""

    def test_returns_atom_mapping(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        mapping = run_rdt.run_rdt_and_map_from_folder(sample_rxn_dir)
        assert isinstance(mapping, run_rdt.AtomMapping)
        assert len(mapping.pairs) > 0

    def test_golden_fbpa_has_20_pairs(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        mapping = run_rdt.run_rdt_and_map_from_folder(sample_rxn_dir)
        # FBPA_h has 20 non-H atom transitions (10 GAP + 10 DHAP = 20,
        # mapping to 20 FBP atoms)
        assert len(mapping.pairs) == 20

    def test_pairs_are_sorted_by_rdt_index(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        mapping = run_rdt.run_rdt_and_map_from_folder(sample_rxn_dir)
        indices = [f.rdt_atom_index for f, _ in mapping.pairs]
        assert indices == sorted(indices)

    def test_summary_matches_pairs(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        mapping = run_rdt.run_rdt_and_map_from_folder(sample_rxn_dir)
        assert str(len(mapping.pairs)) in mapping.summary()

    def test_str_contains_all_species(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        mapping = run_rdt.run_rdt_and_map_from_folder(sample_rxn_dir)
        s = str(mapping)
        assert "M_GAP[h]" in s
        assert "M_DHAP[h]" in s
        assert "M_FBP[h]" in s

    def test_repr_html_contains_species(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        mapping = run_rdt.run_rdt_and_map_from_folder(sample_rxn_dir)
        html = mapping._repr_html_()
        assert "M_GAP[h]" in html
        assert "M_FBP[h]" in html


def _strip_generated(rxn_dir):
    for pattern in ["mapping*", "MOL_*.mdl", "MOL_*.inchi",
                    "MOL_*.inchikey", "MOL_*.rdt_index", "MOL_*.species_id"]:
        for f in rxn_dir.glob(pattern):
            f.unlink()
