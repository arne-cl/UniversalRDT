import consistency
import pytest
import run_rdt
from unittest.mock import patch


def _make_entry(rdt_atom_index, side, species_id, element, element_index):
    return run_rdt.MappingEntry(
        rdt_atom_index=rdt_atom_index,
        side=side,
        species_id=species_id,
        compartment=species_id.split("[")[-1].rstrip("]") if "[" in species_id else "",
        element=element,
        element_index=element_index,
    )


def _make_mol(species_id, side, entries, inchi="InChI=1S/C3H7O6P"):
    return run_rdt.MoleculeProcessingResult(
        mol_num=1,
        rdt_index=[],
        inchi_order=[],
        entries=entries,
        species_id=species_id,
        compartment=species_id.split("[")[-1].rstrip("]") if "[" in species_id else "",
        side=side,
        inchi=inchi,
    )


def _make_rxn(rxn_name, from_mols, to_mols):
    return run_rdt.ReactionProcessingResult(
        rxn_name=rxn_name,
        smiles="",
        from_num=len(from_mols),
        to_num=len(to_mols),
        molecules=from_mols + to_mols,
        mapping_lines_text="",
        mapping_text="",
    )


class TestIdentityExact:
    def test_returns_unchanged(self):
        assert consistency.identity_exact("M_ATP[c]") == "M_ATP[c]"

    def test_no_compartment(self):
        assert consistency.identity_exact("M_ATP") == "M_ATP"

    def test_empty_string(self):
        assert consistency.identity_exact("") == ""


class TestDataModel:
    def test_load_result_fields(self):
        result = consistency.LoadResult(reactions={}, skipped=[("rxn", "error")])
        assert result.reactions == {}
        assert result.skipped == [("rxn", "error")]

    def test_equivalence_class_frozen(self):
        ec = consistency.EquivalenceClass(
            reactions=("A", "B"),
            mapping=frozenset({("C#1", "C#5")}),
        )
        with pytest.raises(AttributeError):
            ec.reactions = ("C",)

    def test_pair_inconsistency_fields(self):
        ec = consistency.EquivalenceClass(
            reactions=("A",), mapping=frozenset()
        )
        pi = consistency.PairInconsistency(
            from_species="M_ATP", to_species="M_ADP", classes=[ec]
        )
        assert pi.from_species == "M_ATP"
        assert len(pi.classes) == 1

    def test_label_violation(self):
        lv = consistency.LabelViolation(
            species_id="M_ATP[c]",
            inchi_strings={"rxn1": "InChI=1", "rxn2": "InChI=2"},
        )
        assert len(lv.inchi_strings) == 2

    def test_consistency_result_fields(self):
        cr = consistency.ConsistencyResult(
            inconsistencies=[],
            label_violations=[],
            excluded_species=set(),
        )
        assert cr.inconsistencies == []
        assert cr.excluded_species == set()


class TestExtractPairMapping:
    def test_simple_pair(self):
        from_entries = [
            _make_entry(1, "from", "M_ATP[c]", "C", 1),
            _make_entry(2, "from", "M_ATP[c]", "C", 2),
        ]
        to_entries = [
            _make_entry(1, "to", "M_ADP[c]", "C", 5),
            _make_entry(2, "to", "M_ADP[c]", "C", 3),
        ]
        from_mol = _make_mol("M_ATP[c]", "from", from_entries)
        to_mol = _make_mol("M_ADP[c]", "to", to_entries)
        rxn = _make_rxn("R1", [from_mol], [to_mol])

        result = consistency._extract_pair_mapping(rxn, "M_ATP[c]", "M_ADP[c]")
        assert result == frozenset({("C#1", "C#5"), ("C#2", "C#3")})

    def test_excludes_hydrogen(self):
        from_entries = [
            _make_entry(1, "from", "M_A[c]", "C", 1),
            _make_entry(2, "from", "M_A[c]", "H", 1),
        ]
        to_entries = [
            _make_entry(1, "to", "M_B[c]", "C", 1),
            _make_entry(2, "to", "M_B[c]", "H", 2),
        ]
        from_mol = _make_mol("M_A[c]", "from", from_entries)
        to_mol = _make_mol("M_B[c]", "to", to_entries)
        rxn = _make_rxn("R1", [from_mol], [to_mol])

        result = consistency._extract_pair_mapping(rxn, "M_A[c]", "M_B[c]")
        assert result == frozenset({("C#1", "C#1")})


class TestExtractAllPairMappings:
    def test_basic_cross_join(self):
        r1_from = _make_mol("M_A[c]", "from", [
            _make_entry(1, "from", "M_A[c]", "C", 1),
        ])
        r1_to = _make_mol("M_B[c]", "to", [
            _make_entry(1, "to", "M_B[c]", "C", 5),
        ])
        rxn = _make_rxn("R1", [r1_from], [r1_to])

        result = consistency._extract_all_pair_mappings(
            {"R1": rxn}, consistency.identity_exact, skip_self_pairs=True,
            excluded_species=set(),
        )
        assert ("M_A[c]", "M_B[c]") in result
        assert result[("M_A[c]", "M_B[c]")]["R1"] == frozenset({("C#1", "C#5")})

    def test_skip_self_pairs(self):
        mol = _make_mol("M_ATP[c]", "from", [
            _make_entry(1, "from", "M_ATP[c]", "C", 1),
        ])
        mol2 = _make_mol("M_ATP[c]", "to", [
            _make_entry(1, "to", "M_ATP[c]", "C", 2),
        ])
        rxn = _make_rxn("R1", [mol], [mol2])

        result = consistency._extract_all_pair_mappings(
            {"R1": rxn}, consistency.identity_exact, skip_self_pairs=True,
            excluded_species=set(),
        )
        assert ("M_ATP[c]", "M_ATP[c]") not in result

    def test_include_self_pairs_when_disabled(self):
        mol = _make_mol("M_ATP[c]", "from", [
            _make_entry(1, "from", "M_ATP[c]", "C", 1),
        ])
        mol2 = _make_mol("M_ATP[c]", "to", [
            _make_entry(1, "to", "M_ATP[c]", "C", 2),
        ])
        rxn = _make_rxn("R1", [mol], [mol2])

        result = consistency._extract_all_pair_mappings(
            {"R1": rxn}, consistency.identity_exact, skip_self_pairs=False,
            excluded_species=set(),
        )
        assert ("M_ATP[c]", "M_ATP[c]") in result

    def test_excluded_species_skipped(self):
        r1_from = _make_mol("M_A[c]", "from", [
            _make_entry(1, "from", "M_A[c]", "C", 1),
        ])
        r1_to = _make_mol("M_B[c]", "to", [
            _make_entry(1, "to", "M_B[c]", "C", 5),
        ])
        rxn = _make_rxn("R1", [r1_from], [r1_to])

        result = consistency._extract_all_pair_mappings(
            {"R1": rxn}, consistency.identity_exact, skip_self_pairs=True,
            excluded_species={"M_A[c]"},
        )
        assert ("M_A[c]", "M_B[c]") not in result

    def test_multiple_reactions_same_pair(self):
        def _build_rxn(name):
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)])
            return _make_rxn(name, [f], [t])

        rxns = {"R1": _build_rxn("R1"), "R2": _build_rxn("R2")}
        result = consistency._extract_all_pair_mappings(
            rxns, consistency.identity_exact, skip_self_pairs=True,
            excluded_species=set(),
        )
        pair_data = result[("M_A[c]", "M_B[c]")]
        assert "R1" in pair_data
        assert "R2" in pair_data


class TestVerifyInchiLabels:
    def test_no_violations_when_consistent(self):
        rxn1 = _make_rxn("R1",
            [_make_mol("M_ATP[c]", "from", [], inchi="InChI=1S/ATP")],
            [_make_mol("M_ADP[c]", "to", [], inchi="InChI=1S/ADP")],
        )
        rxn2 = _make_rxn("R2",
            [_make_mol("M_ATP[c]", "from", [], inchi="InChI=1S/ATP")],
            [_make_mol("M_ADP[c]", "to", [], inchi="InChI=1S/ADP")],
        )
        violations, excluded = consistency._verify_inchi_labels({"R1": rxn1, "R2": rxn2})
        assert violations == []
        assert excluded == set()

    def test_violation_when_inchi_differs(self):
        rxn1 = _make_rxn("R1",
            [_make_mol("M_ATP[c]", "from", [], inchi="InChI=1S/ATP_v1")],
            [],
        )
        rxn2 = _make_rxn("R2",
            [_make_mol("M_ATP[c]", "from", [], inchi="InChI=1S/ATP_v2")],
            [],
        )
        violations, excluded = consistency._verify_inchi_labels({"R1": rxn1, "R2": rxn2})
        assert len(violations) == 1
        assert violations[0].species_id == "M_ATP[c]"
        assert "R1" in violations[0].inchi_strings
        assert "R2" in violations[0].inchi_strings
        assert "M_ATP[c]" in excluded

    def test_species_in_only_one_reaction(self):
        rxn1 = _make_rxn("R1",
            [_make_mol("M_ATP[c]", "from", [], inchi="InChI=1S/ATP")],
            [],
        )
        rxn2 = _make_rxn("R2",
            [_make_mol("M_ADP[c]", "from", [], inchi="InChI=1S/ADP")],
            [],
        )
        violations, excluded = consistency._verify_inchi_labels({"R1": rxn1, "R2": rxn2})
        assert violations == []
        assert excluded == set()

    def test_empty_inchi_not_flagged(self):
        rxn1 = _make_rxn("R1",
            [_make_mol("M_X[c]", "from", [], inchi="")],
            [],
        )
        violations, excluded = consistency._verify_inchi_labels({"R1": rxn1})
        assert violations == []
        assert excluded == set()


class TestFindInconsistencies:
    def _build_consistent_rxns(self):
        def _build(name):
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)])
            return _make_rxn(name, [f], [t])
        return {"R1": _build("R1"), "R2": _build("R2")}

    def test_no_inconsistencies(self):
        result = consistency.find_inconsistencies(self._build_consistent_rxns())
        assert result.inconsistencies == []

    def test_inconsistency_detected(self):
        def _build_rxn1():
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)])
            return _make_rxn("R1", [f], [t])

        def _build_rxn2():
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 3)])
            return _make_rxn("R2", [f], [t])

        result = consistency.find_inconsistencies({"R1": _build_rxn1(), "R2": _build_rxn2()})
        assert len(result.inconsistencies) == 1
        inc = result.inconsistencies[0]
        assert inc.from_species == "M_A[c]"
        assert inc.to_species == "M_B[c]"
        assert len(inc.classes) == 2

    def test_label_violation_excludes_species(self):
        rxn1 = _make_rxn("R1",
            [_make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)], inchi="InChI=1S/A_v1")],
            [_make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)], inchi="InChI=1S/B")],
        )
        rxn2 = _make_rxn("R2",
            [_make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)], inchi="InChI=1S/A_v2")],
            [_make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)], inchi="InChI=1S/B")],
        )
        result = consistency.find_inconsistencies({"R1": rxn1, "R2": rxn2})
        assert len(result.label_violations) == 1
        assert "M_A[c]" in result.excluded_species
        assert result.inconsistencies == []

    def test_equivalence_classes_grouped(self):
        def _build_a():
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)])
            return _make_rxn("R_A", [f], [t])

        def _build_b():
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)])
            return _make_rxn("R_B", [f], [t])

        def _build_c():
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 3)])
            return _make_rxn("R_C", [f], [t])

        result = consistency.find_inconsistencies({
            "R_A": _build_a(), "R_B": _build_b(), "R_C": _build_c(),
        })
        assert len(result.inconsistencies) == 1
        inc = result.inconsistencies[0]
        assert len(inc.classes) == 2
        class_a = [c for c in inc.classes if "R_A" in c.reactions][0]
        assert "R_B" in class_a.reactions
        assert "R_C" not in class_a.reactions

    def test_classes_sorted_by_size_descending(self):
        def _build(name):
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 5)])
            return _make_rxn(name, [f], [t])

        def _build_different(name):
            f = _make_mol("M_A[c]", "from", [_make_entry(1, "from", "M_A[c]", "C", 1)])
            t = _make_mol("M_B[c]", "to", [_make_entry(1, "to", "M_B[c]", "C", 3)])
            return _make_rxn(name, [f], [t])

        result = consistency.find_inconsistencies({
            "R1": _build("R1"), "R2": _build("R2"), "R3": _build("R3"),
            "R4": _build_different("R4"),
        })
        inc = result.inconsistencies[0]
        assert len(inc.classes[0].reactions) >= len(inc.classes[1].reactions)


class TestLoadAllReactions:
    def test_loads_from_directory(self, tmp_path):
        rxn_dir = tmp_path / "reactions" / "R1"
        rxn_dir.mkdir(parents=True)
        (rxn_dir / "ECBLAST_smiles_AAM.rxn").write_text("$RXN\n\n  1  1\n$MOL\n\nxxx\nV2000\n")
        (rxn_dir / "rxn.smiles").write_text("A>>B")
        (rxn_dir / "species_id_inchikey.txt").write_text("KEY1 M_A\n")
        (rxn_dir / "from_species_with_cmp").write_text("M_A[c]\n")
        (rxn_dir / "to_species_with_cmp").write_text("M_B[c]\n")

        mock_result = _make_rxn("R1", [], [])
        with patch("consistency.process_reaction_data", return_value=mock_result):
            result = consistency.load_all_reactions(tmp_path / "reactions")
        assert "R1" in result.reactions
        assert result.skipped == []

    def test_skips_failing_reactions(self, tmp_path):
        rxn_dir = tmp_path / "reactions" / "BadRxn"
        rxn_dir.mkdir(parents=True)
        (rxn_dir / "dummy.txt").write_text("x")

        with patch("consistency.process_reaction_data", side_effect=FileNotFoundError("no RXN")):
            result = consistency.load_all_reactions(tmp_path / "reactions")
        assert result.reactions == {}
        assert len(result.skipped) == 1
        assert result.skipped[0][0] == "BadRxn"
        assert "no RXN" in result.skipped[0][1]

    def test_mixed_success_and_failure(self, tmp_path):
        r1 = tmp_path / "reactions" / "R1"
        r1.mkdir(parents=True)
        r2 = tmp_path / "reactions" / "R2"
        r2.mkdir(parents=True)

        mock_result = _make_rxn("R1", [], [])

        def side_effect(path):
            if path.name == "R1":
                return mock_result
            raise FileNotFoundError("broken")

        with patch("consistency.process_reaction_data", side_effect=side_effect):
            result = consistency.load_all_reactions(tmp_path / "reactions")
        assert "R1" in result.reactions
        assert len(result.skipped) == 1
        assert result.skipped[0][0] == "R2"

    def test_empty_directory(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        result = consistency.load_all_reactions(empty)
        assert result.reactions == {}
        assert result.skipped == []

    def test_progress_print(self, tmp_path, capsys):
        for i in range(55):
            d = tmp_path / "reactions" / f"R{i:03d}"
            d.mkdir(parents=True)

        mock_result = _make_rxn("mock", [], [])
        with patch("consistency.process_reaction_data", return_value=mock_result):
            consistency.load_all_reactions(tmp_path / "reactions")
        output = capsys.readouterr().out
        assert "50" in output
