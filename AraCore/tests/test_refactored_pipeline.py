"""Tests for the refactored pipeline: process_reaction_data and postprocess_reaction wrapper."""
import pytest
import run_rdt


class TestReactionProcessingResult:
    """Dataclass exists and has expected structure."""

    def test_can_instantiate_minimal(self):
        """Minimal smoke test that the dataclass is importable and usable."""
        result = run_rdt.ReactionProcessingResult(
            rxn_name="test",
            smiles="C>>C",
            from_num=1,
            to_num=1,
            molecules=[],
            mapping_lines_text="",
            mapping_text="",
        )
        assert result.rxn_name == "test"
        assert result.all_entries == []


class TestMoleculeProcessingResult:
    """MoleculeProcessingResult exists and has expected fields."""

    def test_can_instantiate(self):
        mol = run_rdt.MoleculeProcessingResult(
            mol_num=1,
            rdt_index=[("C", 1)],
            inchi_order=[1],
            entries=[run_rdt.MappingEntry(1, "from", "M_X[h]", "h", "C", 1)],
            species_id="M_X[h]",
            compartment="h",
            side="from",
        )
        assert mol.mol_num == 1
        assert mol.species_id == "M_X[h]"


class TestProcessReactionData:
    """process_reaction_data returns structured Result or None."""

    def test_returns_none_for_nonexistent_dir(self, tmp_path):
        assert run_rdt.process_reaction_data(tmp_path / "nope") is None

    def test_returns_none_for_empty_dir(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        assert run_rdt.process_reaction_data(d) is None

    @pytest.mark.integration
    def test_returns_reaction_processing_result(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert result is not None
        assert isinstance(result, run_rdt.ReactionProcessingResult)

    @pytest.mark.integration
    def test_rxn_name_from_folder(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert result.rxn_name == "FBPA_h"

    @pytest.mark.integration
    def test_header_counts(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert result.from_num == 2
        assert result.to_num == 1

    @pytest.mark.integration
    def test_three_molecules(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert len(result.molecules) == 3

    @pytest.mark.integration
    def test_molecule_species_ids(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert result.molecules[0].species_id == "M_GAP[h]"
        assert result.molecules[1].species_id == "M_DHAP[h]"
        assert result.molecules[2].species_id == "M_FBP[h]"

    @pytest.mark.integration
    def test_molecule_sides(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert result.molecules[0].side == "from"
        assert result.molecules[1].side == "from"
        assert result.molecules[2].side == "to"

    @pytest.mark.integration
    def test_molecule_inchi_order(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert result.molecules[0].inchi_order == [2, 5, 3, 1, 4, 8, 9, 10, 6, 7]

    @pytest.mark.integration
    def test_molecule_entry_counts(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert len(result.molecules[0].entries) == 10
        assert len(result.molecules[1].entries) == 10
        assert len(result.molecules[2].entries) == 20

    @pytest.mark.integration
    def test_molecule_rdt_index_lengths(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert len(result.molecules[0].rdt_index) == 10
        assert len(result.molecules[1].rdt_index) == 10
        assert len(result.molecules[2].rdt_index) == 20

    @pytest.mark.integration
    def test_all_entries_derived_property(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        manual = [e for mol in result.molecules for e in mol.entries]
        assert result.all_entries == manual
        assert len(result.all_entries) == 40

    @pytest.mark.integration
    def test_mapping_lines_text_has_40_lines(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        lines = result.mapping_lines_text.strip().split("\n")
        assert len(lines) == 40
        for line in lines:
            parts = line.split("\t")
            assert len(parts) == 3
            int(parts[0])
            assert parts[1] in ("from", "to")

    @pytest.mark.integration
    def test_mapping_text_matches_golden(self, sample_rxn_dir, golden_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        golden = (golden_dir / "FBPA_h.mapping.txt").read_text().strip()
        assert result.mapping_text == golden

    @pytest.mark.integration
    def test_smiles_present(self, sample_rxn_dir):
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        assert ">>" in result.smiles


class TestPostprocessReactionWrapper:
    """After refactoring, postprocess_reaction must still produce correct output."""

    @pytest.mark.integration
    def test_still_returns_bool_and_strings(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        success, lines, mapping = run_rdt.postprocess_reaction(sample_rxn_dir)
        assert success is True
        assert isinstance(lines, str)
        assert isinstance(mapping, str)

    @pytest.mark.integration
    def test_still_writes_mapping_txt(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        run_rdt.postprocess_reaction(sample_rxn_dir)
        assert (sample_rxn_dir / "mapping.txt").exists()

    @pytest.mark.integration
    def test_consistent_with_process_reaction_data(self, sample_rxn_dir):
        _strip_generated(sample_rxn_dir)
        result = run_rdt.process_reaction_data(sample_rxn_dir)
        success, lines, mapping = run_rdt.postprocess_reaction(sample_rxn_dir)
        assert success is True
        assert mapping == result.mapping_text
        assert lines == result.mapping_lines_text


def _strip_generated(rxn_dir):
    for pattern in ["mapping*", "MOL_*.mdl", "MOL_*.inchi",
                    "MOL_*.inchikey", "MOL_*.rdt_index", "MOL_*.species_id"]:
        for f in rxn_dir.glob(pattern):
            f.unlink()
