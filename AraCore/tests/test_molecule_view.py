"""Tests for molecule_view.py: build_reaction_view and view dataclasses."""
import pytest


class TestMoleculeAtomView:
    """MoleculeAtomView dataclass exists and is usable."""

    def test_can_instantiate(self):
        from molecule_view import MoleculeAtomView
        atom = MoleculeAtomView(
            rdt_index=1, element="C", mdl_line_number=1,
            inchi_position=1, inchi_element_index=1,
            inchi_label="M_X[h]:C#1=",
            species_id="M_X[h]", compartment="h", side="from",
        )
        assert atom.rdt_index == 1
        assert atom.element == "C"


class TestMoleculeView:
    """MoleculeView dataclass exists and is usable."""

    def test_can_instantiate(self):
        from molecule_view import MoleculeAtomView, MoleculeView
        atom = MoleculeAtomView(
            rdt_index=1, element="C", mdl_line_number=1,
            inchi_position=1, inchi_element_index=1,
            inchi_label="M_X[h]:C#1=",
            species_id="M_X[h]", compartment="h", side="from",
        )
        mol = MoleculeView(
            species_id="M_X[h]", compartment="h", side="from",
            atoms_inchi_order=[atom],
            atoms_mdl_order=[atom],
        )
        assert mol.species_id == "M_X[h]"
        assert len(mol.atoms_inchi_order) == 1

    def test_str_contains_both_columns(self):
        from molecule_view import MoleculeAtomView, MoleculeView
        atom = MoleculeAtomView(
            rdt_index=1, element="C", mdl_line_number=1,
            inchi_position=1, inchi_element_index=1,
            inchi_label="M_X[h]:C#1=",
            species_id="M_X[h]", compartment="h", side="from",
        )
        mol = MoleculeView(
            species_id="M_X[h]", compartment="h", side="from",
            atoms_inchi_order=[atom],
            atoms_mdl_order=[atom],
        )
        s = str(mol)
        assert "RDT (MDL) order" in s
        assert "InChI canonical order" in s
        assert "C#1" in s


class TestReactionView:
    """ReactionView dataclass exists and is usable."""

    def test_can_instantiate(self):
        from molecule_view import MoleculeAtomView, MoleculeView, ReactionView
        atom = MoleculeAtomView(
            rdt_index=1, element="C", mdl_line_number=1,
            inchi_position=1, inchi_element_index=1,
            inchi_label="M_X[h]:C#1=",
            species_id="M_X[h]", compartment="h", side="from",
        )
        mol = MoleculeView(
            species_id="M_X[h]", compartment="h", side="from",
            atoms_inchi_order=[atom],
            atoms_mdl_order=[atom],
        )
        view = ReactionView(
            rxn_name="test", smiles="C>>C",
            molecules=[mol], mapping_text="",
            summary="1 atom pair",
        )
        assert view.rxn_name == "test"
        assert len(view.molecules) == 1

    def test_str_contains_headers(self):
        from molecule_view import MoleculeAtomView, MoleculeView, ReactionView
        atom = MoleculeAtomView(
            rdt_index=1, element="C", mdl_line_number=1,
            inchi_position=1, inchi_element_index=1,
            inchi_label="M_X[h]:C#1=",
            species_id="M_X[h]", compartment="h", side="from",
        )
        mol = MoleculeView(
            species_id="M_X[h]", compartment="h", side="from",
            atoms_inchi_order=[atom],
            atoms_mdl_order=[atom],
        )
        view = ReactionView(
            rxn_name="test", smiles="C>>C",
            molecules=[mol], mapping_text="M_X[h]:C#1=M_Y[h]:C#1",
            summary="1 atom pair (1 C)",
        )
        s = str(view)
        assert view.rxn_name in s
        assert "SMILES" in s
        assert "Transitions" in s


class TestBuildReactionView:
    """Integration tests using real FBPA_h data."""

    @pytest.mark.integration
    def test_returns_reaction_view(self, sample_rxn_dir):
        from molecule_view import build_reaction_view, ReactionView
        view = build_reaction_view(sample_rxn_dir)
        assert isinstance(view, ReactionView)

    @pytest.mark.integration
    def test_reaction_name(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        assert view.rxn_name == "FBPA_h"

    @pytest.mark.integration
    def test_three_molecules(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        assert len(view.molecules) == 3

    @pytest.mark.integration
    def test_molecule_species_and_sides(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        assert view.molecules[0].species_id == "M_GAP[h]"
        assert view.molecules[0].side == "from"
        assert view.molecules[1].species_id == "M_DHAP[h]"
        assert view.molecules[1].side == "from"
        assert view.molecules[2].species_id == "M_FBP[h]"
        assert view.molecules[2].side == "to"

    @pytest.mark.integration
    def test_each_molecule_has_equal_atom_counts(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        for mol in view.molecules:
            assert len(mol.atoms_inchi_order) == len(mol.atoms_mdl_order)
            assert len(mol.atoms_inchi_order) > 0

    @pytest.mark.integration
    def test_atom_counts_per_molecule(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        assert len(view.molecules[0].atoms_inchi_order) == 10  # GAP
        assert len(view.molecules[1].atoms_inchi_order) == 10  # DHAP
        assert len(view.molecules[2].atoms_inchi_order) == 20  # FBP

    @pytest.mark.integration
    def test_molecule_str_contains_species(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        for mol in view.molecules:
            s = str(mol)
            assert mol.species_id in s

    @pytest.mark.integration
    def test_reaction_str_contains_smiles_and_transitions(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        s = str(view)
        assert "FBPA_h" in s
        assert "SMILES" in s
        assert "Transitions" in s

    @pytest.mark.integration
    def test_summary_has_pair_count(self, sample_rxn_dir):
        from molecule_view import build_reaction_view
        view = build_reaction_view(sample_rxn_dir)
        assert "20" in view.summary
        assert "C" in view.summary
