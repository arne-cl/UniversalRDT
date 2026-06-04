"""Molecule-level view of a single reaction for display in Jupyter notebooks.

Builds on `run_rdt.process_reaction_data` — no file I/O beyond reading
the source RXN + species data.
"""
from collections import Counter
from dataclasses import dataclass
from html import escape
from pathlib import Path

from run_rdt import (
    MoleculeProcessingResult,
    ReactionProcessingResult,
    MappingEntry,
    process_reaction_data,
)


# ── Dataclasses ──────────────────────────────────────────────────────


@dataclass
class MoleculeAtomView:
    """One atom as displayed in a molecule table.

    Attributes:
        rdt_index: Global RDT atom index.
        element: Chemical element symbol.
        mdl_line_number: 1-based position in the MDL V2000 atom table.
        inchi_position: 1-based position in InChI canonical order.
        inchi_element_index: Per-element counter in InChI order (C#1, C#2, ...).
        inchi_label: Full InChI label, e.g. "M_GAP[h]:C#1=".
        species_id: Identified species with compartment.
        compartment: Subcellular compartment tag.
        side: "from" (reactant) or "to" (product).
    """
    rdt_index: int
    element: str
    mdl_line_number: int
    inchi_position: int
    inchi_element_index: int
    inchi_label: str
    species_id: str
    compartment: str
    side: str


@dataclass
class MoleculeView:
    """All atoms of one molecule, presented in both orderings.

    Attributes:
        species_id: Identified species with compartment.
        compartment: Subcellular compartment tag.
        side: "from" (reactant) or "to" (product).
        atoms_inchi_order: Atoms sorted by InChI canonical position.
        atoms_mdl_order: Atoms sorted by MDL line number.
    """
    species_id: str
    compartment: str
    side: str
    atoms_inchi_order: list[MoleculeAtomView]
    atoms_mdl_order: list[MoleculeAtomView]

    def __str__(self) -> str:
        return self._format_text()

    def _format_text(self) -> str:
        n = len(self.atoms_inchi_order)
        lines: list[str] = []
        lines.append(f"Molecule: {self.species_id} ({n} atoms)")
        lines.append("")

        header_mdl = "RDT (MDL) order:"
        header_inchi = "InChI canonical order:"
        gap = 4
        lines.append(f"{header_mdl:<{30}}{' ' * gap}{header_inchi}")
        col_hdr = f"{'#':>3}  {'Elem':>4}  {'RDT#':>4}  {'InChI':<10}  {' ' * gap}{'#':>3}  {'Elem':>4}  {'RDT#':>4}  {'Label':<12}"
        lines.append(col_hdr)

        for mdl_atom, inchi_atom in zip(self.atoms_mdl_order, self.atoms_inchi_order):
            left = f"{mdl_atom.mdl_line_number:>3}  {mdl_atom.element:>4}  {mdl_atom.rdt_index:>4}  {mdl_atom.element}#{mdl_atom.inchi_element_index:<5}"
            right = f"{inchi_atom.inchi_position:>3}  {inchi_atom.element:>4}  {inchi_atom.rdt_index:>4}  {inchi_atom.element}#{inchi_atom.inchi_element_index}"
            lines.append(f"{left:<35}{' ' * gap}{right}")

        return "\n".join(lines)

    def _repr_html_(self) -> str:
        return self._format_html()

    def _format_html(self) -> str:
        th = 'style="padding: 2px 8px; text-align: left; border-bottom: 1px solid #ccc;"'
        parts: list[str] = []
        parts.append('<div style="font-family: monospace; margin: 4px 0;">')
        parts.append(
            f"<strong>{escape(self.species_id)}</strong> "
            f"({len(self.atoms_inchi_order)} atoms)"
        )
        parts.append('<table style="border-collapse: collapse; margin: 4px 0;">')
        parts.append(
            f"<tr><th {th}>#</th><th {th}>Elem</th><th {th}>RDT#</th>"
            f"<th {th}>InChI</th>"
            f'<th style="padding: 2px 12px;"></th>'
            f"<th {th}>#</th><th {th}>Elem</th><th {th}>RDT#</th>"
            f"<th {th}>Label</th></tr>"
        )
        for mdl_atom, inchi_atom in zip(self.atoms_mdl_order, self.atoms_inchi_order):
            bg_mdl = _element_color(mdl_atom.element)
            bg_inchi = _element_color(inchi_atom.element)
            td_mdl = f'style="padding: 2px 8px; background: {bg_mdl};"'
            td_inchi = f'style="padding: 2px 8px; background: {bg_inchi};"'
            parts.append("<tr>")
            parts.append(f'<td {td_mdl}>{mdl_atom.mdl_line_number}</td>')
            parts.append(f'<td {td_mdl}>{escape(mdl_atom.element)}</td>')
            parts.append(f'<td {td_mdl}>{mdl_atom.rdt_index}</td>')
            parts.append(f'<td {td_mdl}>{escape(mdl_atom.element)}#{mdl_atom.inchi_element_index}</td>')
            parts.append('<td style="padding: 2px 6px; color: #999;">&harr;</td>')
            parts.append(f'<td {td_inchi}>{inchi_atom.inchi_position}</td>')
            parts.append(f'<td {td_inchi}>{escape(inchi_atom.element)}</td>')
            parts.append(f'<td {td_inchi}>{inchi_atom.rdt_index}</td>')
            parts.append(f'<td {td_inchi}>{escape(inchi_atom.element)}#{inchi_atom.inchi_element_index}</td>')
            parts.append("</tr>")
        parts.append("</table>")
        parts.append("</div>")
        return "".join(parts)


@dataclass
class ReactionView:
    """Complete view of a single reaction.

    Attributes:
        rxn_name: Reaction name (folder name).
        smiles: Reaction SMILES.
        molecules: Per-molecule views, grouped by side.
        mapping_text: Full atom-to-atom mapping string.
        summary: One-line summary string.
    """
    rxn_name: str
    smiles: str
    molecules: list[MoleculeView]
    mapping_text: str
    summary: str

    def __str__(self) -> str:
        return self._format_text()

    def _format_text(self) -> str:
        lines: list[str] = []
        lines.append(f"=== {self.rxn_name} ===")
        lines.append(f"SMILES: {self.smiles}")
        lines.append("")

        reactants = [m for m in self.molecules if m.side == "from"]
        products = [m for m in self.molecules if m.side == "to"]

        if reactants:
            lines.append("--- Reactants ---")
            lines.append("")
            for m in reactants:
                lines.append(str(m))
                lines.append("")

        if products:
            lines.append("--- Products ---")
            lines.append("")
            for m in products:
                lines.append(str(m))
                lines.append("")

        lines.append("--- Transitions ---")
        lines.append(self.summary)
        if self.mapping_text:
            lines.append("")
            for pair in self.mapping_text.split(","):
                lines.append(f"  {pair}")

        return "\n".join(lines)

    def _repr_html_(self) -> str:
        return self._format_html()

    def _format_html(self) -> str:
        parts: list[str] = []
        parts.append('<div style="font-family: monospace; margin: 8px 0;">')
        parts.append(f"<h3>{escape(self.rxn_name)}</h3>")
        parts.append(f"<p>SMILES: <code>{escape(self.smiles)}</code></p>")

        for section_name, side_filter in [("Reactants", "from"), ("Products", "to")]:
            section_mols = [m for m in self.molecules if m.side == side_filter]
            if not section_mols:
                continue
            parts.append(f"<h4>{section_name}</h4>")
            for mol in section_mols:
                parts.append(mol._repr_html_())

        parts.append("<h4>Transitions</h4>")
        parts.append(f"<p><em>{escape(self.summary)}</em></p>")
        if self.mapping_text:
            parts.append('<table style="border-collapse: collapse;">')
            for pair in self.mapping_text.split(","):
                parts.append(
                    '<tr><td style="padding: 2px 10px;">'
                    f"{escape(pair)}</td></tr>"
                )
            parts.append("</table>")

        parts.append("</div>")
        return "".join(parts)


# ── Element colors (moved from run_rdt.py) ────────────────────────────


_ELEMENT_COLORS: dict[str, str] = {
    "C": "#e3f2fd",
    "N": "#e8f5e9",
    "O": "#fce4ec",
    "P": "#fff3e0",
    "S": "#f3e5f5",
    "H": "#f5f5f5",
}
_DEFAULT_ELEMENT_COLOR = "#fafafa"


def _element_color(element: str) -> str:
    return _ELEMENT_COLORS.get(element, _DEFAULT_ELEMENT_COLOR)


# ── Helpers ───────────────────────────────────────────────────────────


def _atom_label(entry: MappingEntry) -> str:
    return f"{entry.element}#{entry.element_index}"


def _build_molecule_view(mol: MoleculeProcessingResult) -> MoleculeView:
    """Build a MoleculeView from a MoleculeProcessingResult."""
    atoms_inchi: list[MoleculeAtomView] = []
    for inchi_pos, (mdl_line, entry) in enumerate(zip(mol.inchi_order, mol.entries), 1):
        atoms_inchi.append(MoleculeAtomView(
            rdt_index=entry.rdt_atom_index,
            element=entry.element,
            mdl_line_number=mdl_line,
            inchi_position=inchi_pos,
            inchi_element_index=entry.element_index,
            inchi_label=entry.label,
            species_id=entry.species_id,
            compartment=entry.compartment,
            side=entry.side,
        ))

    line_to_entry = dict(zip(mol.inchi_order, mol.entries))
    line_to_inchi_pos = {line: i + 1 for i, line in enumerate(mol.inchi_order)}

    atoms_mdl: list[MoleculeAtomView] = []
    for mdl_pos, (element, rdt_idx) in enumerate(mol.rdt_index, 1):
        entry = line_to_entry[mdl_pos]
        inchi_pos = line_to_inchi_pos[mdl_pos]
        atoms_mdl.append(MoleculeAtomView(
            rdt_index=rdt_idx,
            element=element,
            mdl_line_number=mdl_pos,
            inchi_position=inchi_pos,
            inchi_element_index=entry.element_index,
            inchi_label=entry.label,
            species_id=entry.species_id,
            compartment=entry.compartment,
            side=entry.side,
        ))

    return MoleculeView(
        species_id=mol.species_id,
        compartment=mol.compartment,
        side=mol.side,
        atoms_inchi_order=atoms_inchi,
        atoms_mdl_order=atoms_mdl,
    )


def _compute_summary(all_entries: list[MappingEntry]) -> str:
    """One-line summary: pair count + element breakdown (from-side atoms only)."""
    non_h = [e for e in all_entries if e.element != "H"]
    if not non_h:
        return "No atom pairs."
    from_by_idx: dict[int, MappingEntry] = {
        e.rdt_atom_index: e for e in non_h if e.side == "from"
    }
    to_indices = {e.rdt_atom_index for e in non_h if e.side == "to"}
    paired = sorted(from_by_idx.keys() & to_indices)
    counts: Counter[str] = Counter()
    for idx in paired:
        counts[from_by_idx[idx].element] += 1
    elem_parts = [f"{counts[el]} {el}" for el in sorted(counts)]
    return f"{len(paired)} atom pairs ({', '.join(elem_parts)})"


# ── Main entry point ─────────────────────────────────────────────────


def build_reaction_view(rxn_dir: Path) -> ReactionView:
    """Build a ReactionView from a reaction folder.

    Args:
        rxn_dir: Path to a reaction subfolder with RDT output and species data.

    Returns:
        ``ReactionView`` ready for display.

    Raises:
        FileNotFoundError: If the RXN file is missing or empty.
        ValueError: If the RXN file has no ``$MOL`` marker.
        SubprocessError: If obabel fails.
    """
    result = process_reaction_data(rxn_dir)

    molecule_views = [_build_molecule_view(mol) for mol in result.molecules]
    summary = _compute_summary(result.all_entries)

    return ReactionView(
        rxn_name=result.rxn_name,
        smiles=result.smiles,
        molecules=molecule_views,
        mapping_text=result.mapping_text,
        summary=summary,
    )
