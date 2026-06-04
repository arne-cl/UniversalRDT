from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Callable

from run_rdt import (
    process_reaction_data,
    ReactionProcessingResult,
    MappingEntry,
)
from unite_mappings import _resolve_reactions_dir


@dataclass
class LoadResult:
    reactions: dict[str, ReactionProcessingResult]
    skipped: list[tuple[str, str]]


@dataclass(frozen=True)
class EquivalenceClass:
    reactions: tuple[str, ...]
    mapping: frozenset[tuple[str, str]]


@dataclass
class PairInconsistency:
    from_species: str
    to_species: str
    classes: list[EquivalenceClass]


@dataclass
class LabelViolation:
    species_id: str
    inchi_strings: dict[str, str]


@dataclass
class ConsistencyResult:
    inconsistencies: list[PairInconsistency]
    label_violations: list[LabelViolation]
    excluded_species: set[str]


def identity_exact(species_id: str) -> str:
    return species_id


def _extract_pair_mapping(
    result: ReactionProcessingResult,
    from_species: str,
    to_species: str,
) -> frozenset[tuple[str, str]]:
    from_entries = [
        e for e in result.all_entries
        if e.side == "from" and e.species_id == from_species and e.element != "H"
    ]
    to_entries = [
        e for e in result.all_entries
        if e.side == "to" and e.species_id == to_species and e.element != "H"
    ]
    to_by_index = {e.rdt_atom_index: e for e in to_entries}
    pairs = []
    for f in from_entries:
        t = to_by_index.get(f.rdt_atom_index)
        if t is not None:
            pairs.append((f"{f.element}#{f.element_index}", f"{t.element}#{t.element_index}"))
    return frozenset(pairs)


def _extract_all_pair_mappings(
    results: dict[str, ReactionProcessingResult],
    identity_fn: Callable[[str], str],
    skip_self_pairs: bool,
    excluded_species: set[str],
) -> dict[tuple[str, str], dict[str, frozenset[tuple[str, str]]]]:
    pair_mappings: dict[tuple[str, str], dict[str, frozenset[tuple[str, str]]]] = {}
    for rxn_name, result in results.items():
        from_species_ids = set()
        to_species_ids = set()
        for mol in result.molecules:
            sid = identity_fn(mol.species_id)
            if sid in excluded_species:
                continue
            if mol.side == "from":
                from_species_ids.add(sid)
            else:
                to_species_ids.add(sid)
        for f in from_species_ids:
            for t in to_species_ids:
                if skip_self_pairs and f == t:
                    continue
                from_orig = [m.species_id for m in result.molecules if m.side == "from" and identity_fn(m.species_id) == f]
                to_orig = [m.species_id for m in result.molecules if m.side == "to" and identity_fn(m.species_id) == t]
                mapping = frozenset()
                for fo in from_orig:
                    for to_ in to_orig:
                        m = _extract_pair_mapping(result, fo, to_)
                        mapping = mapping | m
                key = (f, t)
                pair_mappings.setdefault(key, {})
                pair_mappings[key][rxn_name] = mapping
    return pair_mappings


def _verify_inchi_labels(
    results: dict[str, ReactionProcessingResult],
) -> tuple[list[LabelViolation], set[str]]:
    species_inchis: dict[str, dict[str, str]] = {}
    for rxn_name, result in results.items():
        for mol in result.molecules:
            if not mol.inchi:
                continue
            species_inchis.setdefault(mol.species_id, {})
            species_inchis[mol.species_id][rxn_name] = mol.inchi
    violations = []
    excluded = set()
    for species_id, inchis_by_rxn in species_inchis.items():
        unique = set(inchis_by_rxn.values())
        if len(unique) > 1:
            violations.append(LabelViolation(
                species_id=species_id,
                inchi_strings=dict(inchis_by_rxn),
            ))
            excluded.add(species_id)
    return violations, excluded


def find_inconsistencies(
    results: dict[str, ReactionProcessingResult],
    identity_fn: Callable[[str], str] = identity_exact,
    skip_self_pairs: bool = True,
) -> ConsistencyResult:
    violations, excluded = _verify_inchi_labels(results)
    pair_mappings = _extract_all_pair_mappings(
        results, identity_fn, skip_self_pairs, excluded,
    )
    inconsistencies = []
    for (from_sp, to_sp), rxn_map in pair_mappings.items():
        by_mapping: dict[frozenset[tuple[str, str]], list[str]] = {}
        for rxn_name, mapping in rxn_map.items():
            by_mapping.setdefault(mapping, [])
            by_mapping[mapping].append(rxn_name)
        if len(by_mapping) > 1:
            classes = []
            for mapping, rxn_names in by_mapping.items():
                classes.append(EquivalenceClass(
                    reactions=tuple(sorted(rxn_names)),
                    mapping=mapping,
                ))
            classes.sort(key=lambda c: -len(c.reactions))
            inconsistencies.append(PairInconsistency(
                from_species=from_sp,
                to_species=to_sp,
                classes=classes,
            ))
    return ConsistencyResult(
        inconsistencies=inconsistencies,
        label_violations=violations,
        excluded_species=excluded,
    )


def load_all_reactions(
    reactions_dir: Path,
) -> LoadResult:
    resolved = _resolve_reactions_dir(reactions_dir)
    reactions: dict[str, ReactionProcessingResult] = {}
    skipped: list[tuple[str, str]] = []
    subdirs = sorted(d for d in resolved.iterdir() if d.is_dir())
    for i, rxn_dir in enumerate(subdirs, 1):
        try:
            result = process_reaction_data(rxn_dir)
            reactions[result.rxn_name] = result
        except Exception as e:
            logging.warning("Skipping %s: %s", rxn_dir.name, e)
            skipped.append((rxn_dir.name, str(e)))
        if i % 50 == 0:
            print(f"Processed {i} reactions...")
    return LoadResult(reactions=reactions, skipped=skipped)
