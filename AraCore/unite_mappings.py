"""Python reimplementation of AraCore/unite_mappings.sh.

Collects per-reaction atom mappings into global mapping tables and
produces nitrogen-specific statistics (counts, histograms).
"""

import argparse
import locale
import sys
import tempfile
import zipfile
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

locale.setlocale(locale.LC_ALL, "")

UNIQ_C_WIDTH = 7


def split_and_sort_mapping(mapping_text: str) -> List[str]:
    """Split a single-line mapping by commas into sorted atom pairs.

    Replaces ``sed 's/,/\\n/g' mapping.txt | sort``.

    Args:
        mapping_text: Contents of a ``mapping.txt`` file (single line).

    Returns:
        Sorted list of individual atom-pair strings.
    """
    text = mapping_text.strip()
    if not text:
        return []
    pairs = text.split(",")
    return sorted(pairs, key=locale.strxfrm)


def collect_all_mappings(reactions_dir: Path) -> Tuple[List[str], Dict[str, List[str]]]:
    """Read every reaction's mapping.txt and collect prefixed pairs.

    Replaces the ``for rxn in …; do sed … | sort … ; cat … | sed … >> all_mapping.txt; done``
    loop.  Also writes per-reaction ``mapping.sorted.txt`` files.

    Args:
        reactions_dir: Path to the ``reaction_intermediates/`` directory.

    Returns:
        Tuple of (all_lines, per_reaction_sorted) where all_lines is a flat
        list of ``"<rxn_name> <pair>"`` strings and per_reaction_sorted maps
        each reaction name to its sorted pairs.
    """
    all_lines: List[str] = []
    per_reaction: Dict[str, List[str]] = {}

    for rxn_dir in sorted(reactions_dir.iterdir()):
        if not rxn_dir.is_dir():
            continue
        rxn_name = rxn_dir.name
        mapping_path = rxn_dir / "mapping.txt"
        if not mapping_path.exists():
            continue
        mapping_text = mapping_path.read_text()
        sorted_pairs = split_and_sort_mapping(mapping_text)
        per_reaction[rxn_name] = sorted_pairs

        sorted_path = rxn_dir / "mapping.sorted.txt"
        sorted_path.write_text("\n".join(sorted_pairs) + "\n" if sorted_pairs else "")

        for pair in sorted_pairs:
            all_lines.append(f"{rxn_name} {pair}")

    return all_lines, per_reaction


def sort_all_mappings(all_lines: List[str]) -> List[str]:
    """Sort all mapping lines lexicographically.

    Replaces ``sort all_mapping.txt > all_mapping.sorted.txt``.
    """
    return sorted(all_lines, key=locale.strxfrm)


def filter_nitrogen_mappings(sorted_lines: List[str]) -> List[str]:
    """Filter sorted mapping lines to only nitrogen atom pairs.

    Replaces ``grep ':N#' all_mapping.sorted.txt > all_mapping.N.sorted.txt``.
    """
    return [line for line in sorted_lines if ":N#" in line]


def count_per_reaction(n_lines: List[str]) -> str:
    """Count N-mapping lines per reaction, formatted like ``uniq -c | sort -n``.

    Replaces ``sed 's/ .*//' … | sort | uniq -c | sort -n``.

    Returns:
        String in ``uniq -c`` format (right-aligned count, space, value),
        sorted numerically by count then lexicographically by name.
    """
    rxn_names = [line.split()[0] for line in n_lines if line.strip()]
    counter = Counter(rxn_names)
    entries = sorted(counter.items(), key=lambda x: (x[1], locale.strxfrm(x[0])))
    lines = [f"{count:{UNIQ_C_WIDTH}d} {name}" for name, count in entries]
    return "\n".join(lines) + "\n" if lines else ""


def make_histogram(count_text: str) -> str:
    """Build a histogram of frequency values from ``uniq -c`` output.

    Replaces ``sed 's/ *//; s/ .*//' … | sort -n | uniq -c``.
    Takes the count column from ``count_per_reaction`` output and produces
    a histogram of how many reactions share each count value.

    Returns:
        String in ``uniq -c`` format.
    """
    counts = []
    for line in count_text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts:
            try:
                counts.append(int(parts[0]))
            except ValueError:
                continue
    if not counts:
        return ""
    counter = Counter(counts)
    entries = sorted(counter.items(), key=lambda x: (x[0], x[1]))
    lines = [f"{freq:{UNIQ_C_WIDTH}d} {value}" for value, freq in entries]
    return "\n".join(lines) + "\n"


def extract_nitrogen_atoms(n_lines: List[str]) -> List[str]:
    """Extract individual atom identifiers from nitrogen mapping pairs.

    Replaces ``sed 's/.* //; s/=/\\n/' all_mapping.N.sorted.txt``.
    Each line ``"<rxn> <from_atom>=<to_atom>"`` yields two atoms.
    """
    atoms: List[str] = []
    for line in n_lines:
        if not line.strip():
            continue
        pair = line.split()[-1]
        atoms.extend(pair.split("="))
    return atoms


def count_atoms(atoms: List[str], unique: bool = False) -> str:
    """Count atom occurrences, formatted like ``sort [-u] | uniq -c | sort -n``.

    Replaces the pipeline for both ``all_atoms.N.sorted.txt`` (unique)
    and ``all_atoms.N.count.txt`` (with counts).

    Args:
        atoms: List of atom identifier strings.
        unique: If True, deduplicate before counting (produces all 1s).

    Returns:
        String in ``uniq -c`` format.
    """
    if unique:
        sorted_atoms = sorted(set(atoms), key=locale.strxfrm)
        lines = [a for a in sorted_atoms]
        return "\n".join(lines) + "\n" if lines else ""
    counter = Counter(sorted(atoms, key=locale.strxfrm))
    entries = sorted(counter.items(), key=lambda x: (x[1], locale.strxfrm(x[0])))
    result_lines = [f"{count:{UNIQ_C_WIDTH}d} {name}" for name, count in entries]
    return "\n".join(result_lines) + "\n" if result_lines else ""


def _resolve_reactions_dir(reactions_dir: Path) -> Path:
    """Return a directory path for *reactions_dir*, extracting from zip if needed."""
    reactions_dir = Path(reactions_dir)
    if reactions_dir.is_dir():
        return reactions_dir
    if reactions_dir.is_file() and reactions_dir.suffix == ".zip":
        tmp = tempfile.TemporaryDirectory(prefix="reaction_intermediates_")
        with zipfile.ZipFile(reactions_dir) as zf:
            zf.extractall(tmp.name)
        top = Path(tmp.name)
        entries = [p for p in top.iterdir() if p.is_dir()]
        if len(entries) == 1:
            subdir = entries[0]
            if subdir.name == reactions_dir.stem:
                tmp_name = tmp.name
                tmp._finalizer.detach()
                return subdir
        return top
    raise FileNotFoundError(reactions_dir)


def unite_mappings(reactions_dir: Path, output_dir: Path) -> None:
    """Run the full unite_mappings pipeline and write all output files.

    Replaces the complete ``unite_mappings.sh`` script.

    Args:
        reactions_dir: Path to ``reaction_intermediates/`` directory or ``.zip`` archive.
        output_dir: Directory where output files are written.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    resolved = _resolve_reactions_dir(reactions_dir)
    all_lines, _ = collect_all_mappings(resolved)

    all_mapping_path = output_dir / "all_mapping.txt"
    all_mapping_path.write_text("\n".join(all_lines) + "\n" if all_lines else "")

    sorted_lines = sort_all_mappings(all_lines)
    sorted_path = output_dir / "all_mapping.sorted.txt"
    sorted_path.write_text("\n".join(sorted_lines) + "\n" if sorted_lines else "")

    n_lines = filter_nitrogen_mappings(sorted_lines)
    n_path = output_dir / "all_mapping.N.sorted.txt"
    n_path.write_text("\n".join(n_lines) + "\n" if n_lines else "")

    rxn_count_text = count_per_reaction(n_lines)
    (output_dir / "all_rxn_N_count.txt").write_text(rxn_count_text)

    rxn_histo_text = make_histogram(rxn_count_text)
    (output_dir / "all_rxn_N_count.histo").write_text(rxn_histo_text)

    atoms = extract_nitrogen_atoms(n_lines)

    atoms_sorted_text = count_atoms(atoms, unique=True)
    (output_dir / "all_atoms.N.sorted.txt").write_text(atoms_sorted_text)

    atoms_count_text = count_atoms(atoms, unique=False)
    (output_dir / "all_atoms.N.count.txt").write_text(atoms_count_text)

    atoms_histo_text = make_histogram(atoms_count_text)
    (output_dir / "all_atoms.N.count.histo").write_text(atoms_histo_text)


def main():
    """CLI entry point: run the unite_mappings pipeline."""
    _script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Combine per-reaction atom mappings into global tables"
    )
    parser.add_argument(
        "--reactions-dir",
        type=Path,
        default=_script_dir / "reaction_intermediates.zip",
        help="Directory or .zip archive containing reaction subfolders",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory where output files are written",
    )
    args = parser.parse_args()
    unite_mappings(args.reactions_dir, args.output_dir)


if __name__ == "__main__":
    main()
