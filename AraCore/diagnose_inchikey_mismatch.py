#!/usr/bin/env python3
"""
Comprehensive InChIKey diagnostic for AraCore.

Tests all combinations of InChIKey generation methods and reports
match rates against reference files.

Methods tested:
  A: SMILES → InChIKey (current recreate_data.sh)
  B: SMILES → InChIKey with -xT/nochg
  C: Existing MDL → InChIKey (keep M CHG)
  D: Existing MDL → InChIKey with -xT/nochg (keep M CHG)
  E: Existing MDL stripped → InChIKey (remove M CHG)
  F: Existing MDL stripped → InChIKey with -xT/nochg (remove M CHG)
  G: Fresh MDL from SMILES → InChIKey
  H: Fresh MDL from SMILES → InChIKey with -xT/nochg
"""

import csv
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set


@dataclass
class SpeciesResult:
    """Results for a single species across all methods."""
    species_id: str
    smiles: str
    reference_key: str
    reference_source: str  # Which reaction folder this came from
    generated_keys: Dict[str, str] = field(default_factory=dict)  # method -> inchikey
    match_status: Dict[str, str] = field(default_factory=dict)  # method -> FULL|PREFIX|NONE|SUFFIX_DIFF


@dataclass
class MethodStats:
    """Statistics for a single method."""
    total_lookups: int = 0
    full_matches: int = 0
    prefix_matches: int = 0
    no_match: int = 0
    suffix_diff: int = 0


def run_obabel(args: List[str]) -> str:
    """Run obabel command and return stdout."""
    cmd = ["obabel"] + args
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30,
                              encoding='utf-8', errors='ignore')
        if result.returncode != 0:
            return ""
        return result.stdout.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return ""


def generate_from_smiles(smiles: str, use_nochg: bool = False) -> Optional[str]:
    """Generate InChIKey from SMILES."""
    args = [f"-:{smiles}", "-oinchikey"]
    if use_nochg:
        args.append("-xT/nochg")
    # Add - to suppress extra output
    args.append("-")
    output = run_obabel(args)
    # InChIKey is the first (and usually only) line
    if output:
        return output.split('\n')[0].strip()
    return None


def generate_mdl_from_smiles(smiles: str, tmp_dir: Path) -> Optional[Path]:
    """Generate MDL file from SMILES, return path to generated file."""
    mdl_path = tmp_dir / "temp.mdl"
    args = [f"-:{smiles}", "-omdl", "-O", str(mdl_path)]
    output = run_obabel(args)
    if mdl_path.exists():
        return mdl_path
    return None


def generate_from_mdl(mdl_path: Path, use_nochg: bool = False, strip_chg: bool = False) -> Optional[str]:
    """Generate InChIKey from MDL file."""
    if strip_chg:
        # Create temp file with M CHG lines removed
        temp_mdl = mdl_path.parent / f"{mdl_path.stem}_nochg.mdl"
        with open(mdl_path, 'r') as f_in, open(temp_mdl, 'w') as f_out:
            for line in f_in:
                if not line.startswith('M  CHG'):
                    f_out.write(line)
        mdl_to_use = temp_mdl
    else:
        mdl_to_use = mdl_path

    args = ["-imdl", str(mdl_to_use), "-oinchikey", "-"]
    if use_nochg:
        args.append("-xT/nochg")

    output = run_obabel(args)
    if output:
        return output.split('\n')[0].strip()
    return None


def load_species_smiles(aracore_dir: Path) -> Dict[str, str]:
    """Load species_id_smiles.txt."""
    smiles_file = aracore_dir / "species_id_smiles.txt"
    species_smiles = {}
    with open(smiles_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                species_id = parts[0].strip()
                smiles = parts[1].strip()
                species_smiles[species_id] = smiles
    return species_smiles


def load_reference_inchikeys(aracore_dir: Path) -> Dict[str, Tuple[str, str]]:
    """
    Load all species_id_inchikey.txt files from reaction folders.

    Returns:
        Dict mapping species_id -> (inchikey, source_folder)
    """
    reference_keys = {}
    reactions_dir = aracore_dir / "reaction_intermediates"

    if not reactions_dir.exists():
        return reference_keys

    for reaction_folder in reactions_dir.iterdir():
        if not reaction_folder.is_dir():
            continue

        inchikey_file = reaction_folder / "species_id_inchikey.txt"
        if not inchikey_file.exists():
            continue

        with open(inchikey_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    species_id = parts[0].strip()
                    inchikey = parts[1].strip()
                    # Store first occurrence
                    if species_id not in reference_keys:
                        reference_keys[species_id] = (inchikey, reaction_folder.name)

    return reference_keys


def classify_match(generated: Optional[str], reference: str) -> str:
    """Classify the match type between generated and reference InChIKeys."""
    if not generated:
        return "NONE"

    if generated == reference:
        return "FULL"

    # Check if only suffix differs (last character)
    if len(generated) == 27 and len(reference) == 27:
        if generated[:14] == reference[:14]:
            if generated[14:-1] == reference[14:-1]:
                # Same connectivity and stereo/charge layers, only last char differs
                return "SUFFIX_DIFF"
            return "PREFIX"
        if generated[:14] == reference[:14]:
            return "PREFIX"

    return "NONE"


def find_existing_mdl_files(aracore_dir: Path) -> Dict[str, Path]:
    """Find existing MDL files for each species."""
    mdl_files = {}
    reactions_dir = aracore_dir / "reaction_intermediates"

    if not reactions_dir.exists():
        return mdl_files

    for reaction_folder in reactions_dir.iterdir():
        if not reaction_folder.is_dir():
            continue

        # Look for species_id file to map MOL files to species
        species_file = reaction_folder / "species_id_inchikey.txt"
        if not species_file.exists():
            continue

        # Load species mapping
        species_mapping = {}
        with open(species_file, 'r') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 1:
                    species_id = parts[0].strip()
                    mol_num = i  # MOL files are in order
                    species_mapping[mol_num] = species_id

        # Find MDL files
        for mol_file in reaction_folder.glob("MOL_??.mdl"):
            # Extract number from filename
            try:
                num = int(mol_file.stem.split('_')[1])
                if num in species_mapping:
                    species_id = species_mapping[num]
                    # Store first found
                    if species_id not in mdl_files:
                        mdl_files[species_id] = mol_file
            except (ValueError, IndexError):
                continue

    return mdl_files


def test_aracore(aracore_dir: Path, tmp_dir: Path) -> Tuple[Dict[str, SpeciesResult], Dict[str, MethodStats]]:
    """Run all diagnostic tests on AraCore dataset."""
    print("Loading AraCore data...", file=sys.stderr)

    species_smiles = load_species_smiles(aracore_dir)
    reference_keys = load_reference_inchikeys(aracore_dir)
    existing_mdls = find_existing_mdl_files(aracore_dir)

    print(f"  - {len(species_smiles)} species with SMILES")
    print(f"  - {len(reference_keys)} species with reference InChIKeys")
    print(f"  - {len(existing_mdls)} species with existing MDL files")

    results = {}
    stats = {f"Method_{chr(65+i)}": MethodStats() for i in range(8)}  # A-H

    # Track which species we've tested
    tested_species = set()

    # Test species that have reference keys
    for species_id, (reference_key, source) in reference_keys.items():
        if species_id not in species_smiles:
            continue

        smiles = species_smiles[species_id]
        result = SpeciesResult(
            species_id=species_id,
            smiles=smiles,
            reference_key=reference_key,
            reference_source=source
        )

        # Method A: SMILES → InChIKey (current)
        key_a = generate_from_smiles(smiles, use_nochg=False)
        result.generated_keys["A"] = key_a or ""
        result.match_status["A"] = classify_match(key_a, reference_key)
        stats["Method_A"].total_lookups += 1

        # Method B: SMILES → InChIKey with nochg
        key_b = generate_from_smiles(smiles, use_nochg=True)
        result.generated_keys["B"] = key_b or ""
        result.match_status["B"] = classify_match(key_b, reference_key)
        stats["Method_B"].total_lookups += 1

        # Methods C-F: Existing MDL files
        if species_id in existing_mdls:
            mdl_path = existing_mdls[species_id]

            # Method C: MDL as-is
            key_c = generate_from_mdl(mdl_path, use_nochg=False, strip_chg=False)
            result.generated_keys["C"] = key_c or ""
            result.match_status["C"] = classify_match(key_c, reference_key)
            stats["Method_C"].total_lookups += 1

            # Method D: MDL as-is + nochg
            key_d = generate_from_mdl(mdl_path, use_nochg=True, strip_chg=False)
            result.generated_keys["D"] = key_d or ""
            result.match_status["D"] = classify_match(key_d, reference_key)
            stats["Method_D"].total_lookups += 1

            # Method E: MDL stripped
            key_e = generate_from_mdl(mdl_path, use_nochg=False, strip_chg=True)
            result.generated_keys["E"] = key_e or ""
            result.match_status["E"] = classify_match(key_e, reference_key)
            stats["Method_E"].total_lookups += 1

            # Method F: MDL stripped + nochg
            key_f = generate_from_mdl(mdl_path, use_nochg=True, strip_chg=True)
            result.generated_keys["F"] = key_f or ""
            result.match_status["F"] = classify_match(key_f, reference_key)
            stats["Method_F"].total_lookups += 1

        # Methods G-H: Fresh MDL from SMILES
        fresh_mdl = generate_mdl_from_smiles(smiles, tmp_dir)
        if fresh_mdl:
            # Method G: Fresh MDL
            key_g = generate_from_mdl(fresh_mdl, use_nochg=False, strip_chg=False)
            result.generated_keys["G"] = key_g or ""
            result.match_status["G"] = classify_match(key_g, reference_key)
            stats["Method_G"].total_lookups += 1

            # Method H: Fresh MDL + nochg
            key_h = generate_from_mdl(fresh_mdl, use_nochg=True, strip_chg=False)
            result.generated_keys["H"] = key_h or ""
            result.match_status["H"] = classify_match(key_h, reference_key)
            stats["Method_H"].total_lookups += 1

            # Clean up
            fresh_mdl.unlink()

        results[species_id] = result
        tested_species.add(species_id)

    # Calculate statistics
    for result in results.values():
        for method, status in result.match_status.items():
            method_key = f"Method_{method}"
            if method_key not in stats:
                continue
            if status == "FULL":
                stats[method_key].full_matches += 1
            elif status == "PREFIX":
                stats[method_key].prefix_matches += 1
            elif status == "SUFFIX_DIFF":
                stats[method_key].suffix_diff += 1
            else:
                stats[method_key].no_match += 1

    return results, stats


def test_metacyc_sample(metacyc_dir: Path, tmp_dir: Path, sample_size: int = 100) -> Tuple[Dict[str, SpeciesResult], Dict[str, MethodStats]]:
    """Run diagnostic tests on a sample of MetaCyc reactions."""
    print(f"\nTesting MetaCyc sample (n={sample_size})...", file=sys.stderr)

    reactions_dir = metacyc_dir / "reaction_intermediates"
    if not reactions_dir.exists():
        print("  MetaCyc reaction_intermediates not found, skipping", file=sys.stderr)
        return {}, {}

    # Sample reaction folders
    reaction_folders = [f for f in reactions_dir.iterdir() if f.is_dir()]
    if len(reaction_folders) > sample_size:
        import random
        random.shuffle(reaction_folders)
        reaction_folders = reaction_folders[:sample_size]

    print(f"  - Sampled {len(reaction_folders)} reaction folders")

    # Collect species from sampled reactions
    species_smiles = {}
    reference_keys = {}

    for reaction_folder in reaction_folders:
        smiles_file = reaction_folder / "species_id_smiles.txt"
        inchikey_file = reaction_folder / "species_inchikey.txt"

        if smiles_file.exists():
            with open(smiles_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('Object'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        species_id = parts[0].strip()
                        smiles = parts[1].strip()
                        if species_id not in species_smiles:
                            species_smiles[species_id] = smiles

        if inchikey_file.exists():
            with open(inchikey_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('Object'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        species_id = parts[0].strip()
                        inchikey = parts[1].strip()
                        if species_id not in reference_keys:
                            reference_keys[species_id] = (inchikey, reaction_folder.name)

    print(f"  - {len(reference_keys)} species with reference InChIKeys")

    results = {}
    stats = {f"MetaCyc_Method_{chr(65+i)}": MethodStats() for i in range(2)}  # Only A, B for MetaCyc

    for species_id, (reference_key, source) in list(reference_keys.items())[:sample_size]:
        if species_id not in species_smiles:
            continue

        smiles = species_smiles[species_id]
        result = SpeciesResult(
            species_id=species_id,
            smiles=smiles,
            reference_key=reference_key,
            reference_source=source
        )

        # Only test Methods A and B for MetaCyc (no MDL files in same structure)
        key_a = generate_from_smiles(smiles, use_nochg=False)
        result.generated_keys["A"] = key_a or ""
        result.match_status["A"] = classify_match(key_a, reference_key)
        stats["MetaCyc_Method_A"].total_lookups += 1

        key_b = generate_from_smiles(smiles, use_nochg=True)
        result.generated_keys["B"] = key_b or ""
        result.match_status["B"] = classify_match(key_b, reference_key)
        stats["MetaCyc_Method_B"].total_lookups += 1

        results[species_id] = result

    # Calculate statistics
    for result in results.values():
        for method in ["A", "B"]:
            status = result.match_status[method]
            method_key = f"MetaCyc_Method_{method}"
            if status == "FULL":
                stats[method_key].full_matches += 1
            elif status == "PREFIX":
                stats[method_key].prefix_matches += 1
            elif status == "SUFFIX_DIFF":
                stats[method_key].suffix_diff += 1
            else:
                stats[method_key].no_match += 1

    return results, stats


def print_report(aracore_stats: Dict[str, MethodStats],
                 metacyc_stats: Dict[str, MethodStats],
                 aracore_results: Dict[str, SpeciesResult]) -> None:
    """Print formatted report to console."""
    print("=" * 80)
    print("InChIKey Diagnostic Report")
    print("=" * 80)
    print()

    # AraCore Results
    print("AraCore Results:")
    print("-" * 80)
    print(f"{'Method':<20} {'Total':>8} {'Full':>8} {'Prefix':>8} {'Suffix':>8} {'None':>8} {'Full%':>8}")
    print("-" * 80)

    method_descriptions = {
        "Method_A": "A: SMILES (current)",
        "Method_B": "B: SMILES + nochg",
        "Method_C": "C: MDL existing",
        "Method_D": "D: MDL existing + nochg",
        "Method_E": "E: MDL stripped",
        "Method_F": "F: MDL stripped + nochg",
        "Method_G": "G: MDL fresh",
        "Method_H": "H: MDL fresh + nochg",
    }

    for method_key in ["Method_A", "Method_B", "Method_E", "Method_F"]:
        if method_key not in aracore_stats:
            continue
        stats = aracore_stats[method_key]
        description = method_descriptions.get(method_key, method_key)
        total = stats.total_lookups
        full_pct = (stats.full_matches / total * 100) if total > 0 else 0

        print(f"{description:<20} {total:>8} {stats.full_matches:>8} "
              f"{stats.prefix_matches:>8} {stats.suffix_diff:>8} {stats.no_match:>8} {full_pct:>7.1f}%")

    print()

    # MetaCyc Results
    if metacyc_stats:
        print("\nMetaCyc Sample Results:")
        print("-" * 80)
        print(f"{'Method':<20} {'Total':>8} {'Full':>8} {'Prefix':>8} {'Suffix':>8} {'None':>8} {'Full%':>8}")
        print("-" * 80)

        for method_key in ["MetaCyc_Method_A", "MetaCyc_Method_B"]:
            if method_key not in metacyc_stats:
                continue
            stats = metacyc_stats[method_key]
            description = method_key.replace("MetaCyc_", "")
            total = stats.total_lookups
            full_pct = (stats.full_matches / total * 100) if total > 0 else 0

            print(f"{description:<20} {total:>8} {stats.full_matches:>8} "
                  f"{stats.prefix_matches:>8} {stats.suffix_diff:>8} {stats.no_match:>8} {full_pct:>7.1f}%")
        print()

    # Top mismatches
    print("\nTop 10 Species with Mismatches (Method A - current):")
    print("-" * 80)

    mismatches = [(s_id, r) for s_id, r in aracore_results.items()
                  if r.match_status.get("A") != "FULL"]
    mismatches.sort(key=lambda x: x[0])

    for i, (species_id, result) in enumerate(mismatches[:10], 1):
        status_a = result.match_status.get("A", "NONE")
        key_a = result.generated_keys.get("A", "")
        key_b = result.generated_keys.get("B", "")
        ref = result.reference_key

        print(f"{i}. {species_id}")
        print(f"   Status: {status_a}")
        print(f"   Reference: {ref}")
        print(f"   Method A:  {key_a}")
        print(f"   Method B:  {key_b}")
        print()

    # Recommendations
    print("\nRecommendations:")
    print("-" * 80)

    # Find best method
    best_method = None
    best_rate = 0
    for method_key, stats in aracore_stats.items():
        if stats.total_lookups == 0:
            continue
        rate = stats.full_matches / stats.total_lookups
        if rate > best_rate:
            best_rate = rate
            best_method = method_key

    if best_method:
        description = method_descriptions.get(best_method, best_method)
        print(f"Best match rate: {description} ({best_rate*100:.1f}% full matches)")

    # Check if nochg helps
    if "Method_A" in aracore_stats and "Method_B" in aracore_stats:
        rate_a = aracore_stats["Method_A"].full_matches / aracore_stats["Method_A"].total_lookups
        rate_b = aracore_stats["Method_B"].full_matches / aracore_stats["Method_B"].total_lookups
        if rate_b > rate_a:
            print(f"  → Adding -xT/nochg improves match rate by {(rate_b-rate_a)*100:.1f}%")
        elif rate_a > rate_b:
            print(f"  → Current method (without -xT/nochg) is better")

    # Check suffix differences
    suffix_diff_count = sum(1 for r in aracore_results.values()
                            if r.match_status.get("A") == "SUFFIX_DIFF")
    if suffix_diff_count > 0:
        print(f"\n  → {suffix_diff_count} species have suffix differences (-M vs -N)")
        print(f"     These are charged molecules affected by charge stripping")

    # Fallback recommendation
    if aracore_stats["Method_A"].no_match > 0:
        total = aracore_stats["Method_A"].total_lookups
        with_prefix = aracore_stats["Method_A"].full_matches + aracore_stats["Method_A"].prefix_matches
        print(f"\n  → Adding 14-char prefix fallback would recover {aracore_stats['Method_A'].prefix_matches} matches")
        print(f"     Total coverage: {with_prefix}/{total} ({with_prefix/total*100:.1f}%)")

    print()
    print("=" * 80)


def export_csv(aracore_results: Dict[str, SpeciesResult],
               aracore_stats: Dict[str, MethodStats],
               metacyc_stats: Dict[str, MethodStats],
               output_path: Path) -> None:
    """Export detailed results to CSV."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)

        # Header
        header = ['species_id', 'reference_inchikey', 'smiles']
        for method in ['A', 'B', 'E', 'F']:
            header.extend([f'method_{method}_key', f'method_{method}_match'])
        writer.writerow(header)

        # Data rows
        for species_id in sorted(aracore_results.keys()):
            result = aracore_results[species_id]
            row = [
                species_id,
                result.reference_key,
                result.smiles
            ]
            for method in ['A', 'B', 'E', 'F']:
                row.append(result.generated_keys.get(method, ""))
                row.append(result.match_status.get(method, "NONE"))
            writer.writerow(row)

    print(f"\nCSV exported to: {output_path}")


def main():
    import tempfile

    # Determine paths
    script_dir = Path(__file__).parent
    aracore_dir = script_dir
    metacyc_dir = script_dir.parent / "MetaCyc"

    # Create temp directory for MDL generation
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Run AraCore tests
        aracore_results, aracore_stats = test_aracore(aracore_dir, tmp_path)

        # Run MetaCyc sample tests
        metacyc_results, metacyc_stats = test_metacyc_sample(metacyc_dir, tmp_path)

        # Print report
        print_report(aracore_stats, metacyc_stats, aracore_results)

        # Export CSV
        csv_path = script_dir / "inchikey_diagnostic_results.csv"
        export_csv(aracore_results, aracore_stats, metacyc_stats, csv_path)


if __name__ == "__main__":
    main()
