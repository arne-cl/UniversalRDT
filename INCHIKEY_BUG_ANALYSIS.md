# InChIKey Mismatch Bug Analysis

## Problem Statement

The `run_rdt.sh` script hangs when processing certain reactions because `grep` receives an empty pattern and waits for stdin. The commit `bee15e6` added guards to prevent the hang, but this only hides the symptom, not the root cause.

## Root Cause (Identified 2026-04-16)

The InChIKey lookup fails due to **inconsistent charge handling** between reference and generated InChIKeys.

### Evidence

| Method | Input | InChIKey |
|--------|-------|----------|
| Reference (`recreate_data.sh`) | Charged SMILES `O=C([O-])c1cc(=O)[nH]c(=O)[nH]1` | `PXQPEWDEAKTCGB-UHFFFAOYSA-M` |
| Generated (`run_rdt.sh`) | Charge-stripped MDL | `PXQPEWDEAKTCGB-UHFFFAOYSA-N` |

**The suffix differs** (`-M` vs `-N`) because:
- `-M` = charge layer included
- `-N` = no charge layer (or charge ignored)

### Code Flow

**1. Reference InChIKey generation** (`AraCore/recreate_data.sh`):
```bash
obabel -:"$smiles" -oinchikey  # Includes charge info
```

**2. RDT processing** (`AraCore/run_rdt.sh` line 24):
```bash
tail -n +2 $fn2 | grep -v '^M  CHG' > ${fn2}.mdl  # REMOVES charge info!
obabel -i mdl ${fn2}.mdl -oinchikey -O ${fn2}.inchikey  # No charge info
```

**3. Lookup fails** (`AraCore/run_rdt.sh` line 36):
```bash
species_id_without_cmp=$(grep "$(cat ${fn2}.inchikey)" species_id_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
# InChIKey not found → species_id_without_cmp is empty
```

### Why the Script Hangs

When `species_id_without_cmp` is empty and unquoted:
```bash
grep $species_id_without_cmp from_species_with_cmp
# Expands to: grep from_species_with_cmp
# grep treats "from_species_with_cmp" as a pattern and reads from stdin → HANGS
```

## Open Babel Version Context

| Version | Release Date | Relevant Changes |
|---------|--------------|------------------|
| 3.0.0 | 2019 | "Code for handling implicit hydrogens and kekulization has been entirely replaced" |
| 3.1.0 | ~2020 | "Fixed tautomer code" (PR #2171) |
| 3.1.1 | March 2024 | Current version in use |

The paper was published in 2022, so original code likely used Open Babel 3.0.0 (2019).

### Version Testing (2026-04-16)

**Tested with Docker (Ubuntu 20.04, Open Babel 3.0.0):**

| Version | Charged SMILES | With `-xT/nochg` | Neutral SMILES |
|---------|----------------|------------------|----------------|
| 3.0.0 (Ubuntu 20.04) | `PXQPEWDEAKTCGB-UHFFFAOYSA-M` | `PXQPEWDEAKTCGB-UHFFFAOYSA-N` | `PXQPEWDEAKTCGB-UHFFFAOYSA-N` |
| 3.1.1 (current) | `PXQPEWDEAKTCGB-UHFFFAOYSA-M` | `PXQPEWDEAKTCGB-UHFFFAOYSA-N` | `PXQPEWDEAKTCGB-UHFFFAOYSA-N` |

**Result**: Identical behavior. The bug is **NOT caused by Open Babel version differences**.

See `test-obabel-300.dockerfile` for the test Dockerfile.

## Test Case

Run `./test_inchikey_mismatch.sh` to reproduce the bug:
```bash
$ ./test_inchikey_mismatch.sh
*** BUG CONFIRMED ***
The grep for InChIKey will fail because:
  - species_id_inchikey.txt contains: PXQPEWDEAKTCGB-UHFFFAOYSA-M
  - Generated from charge-stripped MDL: PXQPEWDEAKTCGB-UHFFFAOYSA-N
```

## Open Questions

1. **Why does `run_rdt.sh` remove `M  CHG` lines?**
   - This is intentional (`grep -v '^M  CHG'`)
   - **Answer from paper**: The authors used `-xT/nochg` because "RDT sometimes changes protons" and they want to ignore these changes for canonical atom identification.
   - However, this creates inconsistency: `-xT/nochg` is used for InChI (canonical ordering) but NOT for InChIKey (species identification)

2. **Does Open Babel 3.0.0 behave differently?**
   - **Answer**: No, identical behavior confirmed via Docker test

## Comprehensive Diagnostic Results (2026-04-16)

Run `python3 AraCore/diagnose_inchikey_mismatch.py` for full diagnostic.

### AraCore Results (225 species tested)

| Method | Description | Total | Full Matches | Suffix Diff | No Match | Full% |
|--------|-------------|-------|--------------|-------------|----------|-------|
| A | SMILES (current recreate_data.sh) | 225 | 225 | 0 | 0 | **100.0%** |
| B | SMILES + `-xT/nochg` | 225 | 57 | 165 | 3 | 25.3% |
| E | MDL stripped (no M CHG) | 223 | 74 | 1 | 148 | 33.2% |
| F | MDL stripped + `-xT/nochg` | 223 | 23 | 1 | 50 | 10.3% |

### Key Findings

1. **Reference files are correct**: `species_id_inchikey.txt` files contain charged InChIKeys (-M suffix)
2. **Current method works perfectly**: Method A (SMILES without `-xT/nochg`) matches 100%
3. **165 species affected by `-xT/nochg`**: These charged molecules get -N suffix instead of -M
4. **Pre-existing MDL files work**: When `M  CHG` lines are present, obabel infers charge from atom records

### The Bug Location

The bug **does not occur** in `recreate_data.sh` or in the reference files. It occurs **during RDT processing** when:

1. RDT generates a new MDL file from reaction SMILES
2. `run_rdt.sh` removes `M  CHG` lines (line 24)
3. If RDT modified protonation, the atom records may have wrong charges
4. obabel generates InChIKey without proper charge info → mismatch

### Why Pre-existing Files Work

The pre-existing `.inchikey` files in `reaction_intermediates/` folders contain correct InChIKeys because:
- The original run preserved charge information correctly
- OR the original run used the correct method (matching reference)

## Comparison with MetaCyc Script

The MetaCyc script (`MetaCyc/run_rdt_metacyc.sh` lines 39-42) has a fallback:
```bash
if [ -z "$species_id" ]
then
species_id=$(grep $(head -c14 ${fn2}.inchikey) species_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
fi
```

This matches only the first 14 characters (connectivity layer) when the full InChIKey fails. The AraCore script lacks this fallback.

## References

- Open Babel releases: https://github.com/openbabel/openbabel/releases
- Tautomer fix PR #2171: commit `a40bfa690` (May 2020)
- InChI format: `-xT/nochg` ignores charge and protonation
