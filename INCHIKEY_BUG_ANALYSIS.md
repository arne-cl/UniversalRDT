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
   - Need to consult the paper to understand the rationale
   - Possible reasons: RDT may modify protonation incorrectly

2. **Does Open Babel 3.0.0 behave differently?**
   - Need to test if 3.0.0 handles the charge-stripped MDL differently
   - Or if the `-xT/nochg` flag behavior changed

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
