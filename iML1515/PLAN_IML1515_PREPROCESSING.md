# Plan: iML1515 Preprocessing Script

## Overview

Create `iML1515/preprocess_iml1515.py` -- a standalone Python script that reads the
iML1515 genome-scale metabolic model (JSON format) and generates the per-reaction
input folders that `run_rdt.py` expects.  The script is purely a preprocessor;
the actual RDT atom mapping is handled by the existing `AraCore/run_rdt.py`.

---

## 1.  Script Location & Output Directory

| Decision | Value |
|----------|-------|
| Script path | `iML1515/preprocess_iml1515.py` |
| Model input | `/home/arne/repos/atnlib/data/BiGG/iML1515.json.gz` (configurable via CLI) |
| Output | `iML1515/reaction_intermediates/<rxn_id>/` folder per reaction (configurable) |
| Integration with `run_rdt.py` | Import from `AraCore/run_rdt.py` (no duplication). The preprocessing step and the RDT step are run separately. |

---

## 2.  CLI Interface

```
python iML1515/preprocess_iml1515.py \
    --model-path /path/to/iML1515.json.gz \
    --output-dir iML1515/reaction_intermediates \
    --cache-file iML1515/smiles_cache.json
```

| Flag | Default | Description |
|------|---------|-------------|
| `--model-path` | hardcoded known path | Path to the gzipped iML1515 JSON model |
| `--output-dir` | `iML1515/reaction_intermediates` | Output directory for per-reaction folders |
| `--cache-file` | `iML1515/smiles_cache.json` | JSON cache for SMILES/InChIKey API results |
| `--verbose` / `-v` | False | Increment verbosity (repeatable) |

---

## 3.  SMILES Source Strategy

Priority chain (stop at first success, with provenance tracking):

1. **ModelSEED Solr API** -- batch-query by `seed.compound` ID
   - Endpoint: `https://modelseed.org/solr/compounds/select?q=id:cpdXXXXX&fl=id,smiles,inchikey&wt=json`
   - Returns both SMILES and InChIKey
   - Covers 1602/1877 metabolites
   - Tested: works, no auth needed

2. **PubChem REST API** -- query by `inchi_key` annotation
   - Endpoint: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/<key>/property/IsomericSMILES,InChIKey/JSON`
   - Returns SMILES and InChIKey
   - Covers 45 additional metabolites (0 unique beyond ChEBI below)
   - Tested: works, no auth needed

3. **ChEBI REST API** -- query by `chebi` annotation
   - Endpoint: `https://www.ebi.ac.uk/chebi/backend/api/public/compound/<id>/`
   - Returns SMILES and InChIKey in `default_structure`
   - Batch endpoint also available: `/chebi/backend/api/public/compounds/?chebi_ids=CHEBI:1,CHEBI:2`
   - Covers 46 metabolites (1 unique beyond PubChem: `moco_c`)
   - Tested: works, no authentication needed, free

4. **Give up** -- 229 metabolites remain (lipids, protein-bound clusters, generic residues).
   Reactions containing these are skipped.

### Caching

All fetched SMILES + InChIKeys are saved to a JSON cache file. On subsequent runs,
the cache is loaded first and only uncached IDs are fetched from APIs.

---

## 4.  InChIKey Source (with Provenance)

Priority chain:

1. **Model annotation** -- use `inchi_key` field from iML1515 JSON if present
   (covers 1507/1877 metabolites)
2. **API response** -- use InChIKey returned by ModelSEED / PubChem / ChEBI
3. **obabel fallback** -- if an API provided SMILES but no InChIKey, compute via
   `obabel -:"<smiles>" -oinchikey`

Each InChIKey is tagged with its source in the provenance log.

---

## 5.  Metabolites Without SMILES

Reactions where **any** metabolite lacks SMILES are **skipped entirely**.
This includes:
- 50 reactions with zero SMILES coverage (all metabolites are "stuck")
- 376 reactions with partial SMILES coverage (some metabolites are "stuck")

Expected yield: ~2278 reactions from ModelSEED alone, ~87 more recovered via
fallback APIs = **~2365+ mappable reactions** out of 2712 total.

---

## 6.  Reaction Filtering

Reactions are excluded before output:

| Category | Count | Reason |
|----------|-------|--------|
| Exchange (prefix `EX_` or single-met boundary) | 337 | No meaningful atom mapping |
| Biomass (`BIOMASS` in ID) | 2 | Pseudoreaction, 70-99 metabolites |
| Any metabolite lacks SMILES | 434 | Cannot build complete SMILES string |

Transport reactions (749) and standard metabolic reactions (1574) are **included**,
matching AraCore's approach. RDT supports transport with the `-b` flag.

---

## 7.  Species ID Format Conversion

iML1515 uses underscore-delimited compartment notation (`dhap_c`, `gtp_e`).
`run_rdt.py` **requires** bracket-enclosed compartment notation (`dhap[c]`, `gtp[e]`)
because `find_species_with_cmp` matches on the `name[` prefix.

Conversion rule: `dhap_c` -> `dhap[c]` (split on last underscore, wrap in brackets)

iML1515 compartments: `c` (cytosol), `e` (extracellular), `p` (periplasm)

AraCore's `_DASH_` convention for dashes in IDs is **not needed** -- iML1515 has
no dashes in any base name. Double underscores (`cysi__L_c`) encode stereochemistry
and are preserved as-is in the base name.

---

## 8.  Stoichiometry Flattening

Non-unit stoichiometry (248 reactions) is flattened by repeating the species,
matching AraCore's `rxn_table_with_cmp_flat.txt` approach:

```
PPA: h2o_c(-1), pi_c(2.0), h_c(1.0), ppi_c(-1.0)
-> flattened: h2o[c] + pi[c] + pi[c] >> h[c] + ppi[c]
```

This is necessary because RDT needs each molecule instance as a separate entry
in the MDL output to produce correct atom mappings.

---

## 9.  Reaction Direction

All reactions use `>>` as the SMILES separator, regardless of reversibility:

```
reactant1.SMILES.reactant2.SMILES>>product1.SMILES.product2.SMILES
```

Reversible (`<=>`) and irreversible (`-->`, `<==`) reactions are all treated the
same way. This matches AraCore's approach.

The iML1515 JSON stores reactions with `lower_bound` and `upper_bound` fields;
reactants are metabolites with negative coefficients, products have positive
coefficients.

---

## 10.  Output Format

Each reaction folder (`reaction_intermediates/<rxn_id>/`) contains 4 files,
exactly matching what `AraCore/run_rdt.py` expects:

### `rxn.smiles`
Single line: SMILES of reactants, `.`-joined, then `>>`, then SMILES of products,
`.`-joined.

Example:
```
O.CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-].CC(C)C(=O)C(=O)[O-]>>[H+].CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS.CC(C)[C@@](O)(CC(=O)[O-])C(=O)[O-]
```

### `from_species_with_cmp`
Sorted, deduplicated list of reactant species IDs with bracket compartments,
one per line (`sort -u` order).

Example:
```
dhap[c]
gtp[c]
h2o[c]
```

### `to_species_with_cmp`
Same format for product species IDs.

### `species_id_inchikey.txt`
Tab-separated, one line per unique species (across both sides):
`<base_name>\t<inchikey>`
The base name is the species ID **without** compartment suffix.

Example:
```
dhap    GNGACRATGGDKBX-UHFFFAOYSA-L
gtp ZKHQWZAMYRWXGA-KQYNXXCUSA-K
h2o XLYOFNOQVPJJNP-UHFFFAOYSA-N
```

---

## 11.  Logging

- **Stdout**: progress summary -- number of metabolites fetched per source,
  number of reactions included/skipped per category, final counts
- **Stderr**: warnings and errors
- **Cache file JSON**: contains all fetched SMILES/InChIKeys with provenance
  metadata (source API, timestamp)

---

## 12.  Expected Impact

| Metric | Count |
|--------|-------|
| Total iML1515 reactions | 2,712 |
| Excluded (exchange + biomass + no-SMILES) | ~347 |
| Expected mappable reactions | ~2,365 |
| Metabolites covered (ModelSEED) | 1,602 |
| Additional via PubChem | 45 |
| Additional via ChEBI | 1 |
| Uncovered (unrecoverable) | 229 |

---

## 13.  Usage Workflow

```bash
# Step 1: Preprocess
python iML1515/preprocess_iml1515.py

# Step 2: Run RDT
python -c "
import sys; sys.path.insert(0, 'AraCore')
from run_rdt import main as rdt_main
import sys; sys.argv = ['run_rdt.py', '--reactions-dir', 'iML1515/reaction_intermediates']
rdt_main()
"

# Step 3: Post-process (process_reaction_data / mapping assembly)
# Same as Step 2 -- run_rdt.py handles both RDT + postprocessing
```
