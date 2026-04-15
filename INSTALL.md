# UniversalRDT — Installation & Running Guide

This document describes how to install and run the UniversalRDT pipeline
using Docker for reproducibility. The pipeline computes **universal
atom-to-atom mappings** for metabolic reactions using the Reaction Decoder
Tool (RDT), then post-processes the output into InChI-ordered atom
mapping tables.

> **Important:** This guide covers the **bash-based atom mapping pipeline**
> only. Downstream analysis scripts (`createAraCoreNumbers.m`, `simMFA.m`,
> `create_S_N.m`) are written in MATLAB and require the COBRA Toolbox. See
> [Downstream MATLAB analysis](#downstream-matlab-analysis) for details.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Quick Start](#quick-start)
3. [Step-by-step: Download external dependencies](#step-by-step-download-external-dependencies)
4. [Step-by-step: Build the Docker image](#step-by-step-build-the-docker-image)
5. [Step-by-step: Prepare input data](#step-by-step-prepare-input-data)
6. [Step-by-step: Run the pipeline](#step-by-step-run-the-pipeline)
7. [Pipeline overview](#pipeline-overview)
8. [Key gotchas & hurdles](#key-gotchas--hurdles)
9. [Downstream MATLAB analysis](#downstream-matlab-analysis)

---

## Prerequisites

- **Docker** (or Podman) installed on the host
- **~130 MB** free disk for the RDT JAR download
- **~250 MB** free disk for the Docker image
- **~500 MB** free disk for extracted reaction intermediates (MetaCyc + AraCore)
- Internet access (to download the RDT JAR; everything else is in the repo)

---

## Quick Start

```bash
# 1. Download RDT JAR (v2.5.0 — the exact version used by this repo)
curl -L -o rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar \
    "https://github.com/asad/ReactionDecoder/releases/download/v2.5.0/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar"

# 2. Build Docker image
docker build -t universal-rdt .

# 3. Extract zipped data (run once)
#    MetaCyc: split zip parts → single zip → extract
cat MetaCyc/reaction_intermediates.zip.01 MetaCyc/reaction_intermediates.zip.02 \
    > MetaCyc/reaction_intermediates.zip
#    AraCore: single zip
#    (both are extracted inside the container below)

# 4. Run the container
docker run --rm -v "$(pwd):/data" -w /data universal-rdt bash -c "\
    unzip -o MetaCyc/reaction_intermediates.zip -d MetaCyc/ && \
    unzip -o MetaCyc/atom_mappings.zip -d MetaCyc/ && \
    unzip -o AraCore/reaction_intermediates.zip -d AraCore/ && \
    echo 'Data extracted.'"
```

---

## Step-by-step: Download external dependencies

### RDT JAR (Reaction Decoder Tool v2.5.0)

The pipeline scripts call RDT as a fat JAR. **You must use exactly version
2.5.0** — the output file naming (`ECBLAST_smiles_AAM.rxn`) and the
mapping format are specific to this version. Later versions (v3.x, v4.0)
changed the output format and are **not compatible**.

```bash
curl -L -o rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar \
    "https://github.com/asad/ReactionDecoder/releases/download/v2.5.0/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar"
```

Verify the download (~30 MB):

```bash
file rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar
# Should say: Java archive data (JAR)
```

No other external downloads are needed. All dataset files
(`species_id_smiles.txt`, `rxn_table.smiles.txt`, SBML model, etc.) are
already in the repository.

---

## Step-by-step: Build the Docker image

```bash
docker build -t universal-rdt .
```

The Dockerfile installs:

| Package | Purpose |
|---------|---------|
| `openjdk-21-jre-headless` | Java runtime for RDT |
| `openbabel` | `obabel` for SMILES ↔ InChI/InChIKey conversion |
| `unzip` | Extract zipped reaction intermediates and atom mappings |

The RDT JAR is placed at `/opt/rdt/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar`
and exposed as the `RDT_JAR` environment variable.

> **Note:** The original scripts reference RDT at the relative path
> `../../../../FluxAndPoolSizeEstimation/ReactionDecoder/target/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar`.
> Inside the container, we use `$RDT_JAR` instead. The run scripts in this
> repo must be invoked with the `RDT_JAR` variable set, or the `java -jar`
> path adjusted. See [Running the pipeline](#step-by-step-run-the-pipeline).

---

## Step-by-step: Prepare input data

Before running the pipeline, the zipped archives must be extracted.

### MetaCyc

The MetaCyc reaction intermediates are split into two zip parts (Git LFS
or large file handling). Reassemble and extract:

```bash
# On the host (outside the container):
cat MetaCyc/reaction_intermediates.zip.01 MetaCyc/reaction_intermediates.zip.02 \
    > MetaCyc/reaction_intermediates.zip

# Extract inside the container:
docker run --rm -v "$(pwd):/data" -w /data universal-rdt bash -c "\
    unzip -o MetaCyc/reaction_intermediates.zip -d MetaCyc/ && \
    unzip -o MetaCyc/atom_mappings.zip -d MetaCyc/"
```

This creates:
- `MetaCyc/reaction_intermediates/<RXN_ID>/` — one directory per reaction,
  each containing a `rxn.smiles` file
- `MetaCyc/atom_mappings/<RXN_ID>` — MetaCyc reference atom mappings

### AraCore

```bash
docker run --rm -v "$(pwd):/data" -w /data universal-rdt bash -c "\
    unzip -o AraCore/reaction_intermediates.zip -d AraCore/"
```

This creates `AraCore/reaction_intermediates/<RXN_ID>/` with:
- `rxn.smiles` — reaction in SMILES notation
- `from_species_with_cmp` — substrate species IDs (with compartment)
- `to_species_with_cmp` — product species IDs (with compartment)
- `species_id_inchikey.txt` — species-to-InChIKey lookup

---

## Step-by-step: Run the pipeline

The pipeline has three phases per dataset: **prepare**, **run RDT**, and
**postprocess**. Below are the commands for each dataset.

> **Important:** All `java -jar` calls in the original scripts use a
> hardcoded relative path. When running inside the container, set
> `RDT_JAR` so you can substitute. The command pattern becomes:
> ```
> java -jar $RDT_JAR -Q SMI -q "$smiles" -g -c -b -j AAM -f TEXT
> ```

### Option A: Run with adjusted environment variable

Start an interactive shell in the container:

```bash
docker run --rm -it -v "$(pwd):/data" -w /data universal-rdt bash
```

Then, inside the container, export the variable used by our wrapper calls:

```bash
export RDT_JAR=/opt/rdt/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar
```

### AraCore Pipeline

#### 1. Prepare (if starting from scratch, not using the pre-generated data)

The `generated/` directory already contains pre-generated files
(`rxn_table.smiles.txt`, `rxn_table_with_cmp_flat.txt`, etc.), so the
`prepare_rdt.sh` step may not be needed. If you need to re-generate:

```bash
cd /data/AraCore
bash prepare_rdt.sh
```

This script:
- Reads `species_id_smiles.txt` (species ID → SMILES mapping)
- Replaces compound IDs in `generated/rxn_table_with_cmp_flat.txt` with SMILES
- Creates `reaction_intermediates/<RXN_ID>/rxn.smiles` for each reaction
- Creates `from_species_with_cmp`, `to_species_with_cmp`, and
  `species_id_inchikey.txt` in each reaction directory

#### 2. Run RDT on all reactions

```bash
cd /data/AraCore

# This runs RDT on every reaction_intermediates/<RXN_ID>/rxn.smiles
# The original script uses a hardcoded path; invoke RDT directly:
for rxn_folder in reaction_intermediates/*; do
    cd "$rxn_folder"
    smiles=$(cat rxn.smiles)
    java -jar $RDT_JAR -Q SMI -q "$smiles" -g -c -b -j AAM -f TEXT
    # ... post-processing happens inline in run_rdt.sh
    cd /data/AraCore
done
```

The original `run_rdt.sh` script does both the RDT call and the
post-processing (splitting the `.rxn` output into MOL files, computing
InChI/InChIKey via `obabel`, and building the mapping). You can run it
directly after adjusting the JAR path inside the script:

```bash
# Quick fix: replace the hardcoded path in the script
sed -i "s|java -jar ../../../../FluxAndPoolSizeEstimation/ReactionDecoder/target/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar|java -jar \$RDT_JAR|g" /data/AraCore/run_rdt.sh

cd /data/AraCore
bash run_rdt.sh
```

#### 3. Unite mappings

```bash
cd /data/AraCore
bash unite_mappings.sh
```

This produces `all_mapping.txt` and `all_mapping.sorted.txt`, and the
nitrogen-specific files `all_mapping.N.sorted.txt` etc.

### MetaCyc Pipeline

The MetaCyc pipeline follows the same pattern but with additional steps
for comparing against MetaCyc's own atom mappings.

#### 1. Prepare

```bash
cd /data/MetaCyc

# Generate InChIKey tables from SMILES
bash create_inchikey_table_from_smiles.sh

# Prepare reaction SMILES and intermediates
bash prepare_rdt_metacyc.sh
```

#### 2. Run RDT

```bash
# Adjust the JAR path
sed -i "s|java -jar ../../../../FluxAndPoolSizeEstimation/ReactionDecoder/target/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar|java -jar \$RDT_JAR|g" /data/MetaCyc/run_rdt_metacyc.sh
sed -i "s|java -jar ../../../../FluxAndPoolSizeEstimation/ReactionDecoder/target/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar|java -jar \$RDT_JAR|g" /data/MetaCyc/postprocess_rdt_metacyc.sh

cd /data/MetaCyc
bash run_rdt_metacyc.sh
```

> **Warning:** MetaCyc has ~14,875 reactions. Running RDT on all of them
> takes many hours. You can test with a smaller subset by creating a
> `rxns_to_correct.txt` with just a few reaction IDs.

#### 3. Post-process and compare mappings

```bash
cd /data/MetaCyc

# Create InChI equivalent atom tables
bash create_inchi_equivalent_atoms.sh

# Convert existing MetaCyc mappings to the same format
bash convert_existing_mappings.sh

# Sort and compare
bash sort_mappings.sh
bash unite_converted_mappings.sh

# Find discrepancies
bash find_mapping_issues.sh
bash diff_converted_mappings.sh
```

---

## Pipeline overview

```
┌─────────────────────────────────────────────────────────────┐
│                      Input Data                             │
│  species_id_smiles.txt  (species → SMILES)                  │
│  rxn_table.smiles.txt   (reaction → SMILES reaction)        │
│  atom_mappings/         (MetaCyc reference, for comparison) │
│  ArabidopsisCoreModel.xml  (SBML model, for downstream)     │
└──────────────────────┬──────────────────────────────────────┘
                       │
          ┌────────────▼────────────┐
          │   1. PREPARE             │
          │   prepare_rdt.sh         │
          │   - Generate rxn.smiles  │
          │   - Build InChIKey table │
          │   - obabel -:"$smiles"   │
          │     -oinchikey           │
          └────────────┬────────────┘
                       │
          ┌────────────▼────────────┐
          │   2. RUN RDT             │
          │   run_rdt.sh             │
          │   java -jar rdt-2.5.0... │
          │   -Q SMI -q "$smiles"    │
          │   -g -c -b -j AAM -f TEXT│
          │   → ECBLAST_smiles_AAM.rxn│
          └────────────┬────────────┘
                       │
          ┌────────────▼────────────┐
          │   3. POST-PROCESS        │
          │   (inline in run_rdt.sh) │
          │   a) csplit .rxn → MOL_* │
          │   b) obabel → InChI+Key  │
          │   c) Identify species    │
          │   d) Build mapping.txt   │
          │      format:             │
          │      from-spec:Elt#n=... │
          └────────────┬────────────┘
                       │
          ┌────────────▼────────────┐
          │   4. UNITE & ANALYZE     │
          │   unite_mappings.sh      │
          │   → all_mapping.sorted   │
          │   → all_mapping.N.sorted │
          └─────────────────────────┘
```

### RDT CLI flags explained

```
java -jar rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar \
    -Q SMI           Input format: SMILES
    -q "$smiles"     Reaction SMILES string (educts>>products)
    -g               Generate PNG image (fails in headless Docker, see gotcha #13)
    -c               Complex mode (use rings)
    -b               Accept reactions with no bond changes (transporters)
    -j AAM           Job type: Atom-Atom Mapping
    -f TEXT          Output format: TEXT (.rxn file)
```

### Output file: `ECBLAST_smiles_AAM.rxn`

RDT produces a single `.rxn` file containing all mapped molecules in MDL
MOL format, separated by `$MOL` headers. The scripts split this file using
`csplit` and then process each molecule individually with `obabel`.

### Output file: `mapping.txt`

The final per-reaction mapping file uses the format:
```
from_SPECIES:ELEMENT#index=to_SPECIES:ELEMENT#index,...
```
Where:
- `from_` / `to_` prefix is stripped, side is indicated by `=` (from→to) and `,` (separator)
- `SPECIES:ELEMENT#index` identifies a specific atom in a metabolite using
  InChI ordering (element-wise, 1-indexed per element)

---

## Key gotchas & hurdles

### 1. RDT version pinning

**Must use v2.5.0.** The output file naming convention
(`ECBLAST_smiles_AAM.rxn`) and the internal `.rxn` file format changed in
v3.0+ (package namespace moved from `uk.ac.ebi` to `com.bioinceptionlabs`).
The post-processing scripts depend on the v2.5.0 `.rxn` structure.

### 2. Split zip files (MetaCyc)

The MetaCyc `reaction_intermediates` archive is split into two parts
(`.zip.01` and `.zip.02`) because it exceeds GitHub's file size limit.
They must be concatenated before extraction:

```bash
cat MetaCyc/reaction_intermediates.zip.01 MetaCyc/reaction_intermediates.zip.02 \
    > MetaCyc/reaction_intermediates.zip
```

### 3. Hardcoded relative JAR path

All scripts reference the RDT JAR at:
```
../../../../FluxAndPoolSizeEstimation/ReactionDecoder/target/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar
```
This is relative to `reaction_intermediates/<RXN_ID>/`, going 4 levels up.
Inside Docker, the JAR is at `/opt/rdt/`. Use `sed` to patch the scripts
or set `$RDT_JAR` and adjust the `java -jar` calls.

### 4. RDT writes output to the current working directory

RDT always creates `ECBLAST_smiles_AAM.rxn` (and `.png`, `.txt`) in the
**current working directory**. The `run_rdt.sh` scripts `cd` into each
reaction's directory before calling RDT, which ensures the output lands
next to the `rxn.smiles` input. Do **not** run RDT from a different
directory.

### 5. `csplit` pattern sensitivity

The scripts split the `.rxn` file on `$MOL` boundaries:
```bash
csplit -f MOL_ ECBLAST_smiles_AAM.rxn '/$MOL/' {*}
```
The `MOL_00` file (before the first `$MOL`) is always removed as it
contains only the header. If RDT produces no output (e.g., for an
unparseable SMILES), `csplit` will fail and the script will continue with
stale files from a previous run. The scripts attempt `rm MOL_*` before
`csplit`, but if no `.rxn` file exists, the loop body will error.

### 6. OpenBabel InChIKey lookup for species identification

After RDT produces mapped molecules, `obabel` computes InChIKeys which
are then looked up in `species_id_inchikey.txt` (or
`species_id_without_cmp_inchikey.txt`) to recover the original metabolite
ID. This lookup can fail if:
- The SMILES contains generic groups like `[R]` (RDT cannot fully map these)
- The InChIKey computed from the RDT output differs from the one computed
  from the input SMILES (stereochemistry handling, protonation states)

The scripts have a fallback: if the full InChIKey doesn't match, they try
matching only the first 14 characters (the connectivity hash, ignoring
stereochemistry).

### 7. Hydrogen atoms are excluded from mappings

The final mapping explicitly filters out hydrogen:
```bash
grep -v ':H#'
```
This is because RDT often adds/removes protons during mapping, and the
hydrogen mapping is typically not biologically meaningful.

### 8. MetaCyc atom mapping files are mostly empty

Many files in `atom_mappings/` are 0 bytes. This means MetaCyc did not
provide an atom mapping for that reaction. The pipeline skips these.

### 9. Equivalent atoms (symmetry)

The `replace_equiv_atoms.sh` script (a massive sed pipeline) handles
molecules with symmetrically equivalent atoms (e.g., the two oxygens in a
carboxylate group). Without this, the mapping comparison would falsely
report mismatches for equivalent atom positions. This is a
manually-curated mapping that is specific to the MetaCyc dataset.

### 10. Empty SMILES reactions

Some reactions in `rxn_table.smiles.txt` have empty SMILES (e.g.,
`1.10.2.2-RXN	>>`). These represent reactions where no SMILES could be
constructed for the metabolites. RDT will fail on these; the scripts
simply skip them or produce an empty mapping.

### 11. Compartment suffixes in AraCore

AraCore species IDs include compartment tags like `[h]`, `[c]`, `[m]`.
The `prepare_rdt.sh` script strips these for SMILES conversion but
preserves them for the flux analysis. The compartment tag affects the
mapping lookup in the downstream MATLAB scripts.

### 12. `recreate_data.sh` regenerates InChIKey tables

This script re-computes `species_id_without_cmp_inchikey.txt` by calling
`obabel` on every entry in `species_id_smiles.txt`. It takes a while for
MetaCyc (~16k entries) but is relatively fast for AraCore (~236 entries).

### 13. `-g` flag (PNG image generation) fails in headless Docker

The `-g` flag tells RDT to generate a PNG image of the atom mapping. This
requires a graphical environment (AWT/X11) and **will throw an exception
in a headless Docker container**. The `.rxn` file is still written
successfully, so this is non-fatal. However, to avoid the error, simply
remove the `-g` flag when running inside Docker:

```bash
# Instead of: java -jar $RDT_JAR -Q SMI -q "$smiles" -g -c -b -j AAM -f TEXT
# Use:        java -jar $RDT_JAR -Q SMI -q "$smiles"    -c -b -j AAM -f TEXT
```

### 14. RDT creates a `.jnati` directory in the user's home

On first run, RDT extracts the JNI-InChI native library to
`~/.jnati/repo/jniinchi/`. This is ~75 KB and is cached across runs. In
Docker, this ends up in `/root/.jnati/` and is ephemeral (not persisted
between container runs unless you mount the home directory).

---

## Downstream MATLAB analysis

The files in `AraCore/` include several MATLAB scripts for simulated
Metabolic Flux Analysis (simMFA):

| File | Purpose |
|------|---------|
| `createAraCoreNumbers.m` | Main script: loads SBML model, builds atom transition matrix, runs simMFA |
| `simMFA.m` | Simulates ¹⁵N labeling enrichment over time using atom tracking |
| `create_S_N.m` | Builds the nitrogen atom stoichiometric matrix from RDT mappings |
| `createMappingFigure.m` | Generates figures for the paper |
| `create_numbers.mat` | Pre-computed number tables |
| `graph_analysis.m` | Network analysis of atom transition graphs |

These scripts depend on:
- **MATLAB** (or GNU Octave with compatible toolboxes)
- **COBRA Toolbox** (`initCobraToolbox`, `readSBML`, `optimizeCbModel`, `fluxVariability`)
- Optionally **Gurobi** solver

To run these, you would need to install MATLAB/Octave and the COBRA
Toolbox separately. The SBML model file is at
`AraCore/ArabidopsisCoreModel.xml`. The atom mapping input for these
scripts is `AraCore/all_mapping.N.sorted.txt`.

---

## Reference

- **RDT (Reaction Decoder Tool):** Rahman, S.A. et al. (2016) *Bioinformatics* 32(12):i205–i213. DOI: [10.1093/bioinformatics/btw096](https://doi.org/10.1093/bioinformatics/btw096)
- **RDT GitHub:** https://github.com/asad/ReactionDecoder
- **MetaCyc:** Caspi, R. et al. (2018) *Nucleic Acids Research* 46(D1):D633–D639.
- **AraCore model:** Cheung, C.Y.M. et al. (2013) *Plant Physiology* 162(3):1630–1642.
