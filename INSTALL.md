# UniversalRDT — Installation & Running Guide

This pipeline computes **universal atom-to-atom mappings** for metabolic
reactions using the Reaction Decoder Tool (RDT), then post-processes the
output into InChI-ordered atom mapping tables.

> **Important:** This guide covers the **bash-based atom mapping pipeline**
> only. Downstream analysis scripts (`createAraCoreNumbers.m`, `simMFA.m`,
> `create_S_N.m`) are written in MATLAB and require the COBRA Toolbox. See
> [Downstream MATLAB analysis](#downstream-matlab-analysis) for details.

---

## Quick Start

```bash
make install                    # install Java, OpenBabel, unzip (Ubuntu 24.04)
make run-aracore                # extract data + run AraCore pipeline (~236 reactions)
# or for MetaCyc (~14,875 reactions, takes hours):
make run-metacyc
```

### Docker alternative

```bash
make docker                     # build the Docker image
make run-aracore-docker         # run AraCore inside Docker
# or:
make run-metacyc-docker         # run MetaCyc inside Docker
```

---

## Makefile targets

| Target | Description |
|--------|-------------|
| `make install` | Install native dependencies (Java 21, OpenBabel, unzip) via apt |
| `make docker` | Build Docker image with all dependencies |
| `make prepare-aracore` | Extract AraCore reaction intermediates |
| `make prepare-metacyc` | Rejoin split zips + extract MetaCyc data |
| `make run-aracore` | Extract data + run RDT + unite mappings for AraCore |
| `make run-metacyc` | Extract data + run RDT for MetaCyc |
| `make postprocess-metacyc` | Run MetaCyc post-processing (compare with reference mappings) |
| `make compare-metacyc` | Run RDT + full post-processing for MetaCyc |
| `make run-aracore-docker` | Run AraCore pipeline inside Docker |
| `make run-metacyc-docker` | Run MetaCyc pipeline inside Docker |
| `make clean` | Remove all extracted data |

---

## Prerequisites

### Native (Ubuntu 24.04)

- `make install` handles everything: `openjdk-21-jre-headless`, `openbabel`, `unzip`
- **~500 MB** free disk for extracted reaction intermediates (MetaCyc + AraCore)

### Docker

- Docker (or Podman) installed
- `make docker` builds the image (~250 MB)

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
          └────────────┬────────────┘
                       │
          ┌────────────▼────────────┐
          │   2. RUN RDT             │
          │   run_rdt.sh             │
          │   java -jar $RDT_JAR     │
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
java -jar $RDT_JAR \
    -Q SMI           Input format: SMILES
    -q "$smiles"     Reaction SMILES string (educts>>products)
    -g               Generate PNG image of the mapping
    -c               Complex mode (use rings)
    -b               Accept reactions with no bond changes (transporters)
    -j AAM           Job type: Atom-Atom Mapping
    -f TEXT          Output format: TEXT (.rxn file)
```

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

### 1. Unicode filename issues in zip archives

The `reaction_intermediates.zip` file contains filenames with special Unicode characters (e.g., `MOL_.mdl`, `MOL_.inchi`) that can cause "mismatching filename" errors during extraction. This happens when the zip file was created on a system with different locale settings than the extraction system.

**Solution:** Use `unzip` with the `-O UTF-8` flag to explicitly specify UTF-8 encoding:
```bash
unzip -o -O UTF-8 reaction_intermediates.zip -d AraCore/
```

The Makefile has been updated to use this flag automatically.

### 2. RDT version pinning

**Must use v2.5.0.** The output file naming convention
(`ECBLAST_smiles_AAM.rxn`) and the internal `.rxn` file format changed in
v3.0+ (package namespace moved from `uk.ac.ebi` to `com.bioinceptionlabs`).
The post-processing scripts depend on the v2.5.0 `.rxn` structure. The JAR
is committed to this repo.

### 2. Split zip files (MetaCyc)

The MetaCyc `reaction_intermediates` archive is split into two parts
(`.zip.01` and `.zip.02`) because it exceeds GitHub's file size limit.
`make prepare-metacyc` handles rejoining them automatically.

### 3. RDT writes output to the current working directory

RDT always creates `ECBLAST_smiles_AAM.rxn` (and `.png`, `.txt`) in the
**current working directory**. The `run_rdt.sh` scripts `cd` into each
reaction's directory before calling RDT. Do **not** run RDT from a different
directory.

### 4. `csplit` pattern sensitivity

The scripts split the `.rxn` file on `$MOL` boundaries:
```bash
csplit -f MOL_ ECBLAST_smiles_AAM.rxn '/$MOL/' {*}
```
If RDT produces no output (e.g., for an unparseable SMILES), `csplit` will
fail. The scripts attempt `rm MOL_*` before `csplit`.

### 5. OpenBabel InChIKey lookup for species identification

After RDT produces mapped molecules, `obabel` computes InChIKeys which are
looked up in `species_id_inchikey.txt` to recover the original metabolite
ID. This lookup can fail if:
- The SMILES contains generic groups like `[R]`
- The InChIKey differs due to stereochemistry or protonation states

The scripts have a fallback: if the full InChIKey doesn't match, they try
matching only the first 14 characters (connectivity hash, ignoring
stereochemistry).

### 6. Hydrogen atoms are excluded from mappings

The final mapping explicitly filters out hydrogen (`grep -v ':H#'`) because
RDT often adds/removes protons during mapping.

### 7. MetaCyc atom mapping files are mostly empty

Many files in `atom_mappings/` are 0 bytes (no mapping provided by MetaCyc).
The pipeline skips these.

### 8. Equivalent atoms (symmetry)

The `replace_equiv_atoms.sh` script handles molecules with symmetrically
equivalent atoms (e.g., the two oxygens in a carboxylate group). Without
this, the mapping comparison would falsely report mismatches.

### 9. Empty SMILES reactions

Some reactions have empty SMILES (e.g., `1.10.2.2-RXN	>>`). RDT will fail
on these; the scripts skip them or produce an empty mapping.

### 10. Compartment suffixes in AraCore

AraCore species IDs include compartment tags like `[h]`, `[c]`, `[m]`.
The `prepare_rdt.sh` script strips these for SMILES conversion but preserves
them for flux analysis.

### 11. PNG generation in Docker

The `-g` flag generates PNG images of atom mappings. This requires a display.
In Docker, `xvfb-run` provides a virtual display so PNG generation works
without a real X server.

---

## Downstream MATLAB analysis

The files in `AraCore/` include MATLAB scripts for simulated Metabolic Flux
Analysis (simMFA):

| File | Purpose |
|------|---------|
| `createAraCoreNumbers.m` | Main script: loads SBML model, builds atom transition matrix, runs simMFA |
| `simMFA.m` | Simulates ¹⁵N labeling enrichment over time using atom tracking |
| `create_S_N.m` | Builds the nitrogen atom stoichiometric matrix from RDT mappings |

These require **MATLAB** (or GNU Octave) with the **COBRA Toolbox**.
The atom mapping input is `AraCore/all_mapping.N.sorted.txt`.

---

## Reference

- **RDT (Reaction Decoder Tool):** Rahman, S.A. et al. (2016) *Bioinformatics* 32(12):i205–i213. DOI: [10.1093/bioinformatics/btw096](https://doi.org/10.1093/bioinformatics/btw096)
- **RDT GitHub:** https://github.com/asad/ReactionDecoder
- **MetaCyc:** Caspi, R. et al. (2018) *Nucleic Acids Research* 46(D1):D633–D639.
- **AraCore model:** Cheung, C.Y.M. et al. (2013) *Plant Physiology* 162(3):1630–1652.
