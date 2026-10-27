# 04 — Disease-relevant regulons

Decides which regulons are disease-relevant under the paper's rule: (1) co-expressed in TCGA, (2) associated with survival
in TCGA (Cox, p ≤ 0.05), (3) co-expressed in the validation cohort, (4) associated with survival in the cohort in the same
direction, and (5) the regulon's TF correlates with it (|rho| ≥ 0.3 within the regulon's coherent patients). It also
computes the hallmark criterion (GO enrichment + Lin similarity to the ten cancer hallmarks) used by the earlier, looser
rule.

| Reads | Writes (to `<cancer>_output/`) |
|---|---|
| `eigengenes.csv`, `regulons.json`, `regulonDf.csv`, `overExpressedMembers.csv`, `underExpressedMembers.csv` (02); `clinical_data.csv` (03); `<cancer>_zscored.csv`; `go-basic.obo`, `goa_human.gaf` (`GO_DIR`); the validation cohort files | `regulon_cox_results.csv`, `disease_relevant_regulons_survival_only.csv`, `disease_relevant_regulons_survival_or_hallmark.csv`, `go_bio_symbol.pkl`, `lin_cache.pkl` (+ `_checkpoint` versions), `<cohort>_zscored.pkl`, `<cohort>_clinical.csv`, `<cohort>_random_regulon_variance.csv`, `<cohort>_regulon_coherence.csv`, `<cohort>_coexp_validated_regulons.csv`, `tcga_regulon_coherence.csv`, `<cohort>_regulon_activity.pkl`, `<cohort>_regulon_coherence_members.csv`, `<cohort>_regulon_cox_results.csv`, `regulator_eigengene_correlation.csv`, `final_disease_relevant_regulons.csv`, `final_disease_relevant_regulons_loose.csv`, `validation_status.txt` |

Settings: `CANCER`, `GO_DIR`, `VALIDATION_COHORTS` (cancer → cohort name and folder; a cancer that is not listed runs
without external validation), `N_JOBS = 16`; GO enrichment `alpha = 0.05`, `methods = ['fdr_bh']`; hallmark Lin
similarity ≥ 0.8; cohort null distribution `N_ITER = 250` gene sets per regulon size; co-expression permutations
`N_PERM = 500`; `MIN_REG_RHO = 0.3`, `MIN_MEMBERS = 5`.

The GO and Lin caches take about an hour to build on a first run and are reused afterwards; delete them after any rerun of
02, since regulon IDs change.

## Functions

### `study_terms(study_genes, assoc)`

Cell: `# Helper fn: tracks GO terms in a particular geneset`

Returns the set of GO term IDs annotated to any gene of a gene set. `study_genes`: list of symbols; `assoc`: dict
{symbol: set of GO IDs}.

### `bh_adjust_with_true_m(pvals, m_true)`

Cell: `# Helper fn: tracks GO terms in a particular geneset`

Benjamini–Hochberg q-values where the number of tests is given explicitly (every GO term in the background) rather than
the number of p-values passed in (only the terms the regulon's genes carry). Returns an array of q-values in the input
order, monotone and capped at 1.

### `regulon_hallmark_sim_fast(reg_id, enriched_go, threshold=0.8)`

Cell: `# Takes one regulon id, its list of sig GO terms and a similarity threshold (here it is 80%)`

Finds which of the ten hallmarks a regulon qualifies for: a hallmark passes if any of the regulon's enriched GO terms has
Lin similarity ≥ `threshold` with any of the hallmark's GO terms. Returns `(reg_id, list of hallmark names)`.

Uses: `lin_cache` (precomputed Lin similarities), `hallmarks_updated` ({hallmark: GO IDs}), `godag`.

### `parse_geo_soft(soft_file)`

Cell: `# Load the validation cohort: expression with one row per gene symbol, plus clinical data if available`

Reads a GEO SOFT family file (plain or `.gz`): the platform table (probe → gene symbol, column found by name) and every
sample's expression table and characteristics. Probes mapping to more than one symbol are dropped; the rest are averaged
per symbol.

Returns `(expression, clinical)`: genes × samples DataFrame indexed `Hugo_Symbol`, and a samples × characteristics
DataFrame indexed `PATIENT_ID` (values as strings, parsed by the cohort-specific code that follows).

In `template/`, defined once in the loader cell and used by every GEO branch; it records a sample's GEO sub-series as
`sub_series` when it differs from `VAL_NAME`, which GSE157011 (LUSC, two sub-series) needs. In `notebooks/<cancer>/` it
is defined inside the cohort branch of COAD, PRAD, HNSC, LGG and LUSC (only LUSC's records the sub-series); BRCA, KIRC,
LUAD and UCEC load their cohorts with inline code, and THCA has no cohort.

The template's loader cell handles every cohort: GSE157011, GSE116918 / GSE39582 / GSE65858, Gravendeel (LGG), GSE72094,
E-MTAB-1980 and cBioPortal studies (METABRIC, CPTAC). Run on the saved cohort files, it reproduces the per-cancer
results exactly for eight cohorts; the E-MTAB-1980 branch is a verbatim copy of KIRC's and was not retested.

### `_null_for_size(A, genes, n, n_iter, seed_base)`

Cell: `# Null distribution:`

Worker process (joblib). Draws `n_iter` random gene sets of size `n` from the cohort expression and returns their variance
and PC1 variance explained (miner's `random_regulon_variance`), as a DataFrame with columns `N`, `Variance`,
`Variance_explained`. The seed is `seed_base + n`, so each size is reproducible.

| Argument | Type | Meaning |
|---|---|---|
| `A` | ndarray (genes × samples) | z-scored cohort expression |
| `genes` | list | row labels of `A` |
| `n`, `n_iter`, `seed_base` | int | set size, number of sets, seed base |

The code is the same in every cancer; only its indentation differs between versions.

### `_coherence_one(X, reg_id, g_idx, m_idx, n_perm, seed)`

Cell: `import os, json, time, gzip` (criterion 1, TCGA co-expression)

Worker process. Co-expression test for one regulon: the share of variance on PC1 of the regulon's genes across its
coherent patients, against two permutation nulls (random gene sets of the same size in the same patients, and the same
genes in random patient sets). Random sets are drawn in memory-bounded chunks.

| Argument | Type | Meaning |
|---|---|---|
| `X` | ndarray (genes × patients) | TCGA expression |
| `reg_id` | int | regulon ID |
| `g_idx`, `m_idx` | ndarray | row positions of the regulon's genes, column positions of its coherent patients |
| `n_perm`, `seed` | int | permutations per null, random seed |

Returns a dict: `regulon`, `n_genes`, `n_coherent_patients`, `pc1_var_explained`, `p_random_genes`, `p_random_patients`
(p = (1 + exceedances) / (n_perm + 1); NaN when fewer than 3 genes or 3 patients).

### `pc1_share(A)`

Defined inside `_coherence_one` (and again at notebook level). Share of variance on the first principal component for a
batch of gene sets, `A` of shape (n_sets, n_genes, n_patients): the top eigenvalue of the smaller Gram matrix divided by
the trace. Returns one value per set.

### `eigengene(M)`

Cell: `import os, json, time, gzip` (criterion 3, cohort co-expression)

First principal component of a regulon's genes across patients, oriented to rise with the genes' mean expression and
scaled to unit standard deviation. `M`: ndarray (genes × patients), z-scored. Returns an array (n_patients,).

## miner functions called

`random_regulon_variance` (inside `_null_for_size`), `generateRegulonActivity` (cohort regulon activity); the notebook
imports `miner.miner` with `*`.
