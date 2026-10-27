# 02 — Regulons, causal flows, states and programs

The core MINER run: mutation matrix from the TCGA MC3 calls, co-expression clusters, regulons, causal flows
(mutation → TF → regulon), transcriptional states and programs, and an optional miRNA layer.

| Reads | Writes |
|---|---|
| `<cancer>_zscored.csv`, `<cancer>_log2tpm.csv` (01); `mc3.v0.2.8.PUBLIC.maf` (shared); `driver_tiers.py`; miner's identifier mappings and TF-target database; `<cancer>_mirna_log2rpm.csv` (optional; the miRNA notebook that builds it is not in the repository) | `<cancer>_data/`: `<cancer>_mutation_matrix.csv`, `<cancer>_mutation_matrix_all.csv`, `driver_genes.csv`. `<cancer>_output/`: `coexpressionDictionary.json`, `coexpressionModules.json`, `coregulationModules.json`, `regulons.json`, `regulonDf.csv`, `mechanisticOutput.json`, `eigengenes.csv`, `completeCausalResults.csv`, `filteredCausalResults.csv`, `overExpressedMembers.csv`, `underExpressedMembers.csv`, `transcriptional_states.json`, `transcriptional_programs.json`, `programs_vs_states.csv/.pdf`, `regulons_activity_heatmap.csv`, `regulon_activity_heatmap.pdf`, `state_coherence_check.csv`, `mutation_filter_log.csv`; miRNA files only when the miRNA matrix exists (not the case for these results) |

Settings: `CANCER`; `RESTRICT_TO_DRIVERS` (True: only driver genes are tested for causal flows; False: every gene mutated
in ≥ 5% of patients, as in the paper); mutation frequency ≥ 5% of patients; miner parameters
`min_number_genes = 6`, `min_correlation = 0.2`, `min_number_regulon_genes = 5`, `num_cores = 5`; causal filter at
p ≤ 0.05 for each test; miRNA layer `MIRNA_MAX_R = -0.3`, `MIRNA_MIN_MEMBERS = 5`.

Each patient's expression is mean-centred after z-scoring, in the cell `# Diagnostic: per-sample mean distribution`.
The paper describes no such step.

## Functions

### `resolve(g)`

Cell: `# Gene symbol -> Ensembl ID conversion with a synonym fallback`

Converts one gene symbol to its Ensembl ID: the primary gene name first, then unambiguous synonyms. Returns the Ensembl ID,
or the symbol unchanged if neither maps (such genes are dropped later).

Uses: `symbol_to_ensembl`, `synonym_to_ensembl` (dicts built in the same cell from miner's identifier mappings).

### `family_members(family_key, available)`

Cell: `# miRNA expression and the miRNA families of miner's target database`

Finds the miRNAs of the expression matrix that belong to one seed family of miner's miRNA target database.

| Argument | Type | Meaning |
|---|---|---|
| `family_key` | str | database key: seed then member names, e.g. `"ACATTCC_MIR1_MIR206"` |
| `available` | set of str | miRNA names in the matrix without copy number, e.g. `"hsa-mir-125b"` |

Returns a sorted list of matching names. The arm suffix (3P / 5P) is dropped, and old names are resolved where possible
(MIR124A → hsa-mir-124, MIR103 → hsa-mir-103a).

### `mirna_causal_flows(mutations, regulator_expr, eigen, regulon_table, significance=0.05)`

Cell: `# Causal flows mutation -> miRNA -> regulon, with the same tests and output columns as miner's causal step`

Tests every (mutation, miRNA) pair the way miner tests (mutation, TF) pairs. A flow is reported when the miRNA differs
between mutated and wild-type patients (Welch t-test), the regulon eigengene does too, and the signs agree:
(mutation effect on miRNA) × (miRNA–regulon correlation) × (mutation effect on regulon) > 0. Mutations with fewer than
3 mutated or 3 wild-type patients are skipped.

| Argument | Type | Meaning |
|---|---|---|
| `mutations` | DataFrame (genes × patients) | 1 = mutated, 0 = wild type |
| `regulator_expr` | DataFrame (miRNAs × patients) | miRNA expression |
| `eigen` | DataFrame (regulons × patients) | regulon eigengenes |
| `regulon_table` | DataFrame, index = regulon ID | columns `Regulator`, `Spearman_R`, `Spearman_p` |
| `significance` | float, default 0.05 | p-value threshold of both t-tests |

Returns a DataFrame with one row per flow and the same columns as miner's `completeCausalResults.csv`; regulons are named
`MR-<id>`.

Self-contained (uses `scipy.stats` as `sstats`, imported in the notebook's first cell).

## miner functions called

`cluster`, `reviseInitialClusters`, `getCoexpressionModules`, `getCoregulationModules`, `getRegulons`,
`regulonDictionary`, `principal_df`, `getEigengenes`, `mechanisticInference`, `background_df`,
`biclusterMembershipDictionary`, `membershipToIncidence`, `matrix_to_dictionary`, `mosaic`, `transcriptionalPrograms`,
`reduceModules`, `programsVsStates`, `read_json`, `write_json`, `read_pkl`; also `miner.causal_inference` and
`miner.subtypes.inferSubtypes` are imported.
