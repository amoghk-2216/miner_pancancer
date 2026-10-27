# 06 — Regulators: TF–TF network and master regulators

Keeps the causal flows that start from a published driver gene, builds a TF–TF network (A → B when A has a binding site in
B's promoter and A's expression predicts B's in a LASSO model), and names up to three master regulators per transcriptional
state: TFs with high out/in degree in that state's network.

| Reads | Writes (to `<cancer>_output/`) |
|---|---|
| `filteredCausalResults.csv`, `eigengenes.csv`, `transcriptional_states.json`, `overExpressedMembers.csv`, `underExpressedMembers.csv` (02); `<cancer>_zscored.csv`, `<cancer>_log2tpm.csv`; `driver_tiers.py`; miner's identifier mappings and TF-binding databases (tfbsdb, tfbsdb2, chea_seq) | `driverCausalResults.csv`, `tf_expression_filter.csv`, `tf_tf_network.csv`, `tf_network_degrees.csv`, `tf_tf_network_per_state.csv`, `master_regulator_state_summary.csv`, `master_regulators_per_state.csv` |

Settings: `MIN_TPM = 1.0`, `MIN_FRAC = 0.20` (a TF must reach 1 TPM in at least 20% of tumours); master regulators
`N_MASTER = 3`, `STATE_P = 0.05`, `STATE_FRAC = 0.5`, `STATE_TF_RULE = "enriched"`, `TF_ENRICH_P = 0.05`,
`LASSO_SAMPLES = "all"`, `MIN_STATE_N = 20`, `CV_FOLDS = 5`.

## Functions

### `pruned_lasso(target, predictors, X, y, cv_folds=5, seed=42)`

Cell: `# Master regulators per state, following the paper (Methods, "Identification of master regulators from TF-TF network")`

LASSO of one target TF's expression on its candidate regulators (TFs with a binding site in its promoter), with
cross-validation over shuffled folds and the one-standard-error rule: the largest penalty whose mean CV error is within one
standard error of the best. Predictors with no variation are left at 0. Runs in worker processes and imports inside.

| Argument | Type | Meaning |
|---|---|---|
| `target` | str | target TF symbol |
| `predictors` | list of str | candidate regulator symbols |
| `X` | ndarray (patients × predictors) | predictor expression |
| `y` | ndarray (patients,) | target expression |
| `cv_folds`, `seed` | int | folds and shuffling seed |

Returns `(target, predictors, coefficients, n_nonzero_at_best)`: the coefficients after pruning, aligned with
`predictors`, and the number of predictors kept at the CV-best penalty before pruning.

## miner functions called

None by name: the notebook imports `miner.miner` with `*` and reads miner's data files (`identifier_mappings.txt`,
`network_dictionaries/*.pkl`) directly.
