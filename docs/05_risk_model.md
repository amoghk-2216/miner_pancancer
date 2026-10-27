# 05 — Risk models

Fits gene, regulon and program risk models that predict each patient's GuanRank survival score with ridge regression,
evaluates them over repeated 70/10/20 splits and 5-fold cross-validation in TCGA, applies the final TCGA models unchanged
to the validation cohort, and tests each program's association with survival in both cohorts.

| Reads | Writes (to `<cancer>_output/`) |
|---|---|
| `clinical_data.csv` (03); `overExpressedMembers.csv`, `underExpressedMembers.csv`, `regulons.json`, `transcriptional_programs.json` (02); `<cancer>_zscored.csv`; `<cohort>_zscored.pkl`, `<cohort>_clinical.csv` (04) | `guan_scores.csv`; `risk_models/`: `tcga_repeated_splits.csv`, `final_models.pkl`, `km_tcga_test.pdf`, `tcga_cv_folds.csv`, `tcga_cv_oof_predictions.csv`, `km_tcga_cv.pdf`, `tcga_cv_summary.csv`, `<cohort>_pred_<model>.csv`, `km_<cohort>.pdf`, `risk_model_summary.csv`; `<cohort>_regulon_activity.pkl`; `program_cox_results.csv`, `<cohort>_program_cox_results.csv`, `program_cox_tcga_vs_validation.csv`, `<cohort>_program_cox_dss_results.csv` (BRCA only) |

Settings: `N_REPEATS = 10` splits; `FINAL_SEED = 0` (the split whose models are applied to the cohort);
`ALPHAS = [0.1, 1, 5, 10, 50, 100, 500, 1000, 5000]`; `K_GRID`: gene 10–2,000 features, regulon 10–2,000 or all,
program 5–100 or all; `PROGRAM_THRESHOLD = 0.5`; `RISK_CUTOFF = 0.5` (predicted GuanScore above it = high risk);
`N_FOLDS = 5`, `CV_SEED = 0`.

## Functions

### `program_activity(reg_act, programs, threshold=PROGRAM_THRESHOLD)`

Cell: `# Builds feature matrix that program risk model is trained on`

Turns discrete regulon activity (1 / 0 / −1) into discrete program activity with miner's `reduceModules`: a program is +1
in a patient when at least `threshold` of its regulons are over-active, −1 when at least that share are under-active.
The same function builds TCGA and cohort features.

| Argument | Type | Meaning |
|---|---|---|
| `reg_act` | DataFrame (regulons × patients) | regulon activity, index = regulon ID (int) |
| `programs` | list of lists | regulon IDs of each program |
| `threshold` | float, default 0.5 | share of regulons that must agree |

Returns a DataFrame (programs × patients); programs with no regulon present are filled with 0.

### `split_70_10_20(seed)`

Cell: `# Train/val/test split`

70% train, 10% validation, 20% test, stratified by event so each part has the same share of events. Returns three lists of
patient IDs. Uses: `patients`, `evt_all`.

### `rank_features(X, y)`

Cell: `# Train/val/test split`

Spearman correlation of every feature column with the GuanScore, computed on ranks in one matrix product; features that
do not vary get 0. Returns feature column positions sorted by |rho|, largest first. `X`: ndarray (n_train × n_features);
`y`: ndarray (n_train,).

### `cindex(pred, ids, time_col='OS_MONTHS', event_col='OS_STATUS', table=None)`

Cell: `# Train/val/test split`

C-index of predicted GuanScores for the given patients (−pred, since a high score means short survival). `table` defaults
to the TCGA survival table `srv`. Returns a float (0.5 = random, 1 = perfect).

### `fit_one_split(model_name, seed)`

Cell: `# Train/val/test split`

One repeated split: ranks features on train, picks k and the ridge penalty on validation, refits on train only and scores
the test set. `model_name`: "gene", "regulon" or "program".

Returns a dict: `model`, `seed`, `k`, `alpha`, `val_cindex` (optimistic, not reported), `test_cindex` (reported),
`features`, `ridge` (fitted model), `split` (train, validation, test IDs).

Uses: `X_tcga` (feature matrices by model), `y_all` (GuanScores), `K_GRID`, `ALPHAS`, `split_70_10_20`, `rank_features`,
`cindex`.

### `km_panel(ax, pred, table, time_col, event_col, title)`

Cell: `# Final models (FINAL_SEED split) and KM curves on the TCGA test set, as in paper Figure 4`

Splits patients at `RISK_CUTOFF` into high and low risk, draws both Kaplan–Meier curves on `ax` and returns
`(C-index, log-rank p)`; p is NaN when everyone falls in one group. Uses: `RISK_CUTOFF`.

### `tune_fit_predict(model_name, train_ids, test_ids, seed)`

Cell: `# Cross-validation: out-of-fold prediction for every TCGA patient`

One cross-validation fold: holds out 12.5% of the training patients (stratified by event) to pick k and alpha, fits on the
rest and predicts the held-out fold, which is never used for fitting or tuning. Returns `(predictions, k, alpha)`.

Uses: `X_tcga`, `y_all`, `evt_all`, `K_GRID`, `ALPHAS`, `rank_features`, `cindex`.

## miner functions called

`guanRank` (GuanScores), `generateRegulonActivity` (cohort regulon activity), `reduceModules` (inside
`program_activity`).
