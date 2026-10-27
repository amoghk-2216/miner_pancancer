# 07 — Candidate regulators (extension, not in gbmMINER)

Ranks transcription factors by Robust Rank Aggregation over five lines of evidence, with a permutation p-value that
reassigns regulons to TFs at random:

- **E1** the TF's regulons are enriched for disease-relevant regulons (hypergeometric);
- **E2** its survival-associated regulons point the same way (two-sided binomial test, risk vs protective);
- **E3** driver-mutation causal flows into its disease-relevant regulons;
- **E4** its activity, the signed mean of its regulons, predicts survival (Cox score test, stratified across TCGA and the
  cohort, optionally adjusted for stromal and immune content);
- **E5** it is a master regulator in 06.

`act_replicated` marks TFs whose activity predicts survival in TCGA and in the cohort separately, in the same direction.

| Reads | Writes (to `<cancer>_output/`) |
|---|---|
| `regulonDf.csv`, `filteredCausalResults.csv` (02); `regulon_cox_results.csv`, `final_disease_relevant_regulons.csv`, `<cohort>_coexp_validated_regulons.csv`, `<cohort>_zscored.pkl`, `<cohort>_clinical.csv` (04); `clinical_data.csv` (03); `master_regulators_per_state.csv`, `tf_network_degrees.csv` (06); `<cancer>_log2tpm.csv`, `<cancer>_zscored.csv`; `open_targets_associations.csv`; `SI_geneset.gmt` (optional) | `stromal_immune_scores.csv` (only with `SI_geneset.gmt`), `candidate_regulators.csv`, `candidate_regulators_top20.png`, `candidate_shortlist.csv` |

Settings: `MIN_TPM = 1.0`, `MIN_FRAC = 0.20` (expression requirement); `EXT_ENDPOINT = "DSS"` (cohort endpoint where
disease-specific survival exists); `N_PERM = 1000`, `PERM_SEED = 0`.

`SI_geneset.gmt` was absent in every run, so E4 is not adjusted for stromal or immune content in the current results.

## Functions

### `bh_fdr(pvals)`

Cell: `# Load the pipeline outputs`

Benjamini–Hochberg q-values for a Series; NaN p-values are ignored and returned as NaN. Returns a Series with the same index.

### `read_id_list(name)`

Cell: `# Load the pipeline outputs`

Reads a header-less regulon list from the output folder. Returns a set of int. Uses: `output_dir`.

### `joined(s)`

Cell: `# Evidence E3: is the TF downstream of a driver mutation of this cancer?`

Distinct values of a Series as one sorted, comma-separated string.

### `ssgsea_scores(expr, gene_set, alpha=0.25)`

Cell: `# Stromal and immune content of each tumour (ESTIMATE approach, Yoshihara et al. 2013)`

Single-sample gene set enrichment score (Barbie et al. 2009), the statistic ESTIMATE uses. `expr`: DataFrame
(genes × patients), any monotonic scale; `gene_set`: list of symbols; `alpha`: rank weight exponent. Returns a Series of
scores per patient (higher = signature genes rank higher in that tumour).

### `regulon_scores(z)`

Cell: `# Evidence E4: the TF's ACTIVITY predicts survival, tested in TCGA and the independent cohort together`

Mean z-score of each regulon's genes per patient, computed with one sparse matrix product. `z`: DataFrame
(symbols × patients). Returns an array (patients × regulons) in the order of `regulon_ids`; 0 where a regulon has fewer
than 2 measured genes. Uses: `regulon_ids`, `reg_genes`.

### `prepare_cohort(name, S, time, event, covariates)`

Same cell. Precomputes what the Cox score test needs for one cohort: the time ordering (longest survival first), the
risk-set positions of the events (Breslow ties), and the projection onto the covariates (intercept, plus standardised
covariates if given). Returns a dict with keys `name`, `S`, `order`, `events`, `last`, `n_risk`, `Z`, `Z_pinv`,
`n_events`.

### `tf_activity(cohort, W, n_regs)`

Same cell. TF activity per patient: the signed mean of its regulon scores, `(W.T @ S.T).T / n_regs`, where `W` is a sparse
regulons × TFs matrix of ±1: the sign of the correlation between the TF's expression and the regulon score
(`sign_TR = np.sign(tf_x @ Sz)`, zeros set to +1), so repressed regulons count with a minus sign. Returns an array (patients × TFs).

### `score_test(A_by_cohort, cohorts_)`

Same cell. Cox partial-likelihood score test of survival on each column of an activity matrix, stratified by cohort: the
score U and information I are summed over cohorts after each column is standardised and adjusted for the covariates.
Returns z = U / √I per TF (> 0: higher activity, higher hazard).

### `e4_stat(owner)`

Same cell. Stratified score-test z for every TF under a given regulon → TF assignment, overall and per cohort. `owner`:
array (n_regulons,) of indices into the TF table. Returns `(z, [z per cohort])`. Called with the real assignment and with
permuted ones.

Uses: `act_col_of_tf`, `sign_TR`, `regulon_ids`, `act_tfs`, `cohorts`, `tf_activity`, `score_test`.

### `rank_aggregation(norm_ranks)`

Cell: `# Combine the evidence and rank: Robust Rank Aggregation (Kolde et al. 2012) score, with a PERMUTATION p-value`

Robust Rank Aggregation rho per item. `norm_ranks`: array (items × lists) of ranks divided by the number of items, 1.0
where an item has no evidence in a list. For each item, the k-th smallest normalised rank is compared with its
Beta(k, n − k + 1) distribution; rho is the smallest of these probabilities times the number of lists, capped at 1.

### `evidence(owner)`

Same cell. The five evidence lines for the ranked TFs under one regulon → TF assignment, returned as normalised ranks
(array n_ranked × 5, 1.0 = no evidence) and a dict of the raw values `E1`–`E5`. Used for the observed ranking and for each
permutation.

Uses: `is_final`, `is_risk`, `is_protect`, `flows_into`, `n_regs_all`, `ranked_pos`, `n_tf_all`, `regulon_ids`,
`act_pos_of_ranked`, `masters`, `e4_stat`.

### `classify(score)`

Cell: `# Annotation only: prior evidence for this cancer in Open Targets (does not change the ranking)`

Labels a gene by its Open Targets association score with the cancer: "Known (strong)" (≥ 0.5), "Known (moderate)"
(≥ 0.1), "Weak prior evidence" (> 0), "No prior association" (0).

### `colour(r)`

Cell: `# Top 20 candidate regulators`

Bar colour for the top-20 plot from a row of the TF table: grey when `act_fdr` > 0.05, orange when higher activity goes
with worse survival (`act_z` > 0), blue otherwise.

## miner functions called

None: the notebook works from the output files of 02–06.
