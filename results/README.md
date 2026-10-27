# Result tables

Copies of the smaller output files of each cancer's run, one folder per cancer (`brca`, `kirc`, `luad`, `lusc`, `ucec`,
`thca`, `prad`, `hnsc`, `lgg`, `coad`). Regulon IDs are those of that cancer's run and are not comparable between
cancers. In `regulonDf.csv`, regulators and genes are Ensembl IDs; the files written by notebooks 06 and 07 use gene
symbols. `<cohort>` is the validation cohort's prefix (for example `gse116918`).

| File | Written by | Contents |
|---|---|---|
| `regulonDf.csv` | 02 | one row per (regulon, regulator, gene): `Regulon_ID`, `Regulator`, `Gene` |
| `transcriptional_states.json` | 02 | {state: patient barcodes} |
| `transcriptional_programs.json` | 02 | {program: regulon IDs} |
| `programs_vs_states.csv` | 02 | program activity in each state |
| `mutation_filter_log.csv` | 02 | genes mutated in ≥ 5% of patients, their driver tier and whether they were kept (not in BRCA, KIRC, LUAD, UCEC, which ran an older notebook 02) |
| `clinical_endpoint.txt` | 03 | endpoint and source study (not in BRCA, KIRC, LUAD, UCEC, which use overall survival) |
| `regulon_cox_results.csv` | 04 | Cox hazard ratio (per unit of the regulon's eigengene, as written by miner) and p-value per regulon in TCGA |
| `tcga_regulon_coherence.csv` | 04 | co-expression of each regulon in TCGA: PC1 share within its coherent patients and two permutation p-values |
| `<cohort>_coexp_validated_regulons.csv` | 04 | regulon IDs co-expressed in the cohort (earlier all-patient test; no header) |
| `<cohort>_regulon_cox_results.csv` | 04 | Cox hazard ratio (per SD of the regulon's eigengene in the cohort) and p-value per regulon in the cohort |
| `final_disease_relevant_regulons.csv` | 04 | regulon IDs passing the paper's rule (no header) |
| `final_disease_relevant_regulons_loose.csv` | 04 | regulon IDs passing the earlier, looser rule (no header) |
| `validation_status.txt` | 04 | which criteria could be applied |
| `program_cox_results.csv` | 05 | Cox hazard ratio and p-value per program in TCGA |
| `program_cox_tcga_vs_validation.csv` | 05 | the same in both cohorts, with significance and replication flags (not in UCEC, THCA) |
| `risk_models/risk_model_summary.csv` | 05 | chosen k and ridge penalty, TCGA test C-index, and cohort C-index and log-rank p per model (not in UCEC) |
| `risk_models/tcga_cv_summary.csv` | 05 | 5-fold cross-validation: out-of-fold C-index, log-rank p, per-fold mean and SD |
| `risk_models/tcga_repeated_splits.csv` | 05 | test C-index of each of the 10 repeated splits |
| `risk_models/km_*.pdf` | 05 | Kaplan–Meier curves of the final models (TCGA test split, TCGA cross-validation, cohort) |
| `driverCausalResults.csv` | 06 | causal flows from a driver mutation: `Mutation`, `Regulator` (symbol), `Regulon`, `Tier` |
| `tf_network_degrees.csv` | 06 | out-degree, in-degree and their ratio for each TF in the TF–TF network |
| `master_regulators_per_state.csv` | 06 | up to three master regulators per state, with degrees and expression contrast |
| `master_regulator_state_summary.csv` | 06 | per state: patients, regulons, TFs, LASSO models and network edges |
| `candidate_regulators.csv` | 07 | one row per TF: evidence E1–E5, activity statistics (`act_z`, `act_fdr`, `act_replicated`), RRA rho, p and FDR, Open Targets annotation |
| `candidate_regulators_top20.png` | 07 | the top 20 TFs by RRA p-value |
| `candidate_shortlist.csv` | 07 | TFs meeting the shortlist rule (LUAD and UCEC only; the other runs predate this file) |

Not included because of size: eigengenes, the over/under-activity matrices, the full causal tables
(`completeCausalResults.csv`, `filteredCausalResults.csv`), z-scored validation cohorts, regulon activity matrices and the
GO and Lin-similarity caches. The notebooks write them when the pipeline is rerun.
