# Function reference

One page per notebook, listing the functions each notebook defines, what they take and return, and which notebook-level
variables they rely on. The pages apply both to `template/` and to the per-cancer copies in `notebooks/<cancer>/`: the
functions are the same in all of them except where a page says otherwise.

| Page | Notebook (template / per cancer) | Functions defined |
|---|---|---|
| [01_expression.md](01_expression.md) | `01_expression` / `01_<cancer>_expr` | `collapse_group` |
| [02_regulons_clusters_states.md](02_regulons_clusters_states.md) | `02_regulons_clusters_states` / `02_<cancer>_…` | `resolve`, `family_members`, `mirna_causal_flows` |
| [03_clinical.md](03_clinical.md) | `03_clinical` / `03_<cancer>_clinical` | none (inline code only) |
| [04_disease_relevance.md](04_disease_relevance.md) | `04_disease_relevance` / `04_<cancer>_…` | `study_terms`, `bh_adjust_with_true_m`, `regulon_hallmark_sim_fast`, `parse_geo_soft`, `_null_for_size`, `_coherence_one`, `pc1_share`, `eigengene` |
| [05_risk_model.md](05_risk_model.md) | `05_risk_model` / `05_<cancer>_risk_model` | `program_activity`, `split_70_10_20`, `rank_features`, `cindex`, `fit_one_split`, `km_panel`, `tune_fit_predict` |
| [06_regulators.md](06_regulators.md) | `06_regulators` / `06_<cancer>_regulators` | `pruned_lasso` |
| [07_candidate_regulators.md](07_candidate_regulators.md) | `07_candidate_regulators` / `07_<cancer>_…` | `bh_fdr`, `read_id_list`, `joined`, `ssgsea_scores`, `regulon_scores`, `prepare_cohort`, `tf_activity`, `score_test`, `e4_stat`, `rank_aggregation`, `evidence`, `classify`, `colour` |
| [driver_tiers.md](driver_tiers.md) | `driver_tiers.py` (module) | `get_driver_tiers`, `resolve_disease`, `fetch_associations`, `_post`, `_apply_min_score` |

## Conventions

- **Cell** names the code cell a function lives in by its first comment line, which is how the notebooks are navigated.
- **Uses** lists notebook-level variables the function reads without taking them as arguments. A function with a Uses line
  cannot be lifted into another notebook without those variables; the others are self-contained.
- Functions whose names start with `_` run in worker processes (joblib) and import what they need inside the function.
- miner functions called by each notebook are listed at the end of its page; their documentation is in the miner package
  ([baliga-lab/miner3](https://github.com/baliga-lab/miner3)).
