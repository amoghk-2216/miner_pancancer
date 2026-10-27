# MINER pan-cancer

The gbmMINER causal and mechanistic regulatory network analysis of glioma
([Turkarslan et al., medRxiv 2024](https://doi.org/10.1101/2024.04.05.24305380)) applied to ten TCGA cancers, with
expression from recount3. The glioblastoma replication that this work extends is in a separate repository,
[miner_replica](https://github.com/amoghk-2216/miner_replica).

For each cancer the pipeline infers regulons (a transcription factor and its co-expressed target genes), causal flows
(somatic mutation → transcription factor → regulon), transcriptional states and programs. It then decides which regulons
are disease-relevant, validates them in an independent cohort where one exists, builds survival risk models, finds master
regulators per state and ranks candidate regulators (notebook 07, an extension that is not part of gbmMINER).

## Repository layout

```
README.md
requirements.txt
driver_tiers.py            published driver genes per cancer (IntOGen via Open Targets); imported by notebooks 02, 06, 07
notebooks/<cancer>/        01_<cancer>_expr … 07_<cancer>_candidate_regulators: the notebooks as each cancer ran, with outputs
template/                  01_expression … 07_candidate_regulators: one generic set for any of the ten cancers (or a new one)
results/<cancer>/          the smaller result tables of each cancer (see results/README.md)
docs/                      function reference: one page per notebook, plus driver_tiers.py
examples/                  example_analysis.ipynb: worked analysis of results/ that needs no raw data
scripts/list_projects.py   lists TCGA projects in recount3 by sample count (how the ten cancers were chosen)
```

`<cancer>` is one of `brca`, `kirc`, `luad`, `lusc`, `ucec`, `thca`, `prad`, `hnsc`, `lgg`, `coad`. Input data,
caches and the larger outputs are not in the repository (see "Rerunning a cancer").

## Where to start

- **The results:** the table below, then `results/<cancer>/` (described in `results/README.md`).
- **A worked example:** `examples/example_analysis.ipynb` compares the ten cancers and follows one candidate regulator
  down to its regulons, hazard ratios, upstream driver mutations and states, using only `results/`. Start Jupyter from
  `examples/`.
- **How a step works:** the notebook in `notebooks/<cancer>/`, and its page in `docs/` for the functions it defines.
- **Running the pipeline yourself:** start from `template/` (see "Rerunning a cancer").

## Cancers and validation cohorts

| Cancer | TCGA project | Patients | Endpoint (events) | Validation cohort |
|---|---|---|---|---|
| Breast | BRCA | 1,093 | OS (151) | METABRIC |
| Kidney clear cell | KIRC | 531 | OS (177) | E-MTAB-1980 |
| Lung adenocarcinoma | LUAD | 516 | OS (184) | GSE72094 |
| Lung squamous | LUSC | 501 | OS (215) | GSE157011 |
| Endometrial | UCEC | 541 | OS (91) | CPTAC 2020 (co-expression only, no survival data) |
| Thyroid | THCA | 505 | PFS (51) | none |
| Prostate | PRAD | 497 | PFS (93) | GSE116918 (biochemical recurrence) |
| Head and neck | HNSC | 502 | OS (219) | GSE65858 |
| Lower-grade glioma | LGG | 514 | OS (125) | Gravendeel (GSE16011), grade II–III only |
| Colon | COAD | 458 | OS (98 ¹) | GSE39582 |

Patients are tumours with expression after notebook 01; events are counted in `clinical_data.csv` (notebook 03).
¹ Among the 435 COAD patients with both expression and clinical data (the cBioPortal study also contains rectal tumours).

## Results summary

| Cancer | Regulons | Programs | States | Disease-relevant regulons | Programs prognostic: TCGA / cohort / both | Risk C-index, TCGA 5-fold CV (gene / regulon / program) | Risk C-index, cohort (gene / regulon / program) | TFs with replicated prognostic activity |
|---|---|---|---|---|---|---|---|---|
| LGG | 8,990 | 278 | 23 | 3,115 | 161 / 110 / 85 | 0.79 / 0.80 / 0.75 | 0.72 / 0.63 / 0.72 | 232 |
| KIRC | 8,995 | 314 | 25 | 1,661 | 184 / 79 / 58 | 0.67 / 0.67 / 0.66 | 0.77 / 0.80 / 0.75 | 130 |
| LUAD | 7,447 | 291 | 24 | 1,068 | 97 / 81 / 43 | 0.59 / 0.62 / 0.59 | 0.59 / 0.67 / 0.66 | 92 |
| PRAD | 8,422 | 238 | 20 | 242 | 84 / 24 / 11 | 0.64 / 0.67 / 0.68 | 0.56 / 0.51 / 0.65 | 13 |
| BRCA | 4,596 | 157 | 25 | 208 | 19 / 54 / 10 | 0.58 / 0.55 / 0.59 | 0.52 / 0.55 / 0.54 | not computed ² |
| HNSC | 7,609 | 284 | 21 | 264 | 53 / 31 / 6 | 0.60 / 0.60 / 0.57 | 0.51 / 0.48 / 0.50 | 18 |
| COAD | 7,739 | 216 | 22 | 146 | 24 / 24 / 5 | 0.55 / 0.54 / 0.51 | 0.52 / 0.52 / 0.53 | 5 |
| LUSC | 6,644 | 257 | 20 | 26 | 26 / 23 / 4 | 0.51 / 0.52 / 0.50 | 0.55 / 0.51 / 0.51 | 2 |
| UCEC | 7,557 | 218 | 27 | 1,208 ³ | – | 0.62 / 0.51 / 0.55 | – | 0 |
| THCA | 9,292 | 256 | 20 | 542 ⁴ | – | 0.53 / 0.57 / 0.58 | – | 0 |

² BRCA's notebook 07 predates the replicated-activity column. ³ Criteria 1–3 and 5: CPTAC has no survival data.
⁴ Criteria 1, 2 and 5: no independent cohort.

The method works and replicates in **KIRC, LUAD and PRAD**. **LGG** gives the strongest, best-replicating signal, most
likely driven by IDH status and grade. It works **partly in BRCA, HNSC and COAD**: single regulons replicate, but the
multivariable risk models transfer weakly or not at all (in COAD, tumour stage alone gives C = 0.70 while the expression
models stay near 0.5). It is **weak in LUSC and THCA**, and **UCEC** cannot be tested for survival outside TCGA. Only five
TFs reach RRA FDR ≤ 0.05 in notebook 07 (FOXA2 and POU3F3 in UCEC, SOX7 and TBX15 in THCA, SOX21 in HNSC), so the set of
TFs whose activity replicates is the more useful output of that notebook.

## The pipeline

| # | Notebook | What it does | Main outputs |
|---|---|---|---|
| 01 | `01_<cancer>_expr` | recount3 counts → TPM (exon-union length) → log2 → z-score; protein-coding genes, primary tumours, one sample per patient | `<cancer>_zscored.csv`, `<cancer>_log2tpm.csv` |
| 02 | `02_<cancer>_regulons_clusters_states` | MC3 mutation matrix, co-expression clusters, regulons, causal flows, states, programs | `regulons.json`, `regulonDf.csv`, `filteredCausalResults.csv`, `transcriptional_states.json`, `transcriptional_programs.json` |
| 03 | `03_<cancer>_clinical` | survival table from cBioPortal (OS, or PFS for THCA and PRAD) | `clinical_data.csv` |
| 04 | `04_<cancer>_disease_relevance` | regulon Cox, GO hallmark similarity, co-expression and survival replication in the cohort | `final_disease_relevant_regulons.csv` |
| 05 | `05_<cancer>_risk_model` | gene, regulon and program ridge models on GuanRank scores; repeated splits, 5-fold CV, cohort | `risk_models/` |
| 06 | `06_<cancer>_regulators` | driver-filtered causal flows, TF–TF network, master regulators per state | `master_regulators_per_state.csv` |
| 07 | `07_<cancer>_candidate_regulators` | TFs ranked by Robust Rank Aggregation over five evidence lines | `candidate_regulators.csv` |

A regulon is disease-relevant under the paper's rule when it is (1) co-expressed in TCGA, (2) associated with survival in
TCGA, (3) co-expressed in the cohort, (4) associated with survival in the cohort in the same direction, and (5) its TF
correlates with it (|rho| ≥ 0.3). Cancers without a cohort with survival data use the criteria that can be applied.

## Rerunning a cancer

Use the notebooks in `template/`. They are the newest version (the one COAD, HNSC, LGG, PRAD, THCA and LUSC ran) with
every cancer's validation cohort added, and their outputs cleared. Set `CANCER` in the first code cell of each notebook;
everything else follows from it. The notebooks in `notebooks/<cancer>/` are the record of what each cancer actually ran
and are kept as they were.

```
pip install -r requirements.txt
```

`isb-miner3` needs pandas below version 3.

The notebooks read and write a working folder outside the repository, set as `base_dir` in each notebook's first code
cell (and the Gene Ontology folder as `GO_DIR` in 04). They are set to `base_dir = r"causal_inf_new"` and
`GO_DIR = r"GBM_new"`, relative to the folder Jupyter is started from; change them to your own folders (absolute paths
work too). The working folder must look like this:

```
<base_dir>/
  driver_tiers.py                  copy from this repository (the notebooks import it from base_dir)
  mc3.v0.2.8.PUBLIC.maf            TCGA MC3 somatic mutations, shared by all cancers
  SI_geneset.gmt                   optional, ESTIMATE stromal/immune signatures used by 07
  <cancer>/
    <cancer>_data/                 written by 01 and 02, plus the files below placed by hand
      <cBioPortal study>/data_clinical_patient.txt
      <cancer>_<cohort>/           validation cohort download
    <cancer>_output/               everything 02–07 write
```

| File placed by hand | Source |
|---|---|
| `mc3.v0.2.8.PUBLIC.maf` | TCGA MC3, GDC PanCanAtlas publications page |
| `data_clinical_patient.txt` | cBioPortal study named in notebook 03 (`CANCER_SETTINGS`) |
| `go-basic.obo`, `goa_human.gaf` (in `GO_DIR`) | Gene Ontology |
| validation cohort | METABRIC and CPTAC 2020 from cBioPortal; E-MTAB-1980 from ArrayExpress; the others as GEO SOFT family files |
| `SI_geneset.gmt` (optional) | ESTIMATE R package, `inst/extdata` |

Run 01 → 07 in order and restart the kernel after editing `driver_tiers.py`. Notebook 02 has a miRNA section that runs
only when `<cancer>_mirna_log2rpm.csv` exists; the miRNA notebook that builds it is not included and was not run for these
results, so the section is skipped. Approximate times on the machine used: 01–03 under an hour together, 04 about 1.5–2
hours on a first run (GO enrichment is then cached), 05 15–30 minutes, 06 and 07 a few minutes each. The first call to
`driver_tiers.py` for a cancer needs internet access to reach Open Targets.

How the template relates to the runs:

- Notebook 04's cohort loader was merged from the per-cancer copies. On the saved cohort files it reproduces the z-scored
  cohort matrix and clinical table of eight cohorts exactly (BRCA, LUAD, LUSC, UCEC, PRAD, HNSC, LGG, COAD). KIRC's
  E-MTAB-1980 branch is copied verbatim from the KIRC notebook but could not be tested here, because it needs miner's
  identifier table.
- BRCA, KIRC, LUAD and UCEC ran older versions of notebooks 03, 04, 05 and 07, and ran 02 with
  `RESTRICT_TO_DRIVERS = False` where the template has `True`, so rerunning them from the template will not give identical
  numbers. The other six cancers ran the template's code apart from the cancer name and their own cohort branch.

## Differences from gbmMINER

1. **TPM gene length** is the union of a gene's exons, not its genomic span.
2. **Per-patient mean centring** after z-scoring (notebook 02); the paper describes no such step.
3. **Driver-gene filter.** Causal flows are limited to IntOGen drivers of each cancer from Open Targets (minimum IntOGen
   score 0.6, a judgement rather than a published threshold), because in heavily mutated cancers most frequently mutated
   genes are passengers. BRCA, KIRC, LUAD and UCEC ran 02 with `RESTRICT_TO_DRIVERS = False` and applied the filter in 06;
   LUSC, THCA, PRAD, HNSC, LGG and COAD applied it in 02. The driver-filtered flows used by 06 and 07 should be the same
   either way, since each mutation is tested separately; this has not been checked by rerunning one cancer both ways.
4. **TF expression filter** (06): a TF must reach 1 TPM in at least 20% of tumours.
5. **Endpoints:** progression-free survival for THCA and PRAD, which have few deaths; PRAD's cohort records biochemical
   recurrence.
6. **Validation cohorts** are one per cancer, on different platforms; UCEC's has no survival data and THCA has none; LGG
   uses Gravendeel grade II–III tumours only.
7. **Risk models** are judged over repeated 70/10/20 splits and 5-fold cross-validation, not one split.
8. **Notebook 07** is an addition. Its activity evidence was not adjusted for stromal and immune content, because
   `SI_geneset.gmt` was absent in every run.

## Caveats

- Causal flows are associations in patient data, confounded by tumour subtype and purity; several recurrent master
  regulators (SPI1, MAFB, SOX7, SOX17) mark immune or endothelial cells rather than tumour cells.
- LGG's results probably reflect the known IDH / 1p19q / grade classes; whether the regulons add to them is untested.
- HNSC's risk models may fail to transfer because GSE65858 holds more HPV-positive tumours than TCGA.
- BRCA has few deaths (151 of 1,095), which limits every survival step.
- `base_dir` and `GO_DIR` are relative folder names (`causal_inf_new`, `GBM_new`); set them before running. Saved outputs
  in `notebooks/<cancer>/` print these same shortened paths.

## References

- Turkarslan et al. 2024, medRxiv, doi:10.1101/2024.04.05.24305380 (gbmMINER)
- Wall et al. 2021, MINER method (multiple myeloma)
- Wilks et al. 2021, recount3
- Ellrott et al. 2018, TCGA MC3 mutation calls
- Liu et al. 2018, TCGA Clinical Data Resource and endpoint recommendations
- Martínez-Jiménez et al. 2020, IntOGen
- Sondka et al. 2018, Cancer Gene Census
- Ochoa et al., Open Targets Platform
- Guan et al., GuanRank
- Kolde et al. 2012, Robust Rank Aggregation
- Yoshihara et al. 2013, ESTIMATE; Barbie et al. 2009, ssGSEA
