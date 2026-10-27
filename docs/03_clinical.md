# 03 — Clinical

Builds the survival table from the cBioPortal clinical file: one row per patient, with the chosen endpoint stored in the
columns `OS_STATUS` (1 = event) and `OS_MONTHS`, whatever the endpoint.

| Reads | Writes (to `<cancer>_output/`) |
|---|---|
| `<cancer>_data/<cBioPortal study>/data_clinical_patient.txt` | `clinical_data.csv`, `clinical_endpoint.txt` (which endpoint and source were used) |

Settings: `CANCER`; `ENDPOINTS = {"THCA": "PFS", "PRAD": "PFS"}` (others use overall survival); `CANCER_SETTINGS` maps
each cancer to its cBioPortal study (COAD uses the combined colorectal study, so rectal patients are present in the
clinical table and drop out when matched to the COAD expression data).

Note: BRCA, KIRC, LUAD and UCEC ran an earlier version of this notebook that did not write `clinical_endpoint.txt`;
they use overall survival.

## Functions

None: the notebook is inline code.

miner functions called: none.
