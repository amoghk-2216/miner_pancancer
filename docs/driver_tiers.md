# driver_tiers.py — published driver genes per cancer

Module imported by notebooks 02, 06 and 07. It returns the driver genes of a cancer so that causal flows are limited to
genes under selection in that tumour type, rather than every gene mutated in at least 5% of patients (in heavily mutated
cancers almost all of those are passengers that track mutation burden).

A gene counts as a driver only if IntOGen calls it one for that cancer (Martínez-Jiménez et al. 2020), as reported by the
Open Targets Platform. The Cancer Gene Census channel (Sondka et al. 2018) then splits the calls into two tiers:

| Tier | Meaning |
|---|---|
| `known` | IntOGen driver call for this cancer and in the Cancer Gene Census |
| `candidate` | IntOGen driver call for this cancer, not in the Cancer Gene Census |
| `provisional` | Open Targets could not be reached; a built-in fallback list was used (not citable, not cached) |

The table is cached as `driver_genes.csv` in the cancer's data folder, together with `open_targets_associations.csv`
(used by notebook 07 for annotation). Delete the cache, or pass `refresh=True`, to fetch again.

Settings: `OT_QUERIES` (Open Targets search term per cancer), `OT_DISEASE_IDS` (fixed disease IDs where the search is
ambiguous; BRCA → `MONDO_0004989`), `MIN_INTOGEN_SCORE = 0.6` (a judgement, not a published threshold: established drivers
score 0.68 or more in BRCA, KIRC and LUAD, while known passengers such as CSMD3 and FAT4 score lower), `PROVISIONAL_DRIVERS`
(fallback lists for LUAD, KIRC, BRCA, UCEC, THCA).

## Functions

### `get_driver_tiers(cancer, data_dir, refresh=False)` (public)

Returns the driver table for a cancer, from the cache if present, otherwise from Open Targets, filtered to IntOGen
score ≥ `MIN_INTOGEN_SCORE`.

| Argument | Type | Meaning |
|---|---|---|
| `cancer` | str | TCGA project code, e.g. `"LUAD"` |
| `data_dir` | str | the cancer's data folder (cache location) |
| `refresh` | bool | ignore the cache and fetch again |

Returns a DataFrame, one row per driver: `Gene`, `Tier`, `cgc_score`, `intogen_score`, `ot_score`, `Source`,
`Disease_ID`. If Open Targets cannot be reached and no cache exists, returns the provisional list, or an empty table for a
cancer without one (then no driver filter can be applied).

Check the disease name it prints on the first call for a cancer.

### `resolve_disease(cancer)`

Finds the Open Targets disease for a cancer: the fixed ID if one is set, otherwise an exact name match, otherwise the
other search hits, each checked to still exist. Returns `(disease_id, disease_name)`; raises `LookupError` if none exists.

### `fetch_associations(disease_id)`

Downloads every target associated with the disease, 500 per page. Returns a DataFrame: `Gene`, `ot_score` (overall),
`cgc_score` (Cancer Gene Census), `intogen_score`.

### `_post(query, variables, retries=4)`

Sends one GraphQL query to the Open Targets API, retrying with a growing pause on HTTP 429, 502, 503 and 504. Returns the
response's `data`; raises on query errors or when every attempt fails.

### `_apply_min_score(table)`

Keeps the rows with `intogen_score ≥ MIN_INTOGEN_SCORE` and prints how many of each tier remain.
