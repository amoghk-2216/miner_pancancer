# 01 — Expression (recount3)

Builds the tumour expression matrix: recount3 gene counts → TPM (exon-union gene length) → log2(TPM + 1) → protein-coding
genes, one row per gene symbol, one primary tumour per patient → zero-variance genes dropped → z-score per gene.

| Reads | Writes (to `<cancer>_data/`) |
|---|---|
| recount3 (downloaded): sample metadata, QC table, TCGA table, gene counts, GTF annotation | `<cancer>_tpm_unfiltered.csv`, `<cancer>_log2tpm.csv`, `<cancer>_zscored.csv` |

Settings: `CANCER`; `AVG_READ_LEN_COL = "star.average_mapped_length"` (converts recount3 coverage into read counts).

## Functions

### `collapse_group(sub_df, corr_threshold=0.8)`

Cell: `# Collapse duplicate gene names into one row each`

Collapses the duplicate rows of one gene symbol into one expression profile. Rows correlated at or above the threshold
are linked, connected groups form clusters, each cluster is averaged, and the cluster whose average varies most across
samples is kept. Duplicates can be unrelated loci sharing a symbol or the same gene annotated twice, which is why rows are
not simply averaged.

| Argument | Type | Meaning |
|---|---|---|
| `sub_df` | DataFrame (n_duplicates × n_samples) | log2(TPM + 1) rows (Ensembl IDs) sharing one symbol |
| `corr_threshold` | float, default 0.8 | minimum correlation for two rows to join one cluster |

Returns `(value, n_clusters, largest_cluster_size)`: the collapsed profile (Series, n_samples), the number of clusters the
rows split into, and the size of the largest cluster. A single row is returned unchanged with `1, 1`.

Note: the docstring says Spearman, but the code calls `DataFrame.corr()` with its default, which is Pearson.

miner functions called: none.
