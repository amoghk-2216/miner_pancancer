
"""
Step 1 of the top-10-cancer pilot: list TCGA projects available in recount3
with their sample counts, so we pick the pilot cohort (largest non-GBM
project) from real data instead of assuming.
 
Run this locally (same machine/environment as your existing notebooks):
 
    pip install recount3
    python step1_list_tcga_projects.py
 
It will print the top TCGA projects by sample count and save the full table
to tcga_projects_recount3.csv in the current directory. Send me that CSV
(or just the printed top 10) and I'll confirm the pilot cancer from it.
"""
 
import recount3
import pandas as pd
 
pd.set_option("display.max_rows", 100)
pd.set_option("display.width", 200)
 
df = recount3.available_projects(organism="human", data_sources="tcga")
 
# One row per TCGA project code (e.g. BRCA, LUAD, GBM, ...) with n_samples
df_sorted = df.sort_values("n_samples", ascending=False).reset_index(drop=True)
 
df_sorted.to_csv("tcga_projects_recount3.csv", index=False)
 
print(f"Total TCGA projects found: {len(df_sorted)}\n")
print("Top 15 by sample count:")
print(df_sorted[["project", "n_samples"]].head(15).to_string(index=False))
 
print("\nWhere does GBM rank (for reference, already done):")
print(df_sorted[df_sorted["project"].str.contains("GBM", case=False, na=False)][["project", "n_samples"]])
 
print("\nSaved full table to tcga_projects_recount3.csv")
 

