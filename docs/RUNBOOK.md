# RUNBOOK — mirna_glioma

Step-by-step instructions to reproduce the analysis from scratch.

---

## Prerequisites

1. **Conda environment**:

   ```bash
   conda activate omics-R
   # Or create project-specific env:
   # conda env create -f env/environment.yml && conda activate mirna_glioma_env
   ```

2. **Additional R packages** (install once):

   ```r
   BiocManager::install(c("ComplexHeatmap", "circlize", "rrvgo", "GO.db", "org.Hs.eg.db"))
   install.packages("rbioapi")
   ```

3. **Input files** (place in `data/intermediate/`):
   - `Gliomas_all_counts_merged.csv` — raw count matrix with annotation columns
   - `Metadatos_gliomas_verificados.csv` — clinical metadata

4. **Configuration**: Review `config/de_specs.csv` for comparison definitions.

---

## Execution Order

All commands run from the project root (`~/bioinfo/projects/mirna_glioma/`).

### Phase 1: QC and Preprocessing

```bash
Rscript scripts/01_inspeccion_counts.R data/intermediate/Gliomas_all_counts_merged.csv
Rscript scripts/02_logcpm_mds.R data/intermediate/Gliomas_all_counts_merged.csv
```

### Phase 1b: Age-Matched Group Design (Optional)

For comparisons requiring age-matched case-control groups (e.g., EGFR_VIII):

```bash
Rscript scripts/02b_matching_EGFR_VIII.R \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --ratio 2 \
  --caliper 15
```

This generates:
- Matching diagnostics in `results/tables/matching/` and `results/figures/matching/`
- List of excluded samples for `config/de_specs.csv`

Then update `config/de_specs.csv` to add the matched comparison with excluded samples.

### Phase 2: Differential Expression

```bash
Rscript scripts/03_edgeR_multiDE.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --outdir results/DE
```

### Phase 3: Visualization

```bash
Rscript scripts/04_heatmaps_sigRNAs.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE

Rscript scripts/05_deg_qc_plots.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE
```

### Phase 4: Survival Analysis

```bash
Rscript scripts/08_KM_DE_candidates.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --fdr_cut 0.1
```

### Phase 5: Functional Enrichment (DE-based)

```bash
Rscript scripts/09_miEAA_GSEA_all_comparisons.R \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --out_root results/tables/MiEAA_GSEA

Rscript scripts/10_miEAA_GSEA_bubble_plots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_bubble \
  --run_tag A_conservative

Rscript scripts/11_miEAA_GSEA_barplots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_QC \
  --run_tag A_conservative \
  --cutoff 0.25

Rscript scripts/12_miEAA_GSEA_categories_pvals.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_categories_pvals \
  --run_tag A_conservative
```

### Phase 6: Functional Enrichment (Survival-based)

```bash
Rscript scripts/13_survGSEA.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --out_root results \
  --run_tag SurvivalRank_CoxZ_miRPathDB

Rscript scripts/14_survGSEA_plots.R \
  --in_root results/tables/SurvivalRank_GSEA \
  --out_root results/figures/SurvivalRank_GSEA \
  --run_tag SurvivalRank_CoxZ_miRPathDB \
  --n_mirpathdb 12 \
  --preset double_col

Rscript scripts/15_survGSEA_enrichment_plot.R \
  --in_root results/tables/SurvivalRank_GSEA \
  --out_root results/figures/SurvivalRank_GSEA \
  --run_tag SurvivalRank_CoxZ_miRPathDB \
  --category_regex "GO Biological process" \
  --enrichment depleted \
  --preset single_col

Rscript scripts/16_survGSEA_reduce_redundancy.R \
  --sim_threshold 0.9 \
  --sim_method Rel \
  --run_tag SurvivalRank_CoxZ_miRPathDB
```

### Phase 7: Pathway-miRNA Reports (Optional)

Generate summary reports with miRNAs annotated by rank and statistic value:

```bash
Rscript scripts/17_generate_pathway_mirna_reports.R \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --gsea_root results/tables/MiEAA_GSEA \
  --surv_root results/tables/SurvivalRank_GSEA \
  --run_tag A_conservative \
  --top_n 12
```

This generates:
- `results/tables/reports/GSEA_DE_top12_annotated_*.csv` — DE-based pathways with miRNA ranks and logFC
- `results/tables/reports/GSEA_survGSEA_top12_annotated_*.csv` — Survival-based pathways with miRNA ranks and z-scores

---

## Reproduction Checklist

- [ ] Conda environment activated (`omics-R` or `mirna_glioma_env`)
- [ ] Input data files present in `data/intermediate/`
- [ ] `config/de_specs.csv` matches intended comparisons
- [ ] R and package versions match those in `docs/METHODS_scripts.md`
- [ ] Scripts executed in order (phases 1-6)
- [ ] Logs in `logs/` checked for errors after each script
- [ ] Output figures in `results/figures/` verified visually

---

## Notes

- Scripts 06 and 07 (exploratory survival: Spearman, Cox univariate) are excluded from the final pipeline.
- All scripts generate timestamped logs in `logs/`.
- Scripts 09 and 13 call the miEAA API (`rbioapi`); internet access required.
- Random seed default: 42 (configurable via `--seed`).
