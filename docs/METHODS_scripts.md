# Computational Methods

This document describes the computational methods implemented in the `scripts/` directory of the `mirna_glioma` project. The focus is on reproducibility: we detail inputs, filtering criteria, statistical analyses, outputs, and assumptions. The pipeline comprises quality control, differential expression analysis, survival modeling, and functional enrichment.

---

## A) Study Design and Data

### Data Sources

The analysis uses small RNA sequencing (RNA-seq) count data from glioma tumor samples. Input files include:

- **Count matrix**: `data/intermediate/Gliomas_all_counts_merged.csv` containing raw read counts with annotation columns (`id`, `type`, `entrez_id`, `HGNC_symbol`).
- **Sample metadata**: `data/intermediate/Metadatos_gliomas_verificados.csv` containing clinical and demographic variables.

### Sample Selection

Post-mortem control samples were excluded from the analysis. Samples with column names beginning with `A` (matching regex `^A`) are systematically removed during data loading to ensure analysis focuses on tumor specimens.

---

## B) Preprocessing and Quality Control

### Initial Count Inspection (Script 01)

**Script**: `01_inspeccion_counts.R`

Quality assessment of the raw count matrix includes:

1. **Data structure validation**: Detection of annotation columns and construction of a numeric count matrix.
2. **Sample exclusion**: Removal of post-mortem control samples (prefix `A`).
3. **Data integrity checks**: Verification for missing values (NA), negative counts, and non-integer entries.
4. **Library size calculation**: Using edgeR's `DGEList` with TMM (trimmed mean of M-values) normalization to compute effective library sizes.
5. **Zero fraction analysis**: Calculation of the proportion of zero counts per sample and total counts per feature.

**Outputs**: Inspection figures (library size histograms, zero fraction distributions, log1p boxplots, barplots, and MDS plots) and an RDS file containing key inspection objects.

### Normalization and Multivariate QC (Script 02)

**Script**: `02_logcpm_mds.R`

Normalization and quality visualization include:

1. **Sample filtering**: Exclusion of samples with library size < 100,000 reads.
2. **Feature filtering**: Retention of features with:
   - Total count > 0 across samples
   - Non-zero variance
   - Sum of counts >= 10
3. **Normalization**: Calculation of log-transformed counts per million (logCPM) before and after TMM normalization (edgeR) using a prior count of 2.
4. **Quality visualization**: Generation of RLE (relative log expression) plots, density distributions, and library size barplots.
5. **Dimensionality reduction**: MDS (multidimensional scaling) analysis using limma with coordinate export.

**Outputs**: Normalized expression matrix (`logCPM_TMM_*.csv`) and library size tables.

---

## C) Differential Expression Analysis

### Multi-comparison Differential Expression (Script 03)

**Script**: `03_edgeR_multiDE.R`

Differential expression analysis uses the edgeR quasi-likelihood (QL) framework:

1. **Sample alignment**: Matching between count matrix and metadata, with export of filtered datasets.
2. **Comparison specification**: Defined via `--vars` parameter or through `config/de_specs.csv`.
3. **Supported comparison modes**:
   - `as_is`: Variables used directly without transformation
   - `binary_cut`: Binarization by threshold
   - `collapse_levels`: Grouping of factor levels
   - `range_vs`: Comparison of extreme groups (e.g., quartiles)
   - `survival_cut`: Binarization based on survival outcomes
   - `continuous`: Continuous covariate analysis
4. **Statistical framework**: TMM normalization, `filterByExpr` for low-expression filtering, quasi-likelihood modeling with `glmQLFit` and `glmQLFTest`.
5. **Multiple testing**: For multinomial factors, omnibus tests and pairwise comparisons against the reference level; for binary comparisons, direct case vs. reference contrasts.

**Outputs per comparison**: Complete results tables, FDR < 0.05 and FDR < 0.10 filtered tables, MD plots, and volcano plots.

### Biotype-restricted Differential Expression (Script 03_2)

**Script**: `03_2_edgeR_multiDE_miRNA_protCoding.R`

A variant analysis restricted to biologically relevant feature types:

1. **Feature filtering**: Selection of features annotated as `mirna` or `protein_coding` in the `type` column.
2. **Analysis pipeline**: Identical quasi-likelihood framework as script 03.

This restriction focuses the analysis on regulatory small RNAs and their potential protein-coding targets.

---

## D) Multiple Testing Correction

False discovery rate (FDR) correction is applied using the Benjamini-Hochberg (BH) procedure throughout the pipeline:

- Differential expression analyses
- Survival association tests
- Pathway enrichment analyses

Standard significance thresholds of FDR < 0.05 (stringent) and FDR < 0.10 (exploratory) are applied depending on the analytical context.

---

## E) Functional Enrichment Analysis

### miEAA GSEA for Differential Expression (Script 09)

**Script**: `09_miEAA_GSEA_all_comparisons.R`

Gene Set Enrichment Analysis (GSEA) for miRNAs using miEAA (miRNA Enrichment Analysis and Annotation):

1. **miRNA extraction**: Selection of features annotated as miRNA from differential expression results.
2. **Ranking metric**: By default, `signed_sqrtF` (sign of logFC multiplied by square root of F-statistic). Alternative metrics: `signed_logp`, `signed_F`, `signed_logFC`.
3. **GSEA execution**: Via `rbioapi::rba_mieaa_enrich` using miRPathDB expert categories:
   - GO Biological Process (mature miRNAs)
   - GO Molecular Function (mature miRNAs)
   - KEGG pathways (mature miRNAs)
   - Reactome pathways (mature miRNAs)
4. **Category validation**: Verification against `rba_mieaa_cats()`.
5. **Significance filtering**: No pre-filtering applied (`sig_level = 1`); post-hoc filtering recommended.

**Outputs**: Ranked miRNA lists, complete GSEA tables (unfiltered), top 50 by Q-value, and Q < 0.05 filtered results per comparison.

### GSEA Visualization: Bubble Plots (Script 10)

**Script**: `10_miEAA_GSEA_bubble_plots.R`

Publication-ready bubble plots for GSEA results:

1. **Data source**: Most recent `MiEAA_GSEA_all_*.tsv` per comparison.
2. **P-value detection**: Priority order: P-adjusted > Q-value > P-value.
3. **Term selection**: Top N terms per database (default: 12).
4. **Visualization**: Separate bubble plots for enriched (green) and depleted (orange) pathways, plus combined figures using patchwork.
5. **Export formats**: PDF (vector, 300 dpi), SVG (editable vector), PNG (600 dpi raster).

**Styling specifications**: Colorblind-safe Okabe-Ito palette, explicit font sizing (8.5pt base), double-column journal preset (180 x 120 mm).

### GSEA Quality Control: Barplots (Script 11)

**Script**: `11_miEAA_GSEA_barplots.R`

Stacked barplots for QC assessment of enrichment results:

1. **Significance detection**: Automatic detection of P-adjusted, Q-value, or P-value columns.
2. **Direction detection**: Classification of enrichment direction (enriched vs. depleted).
3. **Visualization**: Stacked barplots by database and direction for configurable significance cutoffs.

**Outputs**: PDF, SVG, and PNG exports with logs copied to output directories.

### Category-level P-value Distributions (Script 12)

**Script**: `12_miEAA_GSEA_categories_pvals.R`

Distribution analysis of significance across pathway categories:

1. **Metric visualization**: -log10(Q-value) and -log10(P-adjusted) distributions.
2. **Plot type**: Boxplots with jittered points, ordered by median significance.
3. **Category ordering**: Descending by median -log10(metric).

**Outputs**: PDF, SVG, and PNG exports per metric.

---

## F) Survival Analysis

### Exploratory Correlation Analysis (Script 06)

**Script**: `06_survival_exploratory_spearman.R`

Spearman rank correlation between expression and survival time:

1. **Feature selection**: Union of features with FDR < threshold across all differential expression comparisons.
2. **Correlation metrics**:
   - **All samples**: Spearman correlation including both events and censored observations (`rho_all`, `p_all`, `n_all`).
   - **Events only**: Restricted to patients with recorded death (`MUERTE = 1`) to assess within-event trends (`rho_event`, `p_event`, `n_event`).
3. **Visualization**: Scatter plots with filled points for events, open circles for censored observations, linear regression trend line, and LOESS smoothing curve (both for illustration only).

**Note**: This analysis is exploratory and does not account for censoring in the statistical test.

### Cox Proportional Hazards Modeling (Script 07)

**Script**: `07_survival_cox_univariate.R`

Univariate Cox regression for survival association:

1. **Feature selection**: Union of differentially expressed features (FDR < threshold).
2. **Model specification**: `Surv(time, event) ~ expression` using `survival::coxph` with time variable (`MESES_SEGUIMIENTO_`) and event indicator (`MUERTE`).
3. **Expression scaling**: Optional standardization for hazard ratio interpretation per standard deviation.
4. **Multiple testing correction**: FDR adjustment (Benjamini-Hochberg).

**Outputs**: Hazard ratios with 95% confidence intervals, p-values, FDR-adjusted p-values, and membership tables.

### Kaplan-Meier Analysis (Script 08)

**Script**: `08_KM_DE_candidates.R`

Kaplan-Meier survival curves for differentially expressed candidates:

1. **Candidate selection**: Features with FDR < threshold from differential expression analyses.
2. **Group stratification**: Dichotomization by median expression or extreme quantiles.
3. **Statistical test**: Log-rank test for survival differences between high and low expression groups.
4. **Multiple testing correction**: FDR adjustment across all tested features.

**Outputs**: Kaplan-Meier curve PDFs, individual feature PNGs, and summary tables with p-values and FDR.

### Survival-ranked GSEA (Script 13)

**Script**: `13_survGSEA.R`

GSEA based on survival association rather than differential expression:

1. **Data preparation**: Construction of logCPM (TMM-normalized) expression matrix restricted to miRNAs.
2. **Survival modeling**: Univariate Cox proportional hazards regression per miRNA with expression as the sole predictor.
3. **Ranking metric**: Signed Wald z-statistic from Cox models, where positive values indicate higher expression associated with worse survival.
4. **GSEA execution**: `rbioapi::rba_mieaa_enrich` with `test_type = "GSEA"` using the survival-ranked miRNA list.
5. **Pathway databases**: miRPathDB expert categories (GO BP, GO MF, KEGG, Reactome).

**Key parameters**:
- Expression scaling: Optional standardization (default: enabled) for consistent effect size interpretation.
- Minimum events: Analyses require >= 5 events for stable Cox coefficient estimation.
- Significance filtering: No pre-filtering (`sig_level = 1`); post-hoc filtering recommended.

**Outputs**:
- `Cox_univariate_miRNA_all_*.tsv`: Complete Cox regression results with HR, CI, p-values, and FDR.
- `Ranked_miRNAs_by_CoxZ_*.txt`: Ordered miRNA list by survival association.
- `miEAA_GSEA_all_miRPathDB_expert_*.tsv`: Complete GSEA results.
- `miEAA_GSEA_top50_miRPathDB_expert_*.tsv`: Top 50 pathways by Q-value.

### Survival GSEA Visualization (Script 14)

**Script**: `14_survGSEA_plots.R`

Publication-ready visualization of survival-ranked GSEA results:

1. **Input**: `miEAA_GSEA_all_miRPathDB_expert_*.tsv` from script 13.
2. **Plot types**:
   - **Individual plots**: Separate bubble plots for enriched (green) and depleted (orange) pathways.
   - **Combined figure**: Side-by-side comparison using patchwork.
   - **Category distributions**: Boxplots of -log10(Q-value) and -log10(P-adjusted) by pathway database.
3. **Styling**:
   - Preset: `double_col` (180 x 120 mm base, journal double-column format).
   - Colors: Okabe-Ito colorblind-safe palette (enriched: #009E73, depleted: #D55E00).
   - Typography: 8.5pt base, Helvetica font family.
   - Combined figure dimensions: 450 x 336 mm with increased Y-axis spacing for term readability.

**Export formats**: PDF (cairo_pdf, 300 dpi), SVG (svglite, editable), PNG (600 dpi).

**Output naming convention**:
- `SurvivalRank_Bubble_miRPathDB_enriched_*.{pdf,svg,png}`
- `SurvivalRank_Bubble_miRPathDB_depleted_*.{pdf,svg,png}`
- `SurvivalRank_Bubble_miRPathDB_combined_*.{pdf,svg,png}`
- `SurvivalRank_Categories_pvals_Q_*.{pdf,svg,png}`
- `SurvivalRank_Categories_pvals_Padj_*.{pdf,svg,png}`

---

## G) Visualization and QC

### Heatmaps of Significant Features (Script 04)

**Script**: `04_heatmaps_sigRNAs.R`

Clustered heatmaps for differentially expressed features:

1. **Feature selection**: FDR < threshold per comparison.
2. **Normalization**: Z-score transformation of logCPM (TMM) values per feature.
3. **Visualization**: ComplexHeatmap with hierarchical clustering.

**Outputs**: Heatmap figures and feature lists.

### QC Plots per Significant Feature (Script 05)

**Script**: `05_deg_qc_plots.R`

Individual feature-level quality control:

1. **Plot types**: Box/jitter plots for categorical variables, scatter plots for continuous variables.
2. **Feature selection**: Top features ranked by FDR.

**Outputs**: Individual PNGs and optional summary PDF.

---

## H) Software Environment

### Programming Language

- **R**: version 4.4.3 (2025-02-28)

### Core Packages

| Package | Version | Purpose |
|---------|---------|---------|
| edgeR | 4.4.2 | Differential expression (quasi-likelihood) |
| limma | 3.62.2 | Linear modeling, MDS |
| data.table | 1.17.8 | High-performance data manipulation |
| ComplexHeatmap | 2.22.0 | Publication-quality heatmaps |
| circlize | 0.4.17 | Circular visualization, color functions |
| survival | 3.8.3 | Cox regression, Kaplan-Meier |
| rbioapi | 0.8.3 | miEAA API access |
| ggplot2 | 4.0.1 | Grammar of graphics visualization |
| patchwork | 1.3.2 | Multi-panel figure assembly |

### System

- **Operating System**: Linux (WSL2)
- **Conda Environment**: `omics-R`

---

## I) Reproducibility

### Logging

All scripts generate timestamped logs in `logs/` with the format `<script>_<timestamp>.txt`. Logs contain:

- Execution timestamp
- Input file paths
- Parameter values
- Processing steps and counts
- Output file paths

### Output Organization

| Directory | Contents |
|-----------|----------|
| `results/figures/` | Visualization outputs (PDF, SVG, PNG) |
| `results/tables/` | Tabular results (TSV, CSV) |
| `results/DE*/` | Differential expression results |
| `data/processed/` | Processed intermediate files |

### Random Seed

Where applicable (e.g., visualization jitter), random seeds are configurable via `--seed` parameter (default: 42). Scripts explicitly call `set.seed()` to ensure reproducibility.

### Reproducibility Checklist

- [ ] Confirm R and package versions match documented versions
- [ ] Archive `config/de_specs.csv` used for the analysis run
- [ ] Record exact paths to input count matrix and metadata
- [ ] Preserve timestamped logs from each execution
- [ ] Verify differential expression tables are from the intended run (check file modification times)

---

## J) Example Commands

### Quality Control and Preprocessing

```bash
# Initial count inspection
Rscript scripts/01_inspeccion_counts.R data/intermediate/Gliomas_all_counts_merged.csv

# LogCPM normalization and MDS
Rscript scripts/02_logcpm_mds.R data/intermediate/Gliomas_all_counts_merged.csv
```

### Differential Expression

```bash
# Multi-comparison DE (via spec file)
Rscript scripts/03_edgeR_multiDE.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --outdir results/DE

# miRNA + protein-coding restricted DE
Rscript scripts/03_2_edgeR_multiDE_miRNA_protCoding.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --outdir results/DE_miRNA_protCoding
```

### Visualization

```bash
# Heatmaps
Rscript scripts/04_heatmaps_sigRNAs.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE

# QC plots
Rscript scripts/05_deg_qc_plots.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE
```

### Survival Analysis

```bash
# Exploratory Spearman correlation
Rscript scripts/06_survival_exploratory_spearman.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE \
  --fdr 0.1

# Cox univariate regression
Rscript scripts/07_survival_cox_univariate.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE \
  --fdr 0.1

# Kaplan-Meier curves
Rscript scripts/08_KM_DE_candidates.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --fdr_cut 0.1
```

### Functional Enrichment (Differential Expression-based)

```bash
# miEAA GSEA
Rscript scripts/09_miEAA_GSEA_all_comparisons.R \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --out_root results/tables/MiEAA_GSEA

# Bubble plots
Rscript scripts/10_miEAA_GSEA_bubble_plots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_bubble \
  --run_tag A_conservative

# QC barplots
Rscript scripts/11_miEAA_GSEA_barplots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_QC \
  --run_tag A_conservative \
  --cutoff 0.25

# Category p-value distributions
Rscript scripts/12_miEAA_GSEA_categories_pvals.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_categories_pvals \
  --run_tag A_conservative
```

### Functional Enrichment (Survival-based)

```bash
# Survival-ranked GSEA
Rscript scripts/13_survGSEA.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --out_root results \
  --run_tag SurvivalRank_CoxZ_miRPathDB

# Survival GSEA visualization
Rscript scripts/14_survGSEA_plots.R \
  --in_root results/tables/SurvivalRank_GSEA \
  --out_root results/figures/SurvivalRank_GSEA \
  --run_tag SurvivalRank_CoxZ_miRPathDB \
  --n_mirpathdb 12 \
  --preset double_col
```

---

## K) Key Parameters Reference

| Parameter | Description | Default |
|-----------|-------------|---------|
| `drop_regex` | Regex for sample exclusion (post-mortem) | `^A` |
| `min_libsize` | Minimum library size threshold | 100,000 |
| `prior_count` | Prior count for logCPM calculation | 2 |
| `min_expr` | Minimum total count for feature retention | 10 |
| `fdrs` | FDR thresholds for visualization | 0.05, 0.10 |
| `fdr_thr` | FDR threshold for feature selection (survival) | 0.10 |
| `max_features` | Maximum features for heatmap/QC | 200 / 30 |
| `rank_mode` | GSEA ranking metric | `signed_sqrtF` |
| `mieaa_categories` | miRPathDB categories | GO BP, GO MF, KEGG, Reactome |
| `bubble_n_mirpathdb` | Top N terms per database | 12 |
| `qc_cutoff` | Significance cutoff for QC plots | 0.25 |
| `km_cut` | Kaplan-Meier group stratification | median |
| `min_group_n` | Minimum samples per KM group | 5 |
| `scale_expr` | Standardize expression in Cox models | TRUE |
| `preset` | Figure dimension preset | `double_col` |
| `seed` | Random seed for reproducibility | 42 |

---

## L) Script Summary Table

| Script | Purpose | Key Outputs |
|--------|---------|-------------|
| `01_inspeccion_counts.R` | Count matrix QC, library sizes | Inspection figures, RDS |
| `02_logcpm_mds.R` | Normalization, MDS | logCPM matrix, MDS coordinates |
| `03_edgeR_multiDE.R` | Differential expression (all features) | DE tables, volcano/MD plots |
| `03_2_edgeR_multiDE_miRNA_protCoding.R` | DE (miRNA + protein-coding) | DE tables, plots |
| `04_heatmaps_sigRNAs.R` | Clustered heatmaps | Heatmap PDFs, feature lists |
| `05_deg_qc_plots.R` | Feature-level QC | Individual PNGs, summary PDF |
| `06_survival_exploratory_spearman.R` | Expression-survival correlation | Correlation tables, scatter plots |
| `07_survival_cox_univariate.R` | Cox regression | HR tables with CI and FDR |
| `08_KM_DE_candidates.R` | Kaplan-Meier analysis | KM curves, summary tables |
| `09_miEAA_GSEA_all_comparisons.R` | GSEA (DE-based) | GSEA tables, ranked lists |
| `10_miEAA_GSEA_bubble_plots.R` | Bubble plot visualization | PDF/SVG/PNG figures |
| `11_miEAA_GSEA_barplots.R` | QC barplots | PDF/SVG/PNG figures |
| `12_miEAA_GSEA_categories_pvals.R` | Category distributions | PDF/SVG/PNG figures |
| `13_survGSEA.R` | GSEA (survival-ranked) | Cox tables, GSEA results |
| `14_survGSEA_plots.R` | Survival GSEA visualization | PDF/SVG/PNG figures |
