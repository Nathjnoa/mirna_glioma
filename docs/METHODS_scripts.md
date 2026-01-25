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
| `03_edgeR_multiDE.R` | Differential expression | DE tables, volcano/MD plots |
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

---

## M) Figure Interpretation Guide

This section provides technical guidance for interpreting the publication-quality figures generated by this pipeline. Understanding these visualizations is essential for biological interpretation and manuscript preparation.

---

### M.1) Heatmaps (Script 04)

**File pattern**: `results/figures/heatmaps/<comparison>/heatmap_Z_FDR*.png`

#### What the Heatmap Shows

Heatmaps display the expression patterns of differentially expressed features (miRNAs, piRNAs, etc.) across samples. Each row represents a feature, and each column represents a sample.

#### Z-score Transformation

The values displayed are **Z-scores** (standardized values), not raw expression:

```
Z = (x - μ) / σ
```

Where:
- `x` = logCPM value for a specific feature in a specific sample
- `μ` = mean logCPM across all samples for that feature
- `σ` = standard deviation across all samples for that feature

**Interpretation of Z-score values**:

| Z-score | Interpretation |
|---------|----------------|
| +2 to +3 | Expression 2-3 standard deviations above the mean (high) |
| +1 to +2 | Expression 1-2 standard deviations above the mean (moderately high) |
| 0 | Expression at the mean level |
| -1 to -2 | Expression 1-2 standard deviations below the mean (moderately low) |
| -2 to -3 | Expression 2-3 standard deviations below the mean (low) |

**Why Z-scores?** Z-transformation normalizes each feature to have mean = 0 and SD = 1, allowing comparison of expression patterns regardless of absolute expression levels. A feature with logCPM of 15 and another with logCPM of 5 can be compared on the same scale.

#### Color Scale

The heatmap uses a divergent **blue-white-red** color palette (colorblind-safe):

| Color | Z-score | Meaning |
|-------|---------|---------|
| Dark blue (#2166AC) | -2 or lower | Strong relative downregulation |
| White | 0 | Mean expression |
| Dark red (#B2182B) | +2 or higher | Strong relative upregulation |

Values are typically clamped at ±2 or ±3 to prevent outliers from dominating the color scale.

#### Hierarchical Clustering

**Row clustering** (features):
- **Distance metric**: 1 - Pearson correlation (features with similar expression patterns cluster together)
- **Linkage method**: Complete linkage (maximizes compactness)
- **Interpretation**: Features in the same cluster show coordinated expression changes across samples

**Column ordering** (samples):
- For categorical comparisons: Samples ordered by group
- For continuous comparisons: Samples ordered by the continuous variable value
- **Top annotation bar**: Shows group membership or continuous variable value

#### Column Annotations

The colored bar above the heatmap indicates:
- **Group**: For categorical comparisons (e.g., High vs Low KI67)
- **Value**: For continuous comparisons (e.g., age gradient)

#### Reading the Heatmap

1. **Identify sample groups**: Look at the top annotation bar
2. **Find expression blocks**: Clusters of red or blue indicate coordinated up/down-regulation
3. **Compare groups**: If one group is predominantly red for certain features, those features are upregulated in that group
4. **Identify patterns**: Features clustering together may share biological functions or regulatory mechanisms

#### Example Interpretation

> "The heatmap shows 45 miRNAs with FDR < 0.05. A distinct cluster of 12 miRNAs (rows 1-12) shows strong upregulation (Z > 2, red) in the high-KI67 group compared to low-KI67. The hierarchical clustering reveals two major feature groups: one coordinately upregulated and one coordinately downregulated with proliferation status."

---

### M.2) Kaplan-Meier Survival Curves (Script 08)

**File pattern**: `results/figures/KM_DE_candidates/KM_<feature_id>_*.png`

#### What the KM Plot Shows

Kaplan-Meier curves display the probability of survival over time, stratified by expression level of a specific feature (e.g., miRNA).

#### Axes

| Axis | Label | Range | Meaning |
|------|-------|-------|---------|
| X-axis | Time (months) | 0 to max follow-up | Time since diagnosis/treatment |
| Y-axis | Survival Probability | 0 to 1.0 (or 0-100%) | Estimated proportion of patients surviving |

#### Curve Elements

1. **Step function**: Each vertical drop represents one or more death events
2. **Horizontal segments**: Periods with no events (deaths)
3. **Censored observations** (tick marks): Patients lost to follow-up or still alive at last observation
4. **Two curves**:
   - **Blue (#0072B2)**: Low expression group
   - **Orange (#D55E00)**: High expression group

#### Stratification Methods

| Method | Description | Groups |
|--------|-------------|--------|
| `median` | Split at median expression | Low (< median) vs High (≥ median) |
| `terciles_extremos` | Compare extreme tertiles | T1 (lowest 33%) vs T3 (highest 33%) |
| `cuartiles_extremos` | Compare extreme quartiles | Q1 (lowest 25%) vs Q4 (highest 25%) |

The plot title indicates the cutoff value used.

#### Statistical Test

**Log-rank test** (also called Mantel-Cox test):
- Tests the null hypothesis that survival curves are identical
- **p-value interpretation**:
  - p < 0.05: Significant difference in survival between groups
  - p < 0.001: Highly significant difference
  - p ≥ 0.05: No significant difference (curves overlap)

The p-value is displayed in the top-right corner.

#### Risk Table (if present)

Shows the number of patients "at risk" (not yet dead or censored) at each time point for each group. Declining numbers indicate events or censoring.

#### Interpretation Guidelines

**Favorable prognosis marker** (protective):
- High expression → curve above low expression curve
- High expression associated with longer survival

**Unfavorable prognosis marker** (risk factor):
- High expression → curve below low expression curve
- High expression associated with shorter survival

**Median survival**: Point where curve crosses 0.5 (50% survival probability). If curves don't reach 0.5, median survival is not reached.

#### Example Interpretation

> "hsa-miR-17-5p expression significantly stratifies patient survival (log-rank p = 0.012, n = 85). Patients with high expression (n = 43, orange) show worse prognosis with median survival of 18 months, compared to low expressors (n = 42, blue) with median survival not reached. This suggests hsa-miR-17-5p may function as an unfavorable prognostic marker."

#### Caveats

- KM analysis does not adjust for confounders (use Cox regression for multivariate analysis)
- Small group sizes (n < 10) produce unstable estimates with wide confidence intervals
- Crossing curves may indicate non-proportional hazards

---

### M.3) GSEA Bubble Plots (Scripts 10, 14)

**File pattern**: `results/figures/MiEAA_GSEA_bubble/<tag>/bubble_*.png`
**File pattern**: `results/figures/SurvivalRank_GSEA/SurvivalRank_Bubble_*.png`

#### What the Bubble Plot Shows

Bubble plots visualize pathway enrichment results from Gene Set Enrichment Analysis (GSEA). Each bubble represents a pathway or gene set.

#### Axes and Visual Encodings

| Element | Encoding | Meaning |
|---------|----------|---------|
| X-axis | -log10(Q-value) or -log10(P-adjusted) | Statistical significance (higher = more significant) |
| Y-axis | Pathway name | Ordered by significance (most significant at top) |
| Bubble size | Number of observed miRNAs | Larger = more miRNAs contributing to pathway |
| Bubble color | Fixed by direction | Green = enriched, Orange = depleted |

#### Understanding -log10(Q-value)

The X-axis displays the negative log10-transformed Q-value (FDR-adjusted p-value):

| -log10(Q) | Q-value | Interpretation |
|-----------|---------|----------------|
| 1.0 | 0.1 | Marginally significant |
| 1.3 | 0.05 | Significant (standard threshold) |
| 2.0 | 0.01 | Highly significant |
| 3.0 | 0.001 | Very highly significant |

**Key point**: Higher values on the X-axis indicate stronger statistical significance. Bubbles further to the right are more significant.

#### Enriched vs. Depleted Panels

Plots are split into two panels based on the **Normalized Enrichment Score (NES)** direction:

| Panel | Color | NES Sign | Meaning |
|-------|-------|----------|---------|
| Enriched | Green (#009E73) | Positive | Pathways active in upregulated/high-risk miRNAs |
| Depleted | Orange (#D55E00) | Negative | Pathways active in downregulated/low-risk miRNAs |

**For DE-based GSEA** (Script 10):
- **Enriched panel**: Pathways targeted by miRNAs upregulated in the comparison
- **Depleted panel**: Pathways targeted by miRNAs downregulated in the comparison

**For Survival-based GSEA** (Script 14):
- **Enriched panel**: Pathways targeted by miRNAs associated with worse survival (high expression → death)
- **Depleted panel**: Pathways targeted by miRNAs associated with better survival (high expression → protection)

#### Bubble Size: Observed miRNAs

The bubble size represents the number of miRNAs from your ranked list that are annotated to target the pathway:

| Size | Interpretation |
|------|----------------|
| Small bubble | Few miRNAs (2-5) contribute to pathway enrichment |
| Medium bubble | Moderate number (5-15) of miRNAs |
| Large bubble | Many miRNAs (>15) converge on this pathway |

Larger bubbles suggest more robust findings with multiple independent miRNAs.

#### Significance Metrics

| Metric | Description | Threshold |
|--------|-------------|-----------|
| Q-value | FDR-adjusted p-value | < 0.25 (discovery), < 0.05 (validation) |
| P-adjusted | Adjusted p-value (method varies) | < 0.05 |

#### Pathway Databases

| Database | Prefix | Description |
|----------|--------|-------------|
| GO BP | `GO:` | Biological processes |
| GO MF | `GO:` | Molecular functions |
| KEGG | `hsa:` | Metabolic and signaling pathways |
| Reactome | `R-HSA-` | Curated biological pathways |

#### Reading the Bubble Plot

1. **Check panel color**: Green = enriched (upregulated/high-risk), Orange = depleted (downregulated/low-risk)
2. **Read X-axis position**: Further right = more statistically significant
3. **Look at bubble size**: Larger = more miRNAs contributing to the enrichment
4. **Identify top pathways**: Pathways at the top of the Y-axis are the most significant

#### Example Interpretation

> "The enriched panel (green) shows 'Cell cycle regulation' as the most significant pathway (-log10(Q) = 2.5, Q = 0.003) with 18 contributing miRNAs (large bubble). This indicates that miRNAs upregulated in the high-KI67 group collectively target cell cycle pathways. In the depleted panel (orange), 'Neuronal differentiation' shows significance (-log10(Q) = 1.8, Q = 0.016) with 12 miRNAs, suggesting downregulated miRNAs normally target neuronal processes."

#### Caveats

- GSEA Q-value < 0.25 is a discovery threshold, not definitive proof
- Pathway databases have variable curation quality
- Redundant pathways may appear separately (e.g., "Cell cycle" in KEGG and Reactome)
- NES values are not displayed in the plot but are used to determine direction (enriched vs depleted)

---

### M.4) Volcano Plots (Script 03)

**File pattern**: `results/DE/<comparison>/Volcano_*.png`

#### What the Volcano Plot Shows

Volcano plots display the relationship between statistical significance (y-axis) and biological effect size (x-axis) for all features in a differential expression analysis.

#### Axes

| Axis | Value | Range | Meaning |
|------|-------|-------|---------|
| X-axis | log2 Fold Change (logFC) | -∞ to +∞ (typically -5 to +5) | Direction and magnitude of expression change |
| Y-axis | -log10(FDR) or -log10(P-value) | 0 to +∞ | Statistical significance |

#### Log2 Fold Change (logFC)

logFC represents the log2-transformed ratio of expression between conditions:

| logFC | Fold Change | Interpretation |
|-------|-------------|----------------|
| +2 | 4x higher | Expression 4 times higher in case vs reference |
| +1 | 2x higher | Expression doubled |
| 0 | No change | Equal expression |
| -1 | 2x lower | Expression halved |
| -2 | 4x lower | Expression 4 times lower in case vs reference |

**Reference direction**: In our pipeline, positive logFC means higher expression in the "case" group (e.g., High-KI67 vs Low-KI67, the first level is case).

#### Significance Threshold Lines

| Line | Value | Meaning |
|------|-------|---------|
| Horizontal dashed | -log10(0.05) ≈ 1.3 | FDR = 0.05 significance threshold |
| Vertical dotted (right) | logFC = +0.584 | 1.5-fold upregulation threshold |
| Vertical dotted (left) | logFC = -0.584 | 1.5-fold downregulation threshold |

#### Color Categories

| Color | Category | Criteria |
|-------|----------|----------|
| Orange (#D55E00) | Up-regulated | FDR < 0.05 AND logFC > 0.584 |
| Blue (#0072B2) | Down-regulated | FDR < 0.05 AND logFC < -0.584 |
| Pink (#CC79A7) | Significant (|logFC| < 0.584) | FDR < 0.05 but |logFC| ≤ 0.584 |
| Gray (#999999) | Not significant | FDR ≥ 0.05 |

#### Reading the Volcano Plot

1. **Upper-right quadrant**: Significantly upregulated features (interesting candidates)
2. **Upper-left quadrant**: Significantly downregulated features (interesting candidates)
3. **Top of plot**: Most statistically significant features (regardless of effect size)
4. **Far left/right**: Largest effect sizes (regardless of significance)
5. **Labeled points**: Top 10 most significant features (auto-labeled with feature IDs)

#### Subtitle Statistics

The plot subtitle shows counts:
- **Up**: Number of upregulated features (FDR < 0.05, logFC > 0.584)
- **Down**: Number of downregulated features (FDR < 0.05, logFC < -0.584)
- **Sig (|logFC|<0.584)**: Statistically significant but small effect size

#### Example Interpretation

> "The volcano plot for KI67_High vs KI67_Low shows 23 upregulated and 18 downregulated miRNAs (FDR < 0.05, |logFC| > 0.584). hsa-miR-17-5p shows the strongest upregulation (logFC = 2.3, FDR = 1.2e-8), while hsa-miR-137 shows the strongest downregulation (logFC = -1.9, FDR = 3.4e-6). An additional 45 features are statistically significant but with modest effect sizes (|logFC| < 0.584)."

#### Caveats

- Very high -log10(FDR) values may be truncated for visualization
- Features with extreme logFC but low significance (bottom corners) may be biologically relevant in small subsets
- The FDR threshold of 0.05 is conventional but arbitrary

---

### M.5) Color Accessibility

All figures use **colorblind-safe palettes** following the Okabe-Ito recommendations:

| Use Case | Colors | Hex Codes |
|----------|--------|-----------|
| Up/Down or High/Low | Blue vs Orange | #0072B2 vs #D55E00 |
| Enriched/Depleted | Green vs Orange | #009E73 vs #D55E00 |
| Heatmap divergent | Blue-White-Red | #2166AC → #FFFFFF → #B2182B |
| Categorical (6-level) | Okabe-Ito | #0072B2, #D55E00, #009E73, #CC79A7, #F0E442, #56B4E9 |

These palettes remain distinguishable for individuals with deuteranopia (red-green colorblindness) and protanopia.

---

### M.6) Figure Export Specifications

All publication figures are exported with the following specifications:

| Format | Resolution | Use Case |
|--------|------------|----------|
| PNG | 300-600 dpi | Journal submission, presentations |
| PDF | Vector (cairo_pdf) | Publication, post-processing |
| SVG | Vector (svglite) | Editing in Illustrator/Inkscape |

**Dimensions** follow the `double_col` journal preset:
- Base: 180 × 120 mm (7.1 × 4.7 inches)
- Font: 8.5pt base, Helvetica/sans-serif
- Line width: 0.5pt minimum

---

### M.7) Summary: Quick Reference

| Figure Type | X-axis | Y-axis | Key Metric | Significance |
|-------------|--------|--------|------------|--------------|
| Heatmap | Samples | Features | Z-score (-3 to +3) | Selected by FDR |
| KM plot | Time (months) | Survival probability (0-1) | Log-rank p-value | p < 0.05 |
| Bubble plot | -log10(Q-value) | Pathways | Q-value, bubble size | Q < 0.25 |
| Volcano | logFC | -log10(FDR) | FDR, logFC | FDR < 0.05, |logFC| > 0.584 |
