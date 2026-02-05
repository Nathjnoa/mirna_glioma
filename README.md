# mirna_glioma

Analysis of miRNA expression in glioma with quality control, differential expression, survival modeling, and functional enrichment.

## Pipeline Overview

This project implements a comprehensive analysis pipeline for small RNA sequencing data from glioma tumor samples:

1. **Quality Control** (Scripts 01-02): Count matrix inspection, TMM normalization, MDS visualization
2. **Matched Group Design** (Script 02b): Age-matched case-control groups (1:2 ratio) for confounding control
3. **Differential Expression** (Script 03): edgeR quasi-likelihood framework with multiple comparison modes
3. **Visualization** (Scripts 04-05): Heatmaps and feature-level QC plots
4. **Survival Analysis** (Script 08): Kaplan-Meier curves for DE candidates
5. **Functional Enrichment - DE-based** (Scripts 09-12): miEAA GSEA with bubble plots and QC visualizations
6. **Functional Enrichment - Survival-based** (Scripts 13-16): Survival-ranked GSEA, enrichment plots, GO redundancy reduction (rrvgo)

## Repository Structure

```
mirna_glioma/
├── config/           # Analysis parameters (de_specs.csv)
├── docs/             # Methods documentation
│   └── METHODS_scripts.md   # Detailed computational methods
├── env/              # Conda environment definition
├── scripts/          # R analysis pipeline (01-05, 08-16)
├── data/             # Input data (local, gitignored)
├── results/          # Output tables and figures (local, gitignored)
└── logs/             # Timestamped execution logs (local, gitignored)
```

## Quick Start

### 1. Environment Setup

```bash
conda env create -f env/environment.yml
conda activate omics-R
```

### 2. Configure Comparisons

Edit `config/de_specs.csv` to define differential expression comparisons.

### 3. Run Pipeline

```bash
# QC and preprocessing
Rscript scripts/01_inspeccion_counts.R data/intermediate/Gliomas_all_counts_merged.csv
Rscript scripts/02_logcpm_mds.R data/intermediate/Gliomas_all_counts_merged.csv

# Age-matched groups (optional, for specific comparisons)
Rscript scripts/02b_matching_EGFR_VIII.R --meta data/intermediate/Metadatos_gliomas_verificados.csv

# Differential expression
Rscript scripts/03_edgeR_multiDE.R --spec config/de_specs.csv

# Visualization
Rscript scripts/04_heatmaps_sigRNAs.R --spec config/de_specs.csv
Rscript scripts/05_deg_qc_plots.R --spec config/de_specs.csv

# Survival analysis
Rscript scripts/08_KM_DE_candidates.R --fdr_cut 0.1

# Functional enrichment (DE-based)
Rscript scripts/09_miEAA_GSEA_all_comparisons.R
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --run_tag A_conservative
Rscript scripts/11_miEAA_GSEA_barplots.R --run_tag A_conservative
Rscript scripts/12_miEAA_GSEA_categories_pvals.R --run_tag A_conservative

# Functional enrichment (survival-based)
Rscript scripts/13_survGSEA.R
Rscript scripts/14_survGSEA_plots.R
Rscript scripts/15_survGSEA_enrichment_plot.R --preset single_col
Rscript scripts/16_survGSEA_reduce_redundancy.R --sim_threshold 0.9
```

## Software Requirements

- **R**: 4.4.3
- **Key packages**: edgeR (4.4.2), limma (3.62.2), survival (3.8.3), rbioapi (0.8.3), rrvgo (1.18.0), ComplexHeatmap (2.22.0), ggplot2 (4.0.1), patchwork (1.3.2)
- **Environment**: `omics-R` (conda)

## Documentation

- **Computational Methods**: [docs/METHODS_scripts.md](docs/METHODS_scripts.md) — Complete pipeline documentation with parameters, outputs, and figure interpretation guide
- **Runbook**: [docs/RUNBOOK.md](docs/RUNBOOK.md) — Step-by-step reproduction instructions with checklist
- **Outputs Map**: [docs/OUTPUTS.md](docs/OUTPUTS.md) — Output directory structure and naming conventions

## Output Highlights

### Figures

- **Bubble plots**: Enriched/depleted pathways with colorblind-safe palette (PDF/SVG/PNG)
- **Heatmaps**: Z-score normalized expression with hierarchical clustering
- **Survival curves**: Kaplan-Meier with log-rank statistics
- **Combined figures**: Multi-panel layouts using patchwork

### Tables

- Differential expression results with FDR correction
- Cox regression hazard ratios with confidence intervals
- GSEA enrichment tables (full and filtered)

## Reproducibility

- All scripts generate timestamped logs in `logs/`
- Figure exports copy logs to output directories
- Random seeds configurable via `--seed` parameter (default: 42)
- Package versions documented in [docs/METHODS_scripts.md](docs/METHODS_scripts.md)

## License

MIT
