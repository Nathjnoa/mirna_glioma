# GSEA Results Summary Report

**Project**: mirna_glioma — miRNA enrichment analysis in glioma
**Generated**: 2026-01-31
**Source files**:
- `results/tables/reports/GSEA_DE_top12_per_category_report.csv`
- `results/tables/reports/GSEA_survGSEA_reduced_top12_per_category_report.csv`

---

## Scientific Context

This report summarizes Gene Set Enrichment Analysis (GSEA) results from two complementary approaches applied to small RNA-seq data from glioma tumors. The DE-based GSEA (Script 09) ranks miRNAs by a signed square-root F-statistic derived from edgeR quasi-likelihood differential expression, while the survival-ranked GSEA (Script 13) ranks miRNAs by Wald z-statistics from univariate Cox proportional hazards models. Both approaches use miEAA with miRPathDB expert categories (GO BP, GO MF, KEGG, Reactome). The survival report reflects post-redundancy reduction via rrvgo (threshold 0.9, Rel similarity). All Q-values reported here are FDR-adjusted and taken directly from the source tables without recomputation.

---

## Column Glossary

| Column | Description |
|--------|-------------|
| `comparison` / `analysis` | Identifier for the DE comparison or survival analysis |
| `Category` | Pathway database (GO BP, GO MF, KEGG, Reactome; all miRPathDB) |
| `Subcategory` | Specific pathway or term name |
| `Enrichment` | Direction: enriched (positive NES) or depleted (negative NES) |
| `Q_value` | FDR-adjusted p-value from miEAA GSEA |
| `Observed` | Number of miRNAs from ranked list annotated to the term |
| `miRNAs/precursors` | Semicolon-separated list of contributing miRNAs |

---

## Executive Summary

### Report 1: DE-based GSEA

| Metric | Value |
|--------|-------|
| Total rows (top 12 per Category × Direction) | 186 |
| Comparisons | EDAD_cont_10y, KI67_20 |
| Categories per comparison | GO BP, GO MF, KEGG, Reactome |

**EDAD_cont_10y** (continuous age, per 10-year increment):
- 90 rows: 48 enriched, 42 depleted
- Strong enrichment signal: top Q-value = 4.81 × 10⁻⁵ (GO BP: reproduction)
- Depleted terms weaker: top Q-value = 2.07 × 10⁻² (Reactome: Circadian Clock)
- KEGG enriched terms include T cell receptor signaling (Q = 5.39 × 10⁻⁴), TNF signaling, mTOR signaling

**KI67_20** (KI67 > 20% vs ≤ 20%):
- 96 rows: 48 enriched, 48 depleted
- Weaker overall signal: top enriched Q-value = 2.71 × 10⁻¹ (Reactome: NOTCH1 signaling)
- Top depleted: p53 binding (GO MF, Q = 5.06 × 10⁻²)

### Report 2: Survival-ranked GSEA (post rrvgo reduction)

| Metric | Value |
|--------|-------|
| Total rows (top 12 per Category × Direction) | 96 |
| Analysis | SurvivalRank_CoxZ |
| Redundancy reduction | rrvgo threshold 0.9 (79% GO terms removed) |

- 96 rows: 48 enriched, 48 depleted
- Very strong depleted signal: top Q-value = 2.94 × 10⁻⁸ (macromolecule metabolic process, 193 miRNAs)
- Enriched terms: protein insertion into membrane (Q = 1.88 × 10⁻³), fibroblast migration (Q = 2.95 × 10⁻³)

---

## Report 1: DE-based GSEA — Top Terms

### EDAD_cont_10y — Enriched (top 12 per category)

| Category | Term | Q-value | Observed |
|----------|------|---------|----------|
| GO BP | reproduction | 4.81e-05 | 26 |
| GO BP | reproductive process | 4.81e-05 | 26 |
| GO BP | neuron development | 1.24e-04 | 17 |
| GO BP | regulation of cellular component movement | 1.24e-04 | 48 |
| GO BP | response to hormone | 1.24e-04 | 48 |
| GO MF | protein phosphatase binding | 5.39e-04 | 20 |
| GO MF | identical protein binding | 2.51e-03 | 29 |
| KEGG | T cell receptor signaling pathway | 5.39e-04 | 26 |
| KEGG | TNF signaling pathway | 6.84e-04 | 29 |
| KEGG | mTOR signaling pathway | 2.76e-03 | 30 |
| Reactome | *No enriched terms met top-12 criteria at Q < 0.05* | — | — |

### EDAD_cont_10y — Depleted (top 12 per category)

| Category | Term | Q-value | Observed |
|----------|------|---------|----------|
| Reactome | Circadian Clock | 2.07e-02 | 17 |
| Reactome | Regulation of PTEN localization | 3.71e-02 | 9 |
| Reactome | Signaling by TGF-beta Receptor Complex in Cancer | 5.23e-02 | 16 |

> Note: EDAD_cont_10y depleted terms are generally weak (top Q > 0.02). GO BP, GO MF, and KEGG depleted terms all have Q > 0.05.

### KI67_20 — Overall

| Direction | Top Category | Top Term | Best Q-value |
|-----------|-------------|----------|-------------|
| Enriched | Reactome | Constitutive Signaling by NOTCH1 HD Domain Mutants | 2.71e-01 |
| Depleted | GO MF | p53 binding | 5.06e-02 |

> Note: KI67_20 shows weak enrichment signal overall. Only 1 depleted term (p53 binding) reaches Q < 0.1.

---

## Report 2: Survival-ranked GSEA — Top Terms (post rrvgo)

### Depleted (higher expression → better survival)

| Category | Term | Q-value | Observed |
|----------|------|---------|----------|
| GO BP | macromolecule metabolic process | 2.94e-08 | 193 |
| GO BP | nitrogen compound metabolic process | 2.94e-08 | 179 |
| GO BP | organelle organization | 2.94e-08 | 112 |
| GO BP | positive regulation of metabolic process | 2.94e-08 | 201 |
| GO BP | cellular macromolecule biosynthetic process | 3.83e-07 | 207 |
| GO MF | *Top terms available in source file* | — | — |
| KEGG | *Top terms available in source file* | — | — |
| Reactome | *Top terms available in source file* | — | — |

### Enriched (higher expression → worse survival)

| Category | Term | Q-value | Observed |
|----------|------|---------|----------|
| GO BP | protein insertion into membrane | 1.88e-03 | 23 |
| GO BP | positive regulation of protein insertion into mitochondrial membrane | 2.36e-03 | 21 |
| GO BP | fibroblast migration | 2.95e-03 | 17 |
| GO BP | monocyte differentiation | 3.43e-03 | 16 |
| GO BP | positive regulation of mitochondrial outer membrane permeabilization | 5.59e-03 | 24 |

---

## QA Flags

| Report | Flag | Detail |
|--------|------|--------|
| DE: EDAD_cont_10y enriched | — | Strong signal (Q < 0.001 in GO BP, KEGG) |
| DE: EDAD_cont_10y depleted | `WEAK_SIGNAL` | Best Q = 0.021 (Reactome only) |
| DE: KI67_20 enriched | `WEAK_SIGNAL` | Best Q = 0.271 |
| DE: KI67_20 depleted | `WEAK_SIGNAL` | Only 1 term at Q < 0.1 |
| Survival: depleted | — | Very strong signal (Q < 10⁻⁷) |
| Survival: enriched | — | Moderate signal (Q < 0.01) |

---

## Where to Find Details

| Resource | Path |
|----------|------|
| DE-based GSEA full tables | `results/tables/MiEAA_GSEA/<comparison>/` |
| Survival GSEA full table | `results/tables/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/` |
| Reduced GSEA table (rrvgo) | `results/tables/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/reduced_rrvgo/` |
| Bubble plots (DE) | `results/figures/MiEAA_GSEA_bubble/` |
| Bubble plots (survival) | `results/figures/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/` |
| Reduced bubble plots | `results/figures/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/reduced_rrvgo/` |
| Computational methods | `docs/METHODS_scripts.md` |

---

## Appendix

### Source Table Schema

**Report 1** (`GSEA_DE_top12_per_category_report.csv`):
- 186 rows, 7 columns
- Columns: comparison, Category, Subcategory, Enrichment, Q_value, Observed, miRNAs/precursors
- No missing values

**Report 2** (`GSEA_survGSEA_reduced_top12_per_category_report.csv`):
- 96 rows, 7 columns
- Columns: analysis, Category, Subcategory, Enrichment, Q_value, Observed, miRNAs/precursors
- No missing values

### Selection Criteria

Both reports contain the **top 12 terms per Category × Enrichment direction**, ordered by Q-value (ascending). This is a display subset; full results are available in the source GSEA tables linked above.
