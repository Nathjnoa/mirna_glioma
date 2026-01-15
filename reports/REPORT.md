# Results Summary Report

Project: mirna_glioma
Generated: 2026-01-15 13:29

Source files:
- results/DE/DE_summary_all_20260108_114709.tsv
- config/de_specs.csv
- docs/METHODS_scripts.md
- README.md
- results/DE/*

## Scientific context (descriptive)
Este reporte resume comparaciones de expresion diferencial (DE) de miRNA en glioma, basadas en la configuracion de `config/de_specs.csv` y los scripts descritos en `docs/METHODS_scripts.md`.
El pipeline utiliza edgeR con normalizacion TMM, filtrado con `filterByExpr` y pruebas QL (glmQLFit/glmQLFTest) para generar tablas completas y subconjuntos por FDR.
Las comparaciones pueden ser binarias, por rango, continuas o de supervivencia, segun el modo especificado en la configuracion.
Este documento no recalcula estadisticos; solo reorganiza y describe los resultados existentes y enlaza salidas verificadas.
Versiones de software y entorno no se encuentran en los resultados resumidos; se mantienen los placeholders definidos en la documentacion del proyecto.

## Executive summary
- Comparaciones en el resumen: 9
- Comparaciones con algun hit (FDR < 0.05): 2
- Comparaciones con algun hit (FDR < 0.10): 4
- Comparaciones sin hits (FDR 0.05 y 0.10 = 0): 5
- Comparaciones con cero features testeadas: 0
- Rango n_samples: min 12, mediana 21, max 23
- Rango n_features_tested: min 1309, mediana 1409, max 1571
- Rango n_sig_FDR0.05: min 0, mediana 0, max 7
- Rango n_sig_FDR0.1: min 0, mediana 0, max 9
- Rango cut (si aplica): min 10, mediana 20, max 24

## Column glossary
| column | description |
|---|---|
| analysis_id | Identifier for each comparison (from config/de_specs.csv). |
| var | Metadata variable used in the comparison (e.g., clinical factor). |
| mode | Comparison mode (as configured in the DE spec). |
| case | Case/group label used for the contrast (if applicable). |
| ref | Reference/control label used for the contrast (if applicable). |
| cut | Numeric cutoff applied for binary/survival/continuous modes (if applicable). |
| cut_type | Cut rule (e.g., gt/ge/per) when a cutoff is used. |
| collapse_map | Mapping rules for collapsing levels before DE (if applicable). |
| exclude_rule | Rule for excluding samples/levels (if applicable). |
| event_var | Event indicator variable for survival mode (if applicable). |
| n_samples | Number of samples included in the comparison. |
| n_features_tested | Number of features tested after filtering. |
| groups | Group counts as reported by the DE run. |
| n_sig_FDR0.05 | Count of features with FDR < 0.05 (as reported). |
| n_sig_FDR0.1 | Count of features with FDR < 0.10 (as reported). |
| comparison | Label of the contrast in the DE output. |

## Main tables by comparison
### EDAD_40_60_vs_gt70
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| EDAD_40_60_vs_gt70 | EDAD | range_vs | gt70 | age40_60 | NA |  | 40-60=age40_60;71-200=gt70 | drop_unmapped |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 12 | 1319 | age40_60:6;gt70:6 | 0 | 0 | gt70 vs age40_60 | MISSING_FIELDS;ZERO_HITS | all: results/DE/EDAD_40_60_vs_gt70/DE_EDAD_40_60_vs_gt70_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/EDAD_40_60_vs_gt70/DE_EDAD_40_60_vs_gt70_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/EDAD_40_60_vs_gt70/DE_EDAD_40_60_vs_gt70_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/EDAD_40_60_vs_gt70/QC_EDAD_40_60_vs_gt70_summary_20260108_114709.tsv<br>mdplot: results/DE/EDAD_40_60_vs_gt70/MDplot_EDAD_40_60_vs_gt70_20260108_114709.png<br>volcano: results/DE/EDAD_40_60_vs_gt70/Volcano_EDAD_40_60_vs_gt70_gt70_vs_age40_60_20260108_114709.png |

### EDAD_cont_10y
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| EDAD_cont_10y | EDAD | continuous |  |  | 10 |  | NA | drop_na |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 23 | 1571 | all:23 | 4 | 4 | continuous_per10 | MISSING_FIELDS | all: results/DE/EDAD_cont_10y/DE_EDAD_cont_10y_continuous_per10_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/EDAD_cont_10y/DE_EDAD_cont_10y_continuous_per10_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/EDAD_cont_10y/DE_EDAD_cont_10y_continuous_per10_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/EDAD_cont_10y/QC_EDAD_cont_10y_summary_20260108_114709.tsv<br>mdplot: results/DE/EDAD_cont_10y/MDplot_EDAD_cont_10y_continuous_per10_20260108_114709.png<br>volcano: results/DE/EDAD_cont_10y/Volcano_EDAD_cont_10y_continuous_per10_20260108_114709.png |

### GENERO
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| GENERO | GENERO | as_is | 1 | 0 | NA |  |  |  |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 23 | 1309 | 0:12;1:11 | 0 | 0 | 1 vs 0 | MISSING_FIELDS;ZERO_HITS | all: results/DE/GENERO/DE_GENERO_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/GENERO/DE_GENERO_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/GENERO/DE_GENERO_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/GENERO/QC_GENERO_summary_20260108_114709.tsv<br>mdplot: results/DE/GENERO/MDplot_GENERO_20260108_114709.png<br>volcano: results/DE/GENERO/Volcano_GENERO_1_vs_0_20260108_114709.png |

### KI67_>20
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| KI67_>20 | KI67_% | binary_cut | gt20 | le20 | 20 | gt |  | drop_na |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 21 | 1503 | le20:15;gt20:6 | 7 | 9 | gt20 vs le20 | MISSING_FIELDS | all: results/DE/KI67_20/DE_KI67_20_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/KI67_20/DE_KI67_20_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/KI67_20/DE_KI67_20_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/KI67_20/QC_KI67_20_summary_20260108_114709.tsv<br>mdplot: results/DE/KI67_20/MDplot_KI67_20_20260108_114709.png<br>volcano: results/DE/KI67_20/Volcano_KI67_20_gt20_vs_le20_20260108_114709.png |

### OLIG__1_vs_0
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| OLIG__1_vs_0 | OLIG_2 | collapse_levels | 1 | 0 | NA |  | 0=0;1=1 | drop_unmapped |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 16 | 1429 | 0:6;1:10 | 0 | 2 | 1 vs 0 | MISSING_FIELDS;FDR0.1_ONLY | all: results/DE/OLIG_1_vs_0/DE_OLIG_1_vs_0_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/OLIG_1_vs_0/DE_OLIG_1_vs_0_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/OLIG_1_vs_0/DE_OLIG_1_vs_0_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/OLIG_1_vs_0/QC_OLIG_1_vs_0_summary_20260108_114709.tsv<br>mdplot: results/DE/OLIG_1_vs_0/MDplot_OLIG_1_vs_0_20260108_114709.png<br>volcano: results/DE/OLIG_1_vs_0/Volcano_OLIG_1_vs_0_1_vs_0_20260108_114709.png |

### P53_TP53_1_vs_0
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| P53_TP53_1_vs_0 | P53_(TP53) | collapse_levels | 1 | 0 | NA |  | 0=0;1=1 | drop_unmapped |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 19 | 1409 | 0:8;1:11 | 0 | 0 | 1 vs 0 | MISSING_FIELDS;ZERO_HITS | all: results/DE/P53_TP53_1_vs_0/DE_P53_TP53_1_vs_0_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/P53_TP53_1_vs_0/DE_P53_TP53_1_vs_0_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/P53_TP53_1_vs_0/DE_P53_TP53_1_vs_0_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/P53_TP53_1_vs_0/QC_P53_TP53_1_vs_0_summary_20260108_114709.tsv<br>mdplot: results/DE/P53_TP53_1_vs_0/MDplot_P53_TP53_1_vs_0_20260108_114709.png<br>volcano: results/DE/P53_TP53_1_vs_0/Volcano_P53_TP53_1_vs_0_1_vs_0_20260108_114709.png |

### RESECCION_01_vs_2
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| RESECCION_01_vs_2 | RADIOLOGICAL_RESECTION | collapse_levels | 2 | 01 | NA |  | 0|1=01;2=2 | drop_na |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 21 | 1366 | 01:12;2:9 | 0 | 0 | 2 vs 01 | MISSING_FIELDS;ZERO_HITS | all: results/DE/RESECCION_01_vs_2/DE_RESECCION_01_vs_2_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/RESECCION_01_vs_2/DE_RESECCION_01_vs_2_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/RESECCION_01_vs_2/DE_RESECCION_01_vs_2_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/RESECCION_01_vs_2/QC_RESECCION_01_vs_2_summary_20260108_114709.tsv<br>mdplot: results/DE/RESECCION_01_vs_2/MDplot_RESECCION_01_vs_2_20260108_114709.png<br>volcano: results/DE/RESECCION_01_vs_2/Volcano_RESECCION_01_vs_2_2_vs_01_20260108_114709.png |

### SEGUIM_>=24m
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| SEGUIM_>=24m | MESES_SEGUIMIENTO_ | survival_cut | ge24 | lt24 | 24 | ge |  | exclude_if_censored_lt_cut | MUERTE |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 22 | 1401 | lt24:14;ge24:8 | 0 | 1 | ge24 vs lt24 | MISSING_FIELDS;FDR0.1_ONLY | all: results/DE/SEGUIM_24m/DE_SEGUIM_24m_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/SEGUIM_24m/DE_SEGUIM_24m_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/SEGUIM_24m/DE_SEGUIM_24m_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/SEGUIM_24m/QC_SEGUIM_24m_summary_20260108_114709.tsv<br>mdplot: results/DE/SEGUIM_24m/MDplot_SEGUIM_24m_20260108_114709.png<br>volcano: results/DE/SEGUIM_24m/Volcano_SEGUIM_24m_ge24_vs_lt24_20260108_114709.png |

### SINAPTOFISINA_1_vs_0
Design
| analysis_id | var | mode | case | ref | cut | cut_type | collapse_map | exclude_rule | event_var |
|---|---|---|---|---|---|---|---|---|---|
| SINAPTOFISINA_1_vs_0 | SINAPTOFISINA | collapse_levels | 1 | 0 | NA |  | 0=0;1=1 | drop_unmapped |  |
Results
| n_samples | n_features_tested | groups | n_sig_FDR0.05 | n_sig_FDR0.1 | comparison | flags | files |
|---|---|---|---|---|---|---|---|
| 12 | 1451 | 0:4;1:8 | 0 | 0 | 1 vs 0 | MISSING_FIELDS;ZERO_HITS | all: results/DE/SINAPTOFISINA_1_vs_0/DE_SINAPTOFISINA_1_vs_0_all_20260108_114709.tsv<br>sig_FDR0.05: results/DE/SINAPTOFISINA_1_vs_0/DE_SINAPTOFISINA_1_vs_0_sig_FDR0.05_20260108_114709.tsv<br>sig_FDR0.1: results/DE/SINAPTOFISINA_1_vs_0/DE_SINAPTOFISINA_1_vs_0_sig_FDR0.1_20260108_114709.tsv<br>qc: results/DE/SINAPTOFISINA_1_vs_0/QC_SINAPTOFISINA_1_vs_0_summary_20260108_114709.tsv<br>mdplot: results/DE/SINAPTOFISINA_1_vs_0/MDplot_SINAPTOFISINA_1_vs_0_20260108_114709.png<br>volcano: results/DE/SINAPTOFISINA_1_vs_0/Volcano_SINAPTOFISINA_1_vs_0_1_vs_0_20260108_114709.png |

## Where to find details
Comparisons detected under `results/DE` and whether they appear in the summary table (normalized match):
| results/DE directory | in summary |
|---|---|
| EDAD_40_60_vs_gt70 | yes |
| EDAD_66 | no |
| EDAD_cont_10y | yes |
| GENERO | yes |
| KI67_20 | yes |
| OLIG_1_vs_0 | yes |
| P53_TP53_1_vs_0 | yes |
| RESECCION_01_vs_2 | yes |
| SEGUIM_24m | yes |
| SINAPTOFISINA_1_vs_0 | yes |

## Appendix
Discovered context files:
- README.md
- docs/METHODS_scripts.md
- config/de_specs.csv

Original column list:
analysis_id, var, mode, case, ref, cut, cut_type, collapse_map, exclude_rule, event_var, n_samples, n_features_tested, groups, n_sig_FDR0.05, n_sig_FDR0.1, comparison

Missingness per column (count of empty/NA):
| column | missing |
|---|---|
| analysis_id | 0 |
| var | 0 |
| mode | 0 |
| case | 1 |
| ref | 1 |
| cut | 6 |
| cut_type | 7 |
| collapse_map | 4 |
| exclude_rule | 1 |
| event_var | 8 |
| n_samples | 0 |
| n_features_tested | 0 |
| groups | 0 |
| n_sig_FDR0.05 | 0 |
| n_sig_FDR0.1 | 0 |
| comparison | 0 |

Notes on heuristics:
- Flags are derived only from missing fields, zero-tested, and zero-hit counts; no reanalysis performed.
- File links are included only when corresponding files exist under `results/DE/<comparison>/`.
