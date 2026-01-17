# Metodos computacionales (scripts)

Este documento resume lo que hacen los scripts en `scripts/` del proyecto `mirna_glioma`. El enfoque es de metodos reproducibles: se describen entradas, filtros, analisis, salidas y supuestos. Donde faltan detalles, se usan placeholders.

## A) Diseno del estudio y datos

- Tipo de datos: matriz de conteos de RNAs con columnas de muestras y columnas de anotacion (p. ej., `id`, `type`, `entrez_id`, `HGNC_symbol`).
- Entradas principales (por defecto en los scripts):
  - Conteos: `data/intermediate/Gliomas_all_counts_merged.csv`.
  - Metadatos: `data/intermediate/Metadatos_gliomas_verificados.csv`.
- Filtrado de muestras post-mortem: se excluyen columnas de muestras cuyo nombre empieza con `A` (regex `^A`).

## B) Preprocesamiento y QC

### Inspeccion inicial de conteos (`01_inspeccion_counts.R`)

- Lee la matriz de conteos, detecta columnas de anotacion y construye `counts_mat`.
- Excluye muestras con prefijo `A`.
- Coercion robusta a numerico y chequeos de integridad (NA, negativos, no enteros).
- Calcula library size con `edgeR::DGEList` y TMM (effective library size).
- Calcula fraccion de ceros por muestra y conteo total por feature.
- Genera figuras de inspeccion (histogramas de library size, fraccion de ceros, boxplot log1p).
- Agrega barplot de library size (raw counts) y MDS con todas las muestras etiquetadas.
- Guarda un RDS con objetos clave de inspeccion.

### Normalizacion logCPM y QC multivariado (`02_logcpm_mds.R`)

- Filtra muestras con `library_size < 100000` y features con:
  - conteo total > 0,
  - varianza > 0,
  - suma de conteos >= 10.
- Calcula logCPM pre y post TMM (edgeR) con `prior.count=2`.
- Genera RLE pre/post, densidades pre/post, y barplots de library size.
- Calcula MDS (limma) sobre logCPM post TMM y exporta coordenadas.
- Exporta matriz `logCPM_TMM_*.csv` y tablas de library sizes.

## C) Analisis estadistico / modelado

### Diferencial de expresion multicomparacion (`03_edgeR_multiDE.R`)

- Alinea `counts_mat` con `metadata` y exporta bases filtradas.
- Permite ejecucion por lista de variables (`--vars`) o por `spec` (`config/de_specs.csv`).
- Modos soportados (segun `spec`):
  - `as_is`, `binary_cut`, `collapse_levels`, `range_vs`, `survival_cut`, `continuous`.
- Usa `edgeR` con normalizacion TMM, `filterByExpr`, y modelo QL (`glmQLFit`/`glmQLFTest`).
- Para factores multinivel produce omnibus y comparaciones vs referencia; para binarios produce contraste case vs ref.
- Comparaciones binarias adicionales via `collapse_levels` incluyen P53_(TP53), SINAPTOFISINA y OLIG_2 (1 vs 0), excluyendo nivel 2 con `drop_unmapped`.
- Salidas por comparacion: tablas completas, FDR<0.05, FDR<0.1, MD plots y volcano plots.
- Resumen QC por comparacion y resumen global.

### Diferencial de expresion restringido a miRNA/protein_coding (`03_2_edgeR_multiDE_miRNA_protCoding.R`)

- Filtra features por biotipo (`type`) manteniendo `mirna` y `protein_coding`.
- Usa `spec` para definir comparaciones (mismos modos que arriba).
- Ejecuta `edgeR` QL con TMM y `filterByExpr`.
- Salidas y resumen global similares a `03_edgeR_multiDE.R`.

## D) Correccion por multiples pruebas

- Ajuste de FDR con Benjamini-Hochberg (BH) en los analisis de DE y en analisis exploratorios de supervivencia.
- Umbrales reportados frecuentemente: 0.05 y 0.10 (segun tablas y figuras generadas).

## E) Analisis funcional / enriquecimiento

### GSEA con miEAA (`09_miEAA_GSEA_all_comparisons.R`)

- Extrae miRNAs desde tablas DE (campo `type`) y genera ranking segun `rank_mode` (por defecto `signed_sqrtF`, i.e., `sign(logFC) * sqrt(F)`;
  alternativas: `signed_logp`, `signed_F`, `signed_logFC`).
- Ejecuta GSEA via `rbioapi::rba_mieaa_enrich` en categorias miRPathDB expert (GO BP/MF, KEGG, Reactome) y valida IDs contra `rba_mieaa_cats`.
- No aplica filtro de significancia en el request (`sig_level=1`); el filtrado se hace post-hoc.
- Exporta ranking mature, tabla completa sin filtrar (`MiEAA_GSEA_all_*.tsv`), top50 por Q-value y lista Q<0.05 por comparacion.

### Bubble plots de miEAA GSEA (`10_miEAA_GSEA_bubble_plots.R`)

- Lee la tabla completa mas reciente `MiEAA_GSEA_all_*.tsv` por comparacion (run_tag) y detecta columnas de terminos/p-values de forma robusta (prioridad: P-adjusted > Q-value > P-value).
- Para miRPathDB, selecciona top N por base (GO BP/MF, KEGG, Reactome) y genera un bubble plot 2x2 (PDF/SVG/PNG).
- Registra logs en `logs/` y copia del log en `results/figures/MiEAA_GSEA_bubble/<run_tag>/`.

### QC barplots de miEAA GSEA (`11_miEAA_GSEA_barplots.R`)

- Lee la tabla completa mas reciente `MiEAA_GSEA_all_*.tsv` por comparacion; si no existe, usa `MiEAA_GSEA_top50_*.tsv` con warning.
- Detecta columna de significancia (P-adjusted/Q-value/P-value) y direccion (enrichment).
- Genera barplots apilados por base y direccion para un cutoff (y opcionalmente un segundo cutoff).
- Exporta PDF/SVG y registra logs en `logs/` y copia en `results/figures/MiEAA_GSEA_QC/<run_tag>/`.

### Categories p-values miEAA (`12_miEAA_GSEA_categories_pvals.R`)

- Genera distribuciones de -log10(Q-value) y -log10(P-adjusted) por categoria (boxplot + jitter).
- Ordena categorias por mediana de -log10(metric) y exporta PDF/SVG.
- Usa `MiEAA_GSEA_all_*.tsv` si existe; si no, `MiEAA_GSEA_top50_*.tsv` con warning.
- Registra logs en `logs/` y copia en `results/figures/MiEAA_GSEA_categories_pvals/<run_tag>/`.

## F) Analisis de supervivencia

### Exploratorio Spearman (`06_survival_exploratory_spearman.R`)

- Construye la union de features con FDR < umbral en todas las comparaciones DE.
- Calcula correlacion Spearman entre expresion (logCPM) y tiempo de seguimiento.
- Reporta resultados para todas las muestras y solo eventos.
- Genera tabla y un set limitado de figuras.

### Cox univariado (`07_survival_cox_univariate.R`)

- Usa la union de features DE con FDR < umbral.
- Ajusta modelos de Cox univariados por feature con `survival::coxph`.
- Opcion de estandarizar la expresion (HR por 1 SD) o usar logCPM directo.
- Reporta HR, IC95%, p, FDR y tablas con membership.

### Kaplan-Meier de candidatos DE (`08_KM_DE_candidates.R`)

- Selecciona candidatos DE (FDR < umbral) de todas las comparaciones.
- Divide en alta vs baja expresion (mediana o cuantiles extremos) y evalua log-rank.
- Genera PDF con curvas KM y tabla resumen con p-values y FDR.

## G) Visualizacion y QC de DE

### Heatmaps de RNAs significativos (`04_heatmaps_sigRNAs.R`)

- Toma la tabla DE mas reciente por comparacion y selecciona features con FDR < umbral.
- Usa logCPM (TMM) y transforma a Z-score por feature.
- Genera heatmaps con `ComplexHeatmap` y guarda listas de features.

### QC por feature significativo (`05_deg_qc_plots.R`)

- Genera box/jitter (categoricos) o scatter (continuos) para top features por FDR.
- Produce PNGs y PDF opcional con resumen de QC visual.

## H) Software y entorno computacional

- Lenguaje principal: R ({R_VERSION}).
- Paquetes clave: edgeR ({EDGER_VERSION}), limma ({LIMMA_VERSION}), data.table ({DATATABLE_VERSION}),
  ComplexHeatmap ({COMPLEXHEATMAP_VERSION}), circlize ({CIRCLIZE_VERSION}), survival ({SURVIVAL_VERSION}),
  survminer ({SURVMINER_VERSION, opcional}), rbioapi ({RBIOAPI_VERSION}).
- Sistema operativo / entorno: {OS_DISTRIBUCION}, {HARDWARE_RESUMEN}.

## Reproducibilidad

- Los scripts generan logs con timestamp en `logs/`.
- Las salidas se almacenan en `results/figures/`, `results/tables/`, `results/DE*/` y `data/processed/`.
- No se fijan semillas aleatorias; si se requiere determinismo estricto, definir `set.seed({SEED})`.

### Checklist de reproducibilidad

- [ ] Confirmar versiones de R y paquetes.
- [ ] Guardar copia de `config/de_specs.csv` usada en el run.
- [ ] Registrar rutas exactas de inputs (conteos y metadatos).
- [ ] Conservar logs con timestamp del run.
- [ ] Verificar que las tablas DE usadas son las mas recientes (mtime).

### Plantilla de comandos (ejemplos)

```bash
# 01: inspeccion inicial
Rscript scripts/01_inspeccion_counts.R data/intermediate/Gliomas_all_counts_merged.csv

# 02: logCPM + MDS
Rscript scripts/02_logcpm_mds.R data/intermediate/Gliomas_all_counts_merged.csv

# 03: DE multi (por spec)
Rscript scripts/03_edgeR_multiDE.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --outdir results/DE

# 03_2: DE miRNA+protein_coding
Rscript scripts/03_2_edgeR_multiDE_miRNA_protCoding.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --outdir results/DE_miRNA_protCoding

# 04: heatmaps
Rscript scripts/04_heatmaps_sigRNAs.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE

# 05: QC plots por DEG
Rscript scripts/05_deg_qc_plots.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE

# 06: Spearman exploratorio
Rscript scripts/06_survival_exploratory_spearman.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE \
  --fdr 0.1

# 07: Cox univariado
Rscript scripts/07_survival_cox_univariate.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_dir results/DE \
  --fdr 0.1

# 08: KM para candidatos DE
Rscript scripts/08_KM_DE_candidates.R \
  --counts data/intermediate/Gliomas_all_counts_merged.csv \
  --meta data/intermediate/Metadatos_gliomas_verificados.csv \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --fdr_cut 0.1

# 09: miEAA GSEA
Rscript scripts/09_miEAA_GSEA_all_comparisons.R \
  --spec config/de_specs.csv \
  --de_root results/DE \
  --out_root results/tables/MiEAA_GSEA

# 10: bubble plots miEAA
Rscript scripts/10_miEAA_GSEA_bubble_plots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_bubble \
  --run_tag A_conservative

# 11: QC barplots miEAA
Rscript scripts/11_miEAA_GSEA_barplots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_QC \
  --run_tag A_conservative \
  --cutoff 0.25

# 12: categories p-values miEAA
Rscript scripts/12_miEAA_GSEA_categories_pvals.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_categories_pvals \
  --run_tag A_conservative
```

### Tabla de parametros clave

| Parametro | Significado | Valor (por defecto) |
|---|---|---|
| `drop_regex` | Patron para excluir muestras (post-mortem) | `^A` |
| `min_libsize` | Filtro global de library size | `100000` |
| `prior_count` | Prior para logCPM | `2` |
| `min_expr` | Filtro por suma de conteos (script 02) | `10` |
| `fdrs` | Umbrales para figuras QC/heatmaps | `0.05,0.1` |
| `fdr_thr` | Umbral union DE para supervivencia | `0.1` |
| `max_features` | Max features en heatmap/QC | `200` (heatmap), `30` (QC) |
| `rank_mode` | Score para ranking miEAA | `signed_sqrtF` |
| `mieaa_categories` | Categorias miEAA (default) | `miRPathDB_GO_Biological_process_mature; miRPathDB_GO_Molecular_function_mature; miRPathDB_KEGG_mature; miRPathDB_Reactome_mature` |
| `mieaa_sig_level` | Filtro en request miEAA | `1` |
| `mieaa_min_hits` | Min hits para miEAA | `5` |
| `bubble_n_mirpathdb` | Top N por base (miRPathDB) | `12` |
| `qc_cutoff` | Umbral de significancia para QC | `0.25` |
| `qc_also_cutoff` | Umbral adicional QC (opcional) | `0.05` |
| `categories_wrap_width` | Wrap labels categories p-values | `35` |
| `km_cut` | Regla de corte KM | `median` |
| `min_group_n` | Minimo por grupo en KM | `5` |
| `standardize` | HR por 1 SD (Cox) | `FALSE` |

## Resumen por script (trazabilidad)

- `scripts/01_inspeccion_counts.R`: inspeccion de conteos, library size, ceros, barplot de library size y MDS con etiquetas.
- `scripts/02_logcpm_mds.R`: logCPM TMM, RLE, densidades, MDS/PCA y exportes.
- `scripts/03_edgeR_multiDE.R`: DE con edgeR (QL), multiples modos via `spec` o `--vars`.
- `scripts/03_2_edgeR_multiDE_miRNA_protCoding.R`: DE restringido a miRNA y protein_coding.
- `scripts/04_heatmaps_sigRNAs.R`: heatmaps Z-score de features DE significativas.
- `scripts/05_deg_qc_plots.R`: QC visual por feature (box/jitter o scatter continuo).
- `scripts/06_survival_exploratory_spearman.R`: Spearman expr vs seguimiento en union DE.
- `scripts/07_survival_cox_univariate.R`: Cox univariado en union DE con FDR.
- `scripts/08_KM_DE_candidates.R`: KM alto/bajo para candidatos DE.
- `scripts/09_miEAA_GSEA_all_comparisons.R`: GSEA miEAA desde rankings por comparacion (default `signed_sqrtF`).
- `scripts/10_miEAA_GSEA_bubble_plots.R`: bubble plots miEAA (miRPathDB 2x2) desde tablas completas.
- `scripts/11_miEAA_GSEA_barplots.R`: barplots QC miEAA por base y direccion (cutoffs configurables).
- `scripts/12_miEAA_GSEA_categories_pvals.R`: distribuciones p-values por categoria (Q-value y P-adjusted).
