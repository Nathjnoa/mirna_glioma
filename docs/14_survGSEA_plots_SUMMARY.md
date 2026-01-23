# 🎉 SurvivalRank GSEA Plots - Resumen de Mejoras

**Fecha**: 2026-01-22
**Script**: `scripts/14_survGSEA_plots.R`
**Estado**: ✅ **COMPLETADO Y PROBADO**

---

## 📊 Resultados de Ejecución

### Figuras Generadas

Cada ejecución genera **5 tipos de figuras**:

1. **Enriched** (🟢 verde): Procesos enriquecidos (top N por DB)
2. **Depleted** (🟠 naranja): Procesos depletados (top N por DB)
3. **Combined** (🟢 + 🟠): Ambos lado a lado en una sola figura
4. **Categories Q-value**: Boxplot de distribución de Q-values por categoría
5. **Categories P-adjusted**: Boxplot de distribución de P-adjusted por categoría

### Formatos de Export

Cada figura se exporta en **3 formatos**:

- **PDF** (vector, 300 dpi): Para journals y publicaciones
- **SVG** (vector, editable): Para edición en Illustrator/Inkscape
- **PNG** (raster, 600 dpi): Para presentaciones y supplements

### Total de Archivos Generados por Ejecución

```
✅  3 bubble plots individuales PDF (enriched + depleted + combined)
✅  3 bubble plots individuales SVG (enriched + depleted + combined)
✅  3 bubble plots individuales PNG (enriched + depleted + combined)
✅  2 categories plots PDF (Q-value + P-adjusted)
✅  2 categories plots SVG (Q-value + P-adjusted)
✅  2 categories plots PNG (Q-value + P-adjusted)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  15 archivos totales por ejecución
```

---

## 🆕 Mejoras Implementadas

### Basadas en Script 10 (miEAA Bubble Plots)

✅ **Separación enriched/depleted**: Dos bubble plots independientes por dirección
✅ **Lint-harden mejorado**: Validación de inputs, logging estructurado, reproducibilidad (seed)
✅ **Publication-ready presets**: Dimensiones explícitas (180x120 mm, double_col) con DPI explícitos
✅ **Paleta colorblind-safe**: Okabe-Ito (verde #009E73, naranja #D55E00)
✅ **Exports determinísticos**: PDF + SVG + PNG con dimensiones y DPI explícitos
✅ **SVG support**: svglite para exports vectoriales editables
✅ **Tamaño de letra aumentado**: `axis.text.y` usa `p$axis` (8.5pt) en lugar de `p$ticks` (7.5pt)
✅ **Figuras combinadas**: patchwork para lado a lado enriched + depleted

### Mejoras Adicionales Específicas

✅ **Parámetro seed**: `--seed` (default: 42) para reproducibilidad de jitter
✅ **Categories plots mejorados**: Color scheme consistente con Okabe-Ito
✅ **Logging detallado**: Exports con dimensiones y DPI reportados
✅ **Validación de directorio**: Verifica que `in_root` existe antes de procesar

---

## 📁 Estructura de Archivos

### Nomenclatura

```
# Bubble plots individuales
SurvivalRank_Bubble_miRPathDB_enriched_<label>_<timestamp>.{pdf,svg,png}
SurvivalRank_Bubble_miRPathDB_depleted_<label>_<timestamp>.{pdf,svg,png}

# Bubble plot combinado
SurvivalRank_Bubble_miRPathDB_combined_<label>_<timestamp>.{pdf,svg,png}

# Categories p-values
SurvivalRank_Categories_pvals_Q_<label>_<timestamp>.{pdf,svg,png}
SurvivalRank_Categories_pvals_Padj_<label>_<timestamp>.{pdf,svg,png}
```

### Ejemplo de Archivos Generados

```
results/figures/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/
├── SurvivalRank_Bubble_miRPathDB_enriched_..._20260122_192312.pdf   (29 KB)
├── SurvivalRank_Bubble_miRPathDB_enriched_..._20260122_192312.svg   (33 KB)
├── SurvivalRank_Bubble_miRPathDB_enriched_..._20260122_192312.png   (1.2 MB)
├── SurvivalRank_Bubble_miRPathDB_depleted_..._20260122_192312.pdf   (28 KB)
├── SurvivalRank_Bubble_miRPathDB_depleted_..._20260122_192312.svg   (30 KB)
├── SurvivalRank_Bubble_miRPathDB_depleted_..._20260122_192312.png   (971 KB)
├── SurvivalRank_Bubble_miRPathDB_combined_..._20260122_192312.pdf   (34 KB)
├── SurvivalRank_Bubble_miRPathDB_combined_..._20260122_192312.svg   (63 KB)
├── SurvivalRank_Bubble_miRPathDB_combined_..._20260122_192312.png   (1.9 MB)
├── SurvivalRank_Categories_pvals_Q_..._20260122_192312.pdf          (118 KB)
├── SurvivalRank_Categories_pvals_Q_..._20260122_192312.svg          (405 KB)
├── SurvivalRank_Categories_pvals_Q_..._20260122_192312.png          (969 KB)
├── SurvivalRank_Categories_pvals_Padj_..._20260122_192312.pdf       (119 KB)
├── SurvivalRank_Categories_pvals_Padj_..._20260122_192312.svg       (405 KB)
└── SurvivalRank_Categories_pvals_Padj_..._20260122_192312.png       (973 KB)
```

---

## 🎨 Características Técnicas

### Dimensiones y Styling

| Aspecto | Valor |
|---------|-------|
| **Preset** | `double_col` (journal double column) |
| **Individuales** | 225 x 240 mm (180 × 1.25 width, 120 × 2 height) |
| **Combinadas** | 450 x 240 mm (180 × 2.5 width, 120 × 2 height) |
| **Categories** | 180 x 120 mm (base preset) |
| **Base font** | 8.5 pt Helvetica |
| **Títulos** | 10 pt bold |
| **Ejes** | 8.5 pt |
| **Términos (Y)** | 8.5 pt (↑ mejorado desde 7.5 pt) |
| **Ticks (X)** | 7.5 pt |
| **Leyenda** | 7.5 pt bold title, 7.5 pt text |

### Exports y Formatos

| Formato | DPI | Device | Tamaño típico | Uso |
|---------|-----|--------|---------------|-----|
| **PDF** | 300 | cairo_pdf | 28-119 KB | Journals, publicaciones |
| **SVG** | 300 | svglite | 30-405 KB | Edición vectorial |
| **PNG** | 600 | png | 1-2 MB | Presentaciones, web |

### Paleta de Colores

```
Enriched:  #009E73 (verde Okabe-Ito)
Depleted:  #D55E00 (naranja Okabe-Ito)
Other:     #999999 (gris neutro)
Boxplot fill: #009E73 (verde, consistente)
```

✅ Colorblind-safe (deuteranopia, protanopia)
✅ Contraste óptimo para impresión B/N
✅ Colores únicos por plot (no leyenda de color necesaria en individuales)

---

## 🚀 Uso del Script

### Ejecución Básica

```bash
cd ~/bioinfo/projects/mirna_glioma
conda activate omics-R
Rscript scripts/14_survGSEA_plots.R
```

### Parámetros Disponibles

```bash
Rscript scripts/14_survGSEA_plots.R \
  --in_root results/tables/SurvivalRank_GSEA \
  --out_root results/figures/SurvivalRank_GSEA \
  --run_tag SurvivalRank_CoxZ_miRPathDB \
  --label "Survival GSEA Analysis" \
  --n_mirpathdb 15 \              # Top N términos por DB (default: 12)
  --wrap_width 35 \                # Ancho de wrap para labels (default: 40)
  --preset double_col \            # single_col | double_col | presentation | poster
  --stamp true \                   # Añadir timestamp a nombres (default: true)
  --seed 42                        # Seed para reproducibilidad (default: 42)
```

### Ejemplos de Uso

```bash
# Mostrar top 15 términos en lugar de 12
Rscript scripts/14_survGSEA_plots.R --n_mirpathdb 15

# Usar preset single column (más pequeño)
Rscript scripts/14_survGSEA_plots.R --preset single_col

# Preset presentation (para slides)
Rscript scripts/14_survGSEA_plots.R --preset presentation

# Labels más cortos (wrap a 30 caracteres)
Rscript scripts/14_survGSEA_plots.R --wrap_width 30

# Sin timestamp en nombres de archivo
Rscript scripts/14_survGSEA_plots.R --stamp false
```

---

## 📋 Comparación Antes vs Ahora

| Aspecto | Antes | Ahora |
|---------|-------|-------|
| **Plots por ejecución** | 3 (2x2 mixed, enriched-filtered, 2x categories) | 5 (enriched, depleted, combined, 2x categories) |
| **Formatos de export** | PDF + SVG + PNG | PDF + SVG + PNG (mejorados) |
| **SVG disponible** | ✅ Sí | ✅ Sí (svglite) |
| **DPI explícitos** | ❌ No (hardcoded 300/600) | ✅ Sí (en presets) |
| **Tamaño letra términos** | 7.5pt (pequeño) | 8.5pt (más legible) |
| **Figura combinada** | ❌ No | ✅ Sí (patchwork) |
| **Ancho figura combinada** | N/A | 450mm (2.5x base) |
| **Colorblind-safe** | ✅ Sí | ✅ Sí (sin cambios) |
| **Reproducibilidad** | ❌ No | ✅ Sí (seed configurable) |
| **Logging detallado** | Básico | Estructurado con dimensiones/DPI |
| **Validación inputs** | Mínima | Mejorada (directorios, parámetros) |
| **Categories plots** | Color azul | Color verde Okabe-Ito (consistente) |
| **Total archivos** | 9 (3 plots × 3 formatos) | 15 (5 plots × 3 formatos) |

---

## 🔧 Cambios Específicos en el Código

### 1. Añadido parámetro seed y validación

```r
seed <- as.integer(get_arg("--seed", "42"))

# Input validation
if (!nzchar(run_tag)) stop("run_tag cannot be empty")
if (n_mirpathdb < 1) stop("n_mirpathdb must be >= 1")
if (wrap_width < 10) stop("wrap_width must be >= 10")
```

### 2. Importado patchwork y set seed

```r
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)  # For multi-panel combined figures
})

set.seed(seed)
```

### 3. Presets con DPI explícitos

```r
double_col = list(
  width_mm = 180, height_mm = 120, base = 8.5, title = 10, axis = 8.5,
  legend = 7.5, ticks = 7.5, line = 0.7, point = 1.8, spacing_mm = 2.5,
  margin_mm = 7, dpi_vector = 300, dpi_raster = 600
)
```

### 4. Theme mejorado con axis.text.y aumentado

```r
axis.text.y = element_text(size = p$axis, color = "black"),  # 8.5pt instead of 7.5pt
```

### 5. save_plot mejorado con logging

```r
save_plot <- function(plot, out_base, n_rows = 1, width_mult = 1, height_mult = 1) {
  width_mm <- p$width_mm * width_mult
  height_mm <- p$height_mm * height_mult * n_rows

  # ... exports with detailed logging ...
  cat("  Exported:", pdf_path, sprintf("(%d x %d mm)\n", round(width_mm), round(height_mm)))
}
```

### 6. Plots separados por dirección

```r
for (dir_type in c("Enriched", "Depleted")) {
  df_dir <- df_mirpath %>%
    filter(direction == dir_type) %>%
    group_by(db) %>%
    arrange(pval_use, .by_group = TRUE) %>%
    slice_head(n = n_mirpathdb) %>%
    ungroup()

  # ... create and save individual plot ...
  plots_list[[dir_type]] <- p_dir
}
```

### 7. Figura combinada con patchwork

```r
p_combined <- plots_list[["Enriched"]] + plots_list[["Depleted"]] +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("miRPathDB Enrichment Summary: ", plot_label),
    theme = theme(plot.title = element_text(size = p$title + 1, face = "bold", hjust = 0.5))
  )

save_plot(p_combined, out_base_combined, n_rows = 2, width_mult = 2.5)
```

---

## 📦 Dependencias

### R packages (todos instalados en omics-R)

```r
library(ggplot2)      # Plotting
library(dplyr)        # Data manipulation
library(stringr)      # String wrapping
library(patchwork)    # Multi-panel figures (ya estaba instalado)
library(svglite)      # SVG exports (ya estaba instalado)
```

---

## ✅ Estado Final

**Todo completado exitosamente**:

- ✅ Separación de gráficas enriched vs depleted
- ✅ Tamaño de letra de términos aumentado (8.5pt)
- ✅ Figuras combinadas implementadas con patchwork
- ✅ DPI explícitos en presets (300 vector, 600 raster)
- ✅ Seed configurable para reproducibilidad
- ✅ Logging estructurado mejorado
- ✅ save_plot con dimensiones y DPI reportados
- ✅ Validación de inputs mejorada
- ✅ Categories plots con color scheme consistente
- ✅ 15 archivos generados (5 plots × 3 formatos)
- ✅ Script probado y funcionando

**Listo para publicación** 🎉

---

## 📚 Documentación Relacionada

- **Script original**: `scripts/14_survGSEA_plots.R`
- **Script de referencia**: `scripts/10_miEAA_GSEA_bubble_plots.R`
- **Documentación script 10**: `docs/10_bubble_plots_FINAL_SUMMARY.md`
- **Logs de ejecución**: `logs/survGSEA_plots_*.txt`
