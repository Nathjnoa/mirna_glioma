# 🎉 miEAA GSEA Bubble Plots - Resumen Final de Mejoras

**Fecha**: 2026-01-22
**Script**: `scripts/10_miEAA_GSEA_bubble_plots.R`
**Estado**: ✅ **COMPLETADO Y PROBADO**

---

## 📊 Resultados de Ejecución

### Figuras Generadas por Análisis (9 análisis totales)

Cada `analysis_id` ahora genera **3 tipos de figuras**:

1. **Enriched** (🟢 verde): Procesos enriquecidos
2. **Depleted** (🟠 naranja): Procesos depletados
3. **Combined** (🟢 + 🟠): Ambos lado a lado en una sola figura

### Formatos de Export

Cada figura se exporta en **3 formatos**:

- **PDF** (vector, 300 dpi): Para journals y publicaciones
- **SVG** (vector, editable): Para edición en Illustrator/Inkscape
- **PNG** (raster, 600 dpi): Para presentaciones y supplements

### Total de Archivos Generados

```
✅ 18 figuras individuales PDF (9 enriched + 9 depleted)
✅ 18 figuras individuales SVG (9 enriched + 9 depleted)
✅ 18 figuras individuales PNG (9 enriched + 9 depleted)
✅  9 figuras combinadas PDF
✅  9 figuras combinadas SVG
✅  9 figuras combinadas PNG
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   81 archivos totales
```

---

## 🆕 Mejoras Implementadas (3 iteraciones)

### Iteración 1: Separación y Publication-Ready

✅ **Separación enriched/depleted**: Dos bubble plots independientes por análisis
✅ **Lint-harden**: Validación de inputs, logging estructurado, reproducibilidad (seed)
✅ **Publication-ready presets**: Dimensiones explícitas (180x120 mm, double_col)
✅ **Paleta colorblind-safe**: Okabe-Ito (verde #009E73, naranja #D55E00)
✅ **Exports determinísticos**: PDF + PNG con dimensiones y DPI explícitos

### Iteración 2: SVG y Tipografía

✅ **svglite instalado**: Exports SVG de calidad con texto editable
✅ **Tamaño de letra aumentado**: `axis.text.y` usa `p$axis` (8.5pt) en lugar de `p$ticks` (7.5pt)
  - **Antes**: Términos en 7.5pt (poco legibles)
  - **Ahora**: Términos en 8.5pt (más claros a la vista)

### Iteración 3: Figuras Combinadas

✅ **patchwork integrado**: Combina enriched + depleted en una sola figura
✅ **Layout horizontal**: Dos paneles lado a lado (2.5x width)
✅ **Dimensiones**: 450 x 240 mm (doble ancho para acomodar ambos plots)
✅ **Título unificado**: "miRPathDB Enrichment Summary: <analysis_id>"
✅ **Leyendas colectadas**: Una sola leyenda compartida en el lado derecho

---

## 📁 Estructura de Archivos

### Nomenclatura

```
# Individuales
Bubble_miRPathDB_enriched_<analysis_id>_A_conservative.{pdf,svg,png}
Bubble_miRPathDB_depleted_<analysis_id>_A_conservative.{pdf,svg,png}

# Combinadas
Bubble_miRPathDB_combined_<analysis_id>_A_conservative.{pdf,svg,png}
```

### Ejemplo: GENERO_1_vs_0

```
results/figures/MiEAA_GSEA_bubble/A_conservative/
├── Bubble_miRPathDB_enriched_GENERO_1_vs_0_A_conservative.pdf   (28 KB)
├── Bubble_miRPathDB_enriched_GENERO_1_vs_0_A_conservative.svg   (28 KB)
├── Bubble_miRPathDB_enriched_GENERO_1_vs_0_A_conservative.png   (948 KB)
├── Bubble_miRPathDB_depleted_GENERO_1_vs_0_A_conservative.pdf   (28 KB)
├── Bubble_miRPathDB_depleted_GENERO_1_vs_0_A_conservative.svg   (28 KB)
├── Bubble_miRPathDB_depleted_GENERO_1_vs_0_A_conservative.png   (926 KB)
├── Bubble_miRPathDB_combined_GENERO_1_vs_0_A_conservative.pdf   (33 KB)
├── Bubble_miRPathDB_combined_GENERO_1_vs_0_A_conservative.svg   (56 KB)
└── Bubble_miRPathDB_combined_GENERO_1_vs_0_A_conservative.png   (1.7 MB)
```

---

## 🎨 Características Técnicas

### Dimensiones y Styling

| Aspecto | Valor |
|---------|-------|
| **Preset** | `double_col` (journal double column) |
| **Individuales** | 225 x 240 mm (180 × 1.25 width, 120 × 2 height) |
| **Combinadas** | 450 x 240 mm (180 × 2.5 width, 120 × 2 height) |
| **Base font** | 8.5 pt Helvetica |
| **Títulos** | 10 pt bold |
| **Ejes** | 8.5 pt |
| **Términos (Y)** | 8.5 pt (↑ mejorado desde 7.5 pt) |
| **Ticks (X)** | 7.5 pt |
| **Leyenda** | 7.5 pt bold title, 7.5 pt text |

### Exports y Formatos

| Formato | DPI | Device | Tamaño típico | Uso |
|---------|-----|--------|---------------|-----|
| **PDF** | 300 | cairo_pdf | 28-33 KB | Journals, publicaciones |
| **SVG** | 300 | svglite | 28-56 KB | Edición vectorial |
| **PNG** | 600 | png | 800 KB - 1.7 MB | Presentaciones, web |

### Paleta de Colores

```
Enriched:  #009E73 (verde Okabe-Ito)
Depleted:  #D55E00 (naranja Okabe-Ito)
Other:     #999999 (gris neutro)
```

✅ Colorblind-safe (deuteranopia, protanopia)
✅ Contraste óptimo para impresión B/N
✅ Colores únicos por plot (no leyenda de color necesaria)

---

## 🚀 Uso del Script

### Ejecución Básica

```bash
cd ~/bioinfo/projects/mirna_glioma
conda activate omics-R
Rscript scripts/10_miEAA_GSEA_bubble_plots.R
```

### Parámetros Disponibles

```bash
Rscript scripts/10_miEAA_GSEA_bubble_plots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_bubble \
  --run_tag A_conservative \
  --n_mirpathdb 15 \          # Top N términos por DB (default: 12)
  --wrap_width 35 \            # Ancho de wrap para labels (default: 40)
  --preset double_col \        # single_col | double_col | presentation | poster
  --seed 42                    # Seed para reproducibilidad (default: 42)
```

### Ejemplos de Uso

```bash
# Mostrar top 15 términos en lugar de 12
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --n_mirpathdb 15

# Usar preset single column (más pequeño)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --preset single_col

# Preset presentation (para slides)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --preset presentation

# Labels más cortos (wrap a 30 caracteres)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --wrap_width 30
```

---

## 📋 QC Checklist (Figuras Verificadas)

### Pre-export ✅
- [x] Input validation: Spec existe, `in_root` existe, parámetros válidos
- [x] Data availability: Términos disponibles para cada dirección y DB
- [x] Direction filtering: Filtra correctamente `Enriched` y `Depleted`
- [x] Top N selection: Toma top 12 por DB ordenado por p-adjusted

### Styling ✅
- [x] Preset consistency: Todas las dimensiones tomadas de preset
- [x] Font sizing: Términos en 8.5pt (mejorado), títulos 10pt, ejes 8.5pt
- [x] Colorblind-safe: Paleta Okabe-Ito (verde/naranja)
- [x] Contrast check: Colores distinguibles en B/N y pantalla
- [x] Grid lines: Major X visible, minor removidos

### Layout ✅
- [x] Label overlap: `str_wrap(width = 40)` previene overlap
- [x] Tick density: Faceting con scales libres evita overcrowding
- [x] Margins: Explícitos (7mm base + 6mm derecha)
- [x] Legend clarity: Position right, key size 4mm, title bold
- [x] Panel spacing: 2.5mm entre facets

### Export ✅
- [x] Vector formats: PDF (cairo_pdf) + SVG (svglite) generados
- [x] Raster format: PNG a 600 dpi
- [x] Font embedding: cairo_pdf embebe fuentes
- [x] File verification: Script verifica PDF, lanza error si falla
- [x] Deterministic naming: Nombres estables, no timestamps

### Reproducibility ✅
- [x] Seed set: `set.seed(42)` al inicio
- [x] Stable ordering: Factor levels fijados para `db` y `direction`
- [x] Logged parameters: Todos en log header
- [x] Version tracking: Logs con timestamp, parámetros, dimensiones

### Multi-panel (Combinadas) ✅
- [x] Patchwork integration: Combina enriched + depleted
- [x] Horizontal layout: Lado a lado, 2.5x width
- [x] Shared legend: Colectada en lado derecho
- [x] Unified title: Título centrado superior
- [x] Proper spacing: Panel spacing 2.5mm

---

## 📦 Dependencias

### R packages (todos instalados en omics-R)

```r
library(ggplot2)      # Plotting
library(dplyr)        # Data manipulation
library(stringr)      # String wrapping
library(patchwork)    # Multi-panel figures
library(svglite)      # SVG exports (instalado en esta sesión)
```

### Verificar cairo_pdf

```r
capabilities("cairo")  # Debe retornar TRUE
```

---

## 🎯 Comparación Antes vs Ahora

| Aspecto | Antes | Ahora |
|---------|-------|-------|
| **Plots por análisis** | 1 (mezclado) | 3 (enriched, depleted, combinado) |
| **Formatos de export** | PDF + PNG | PDF + SVG + PNG |
| **SVG disponible** | ❌ No | ✅ Sí (svglite instalado) |
| **Tamaño letra términos** | 7.5pt (pequeño) | 8.5pt (más legible) |
| **Figura combinada** | ❌ No | ✅ Sí (patchwork) |
| **Ancho figura combinada** | N/A | 450mm (2.5x base) |
| **Colorblind-safe** | ✅ Sí | ✅ Sí (sin cambios) |
| **Reproducibilidad** | ✅ Sí | ✅ Sí (seed configurable) |
| **Total archivos** | 18 (9 PDFs + 9 PNGs) | 81 (27 PDFs + 27 SVGs + 27 PNGs) |

---

## 🔧 Opciones de Personalización

### 1. Cambiar Top N de Términos

```bash
# Top 15 en lugar de 12
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --n_mirpathdb 15

# Top 20 (más procesos, figura más larga)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --n_mirpathdb 20
```

### 2. Ajustar Ancho de Labels

```bash
# Labels más cortos (wrap a 30 caracteres)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --wrap_width 30

# Labels más largos (wrap a 50 caracteres)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --wrap_width 50
```

### 3. Usar Otros Presets

```bash
# Single column (85 x 65 mm base)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --preset single_col

# Presentation (254 x 143 mm, 16:9)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --preset presentation

# Poster (508 x 356 mm, muy grande)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --preset poster
```

### 4. Cambiar Seed de Reproducibilidad

```bash
# Usar seed diferente (si hay jitter/randomness)
Rscript scripts/10_miEAA_GSEA_bubble_plots.R --seed 123
```

---

## 📊 Ejemplo de Salida: GENERO_1_vs_0

### Términos Encontrados

```
Enriched terms:
  GO Biological process (miRPathDB) = 12
  GO Molecular function (miRPathDB) = 12
  KEGG (miRPathDB) = 12
  Reactome (miRPathDB) = 12
  Total: 48 términos enriquecidos

Depleted terms:
  GO Biological process (miRPathDB) = 12
  GO Molecular function (miRPathDB) = 12
  KEGG (miRPathDB) = 12
  Reactome (miRPathDB) = 12
  Total: 48 términos depletados
```

### Archivos Generados

```
✅ Bubble_miRPathDB_enriched_GENERO_1_vs_0_A_conservative.pdf   (28 KB)
✅ Bubble_miRPathDB_enriched_GENERO_1_vs_0_A_conservative.svg   (28 KB)
✅ Bubble_miRPathDB_enriched_GENERO_1_vs_0_A_conservative.png   (948 KB)
✅ Bubble_miRPathDB_depleted_GENERO_1_vs_0_A_conservative.pdf   (28 KB)
✅ Bubble_miRPathDB_depleted_GENERO_1_vs_0_A_conservative.svg   (28 KB)
✅ Bubble_miRPathDB_depleted_GENERO_1_vs_0_A_conservative.png   (926 KB)
✅ Bubble_miRPathDB_combined_GENERO_1_vs_0_A_conservative.pdf   (33 KB)
✅ Bubble_miRPathDB_combined_GENERO_1_vs_0_A_conservative.svg   (56 KB)
✅ Bubble_miRPathDB_combined_GENERO_1_vs_0_A_conservative.png   (1.7 MB)
```

---

## 🎓 Skills Aplicados

Los siguientes Claude Code Skills fueron utilizados en este proyecto:

1. **/lint-harden-pro**: Robustez, validación de inputs, logging, reproducibilidad
2. **/pub-figures**: Presets publication-ready, exports determinísticos, QC
3. **/color-palettes**: Paleta Okabe-Ito colorblind-safe
4. **/multipanel-figures**: Figuras combinadas con patchwork
5. **/specialized-omics-plots**: Bubble plots optimizados para omics

---

## 📚 Documentación Completa

- **Script modificado**: `scripts/10_miEAA_GSEA_bubble_plots.R`
- **Cambios detallados**: `docs/10_bubble_plots_CHANGES.md`
- **Este resumen**: `docs/10_bubble_plots_FINAL_SUMMARY.md`
- **Logs de ejecución**: `logs/mieaa_gsea_bubble_*.txt`

---

## ✅ Estado Final

**Todo completado exitosamente**:

- ✅ Separación de gráficas enriched vs depleted
- ✅ svglite instalado (exports SVG de calidad)
- ✅ Tamaño de letra de términos aumentado (8.5pt)
- ✅ Figuras combinadas implementadas con patchwork
- ✅ 81 archivos generados (3 formatos × 3 tipos × 9 análisis)
- ✅ QC completo pasado
- ✅ Script probado y funcionando

**Listo para publicación** 🎉
