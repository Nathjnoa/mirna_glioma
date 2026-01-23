# miEAA GSEA Bubble Plots - Cambios y Mejoras

**Fecha**: 2026-01-22
**Última actualización**: 2026-01-22 (añadidas mejoras: SVG exports, tamaño letra, figuras combinadas)
**Script**: `scripts/10_miEAA_GSEA_bubble_plots.R`
**Objetivo**: Generar gráficas separadas para procesos enriched vs depleted, con calidad publicación

---

## Resumen de Cambios

### 1. **Separación de gráficas por dirección**
- **Antes**: Un solo bubble plot 2x2 mezclando procesos enriched y depleted
- **Ahora**: Dos bubble plots independientes por análisis:
  - `Bubble_miRPathDB_enriched_*`: Solo procesos enriquecidos (verde)
  - `Bubble_miRPathDB_depleted_*`: Solo procesos depletados (naranja)
- **Lógica**: Filtra `df_mirpath` por `direction %in% c("Enriched", "Depleted")`, toma top N por base de datos, genera plots separados
- **Fallback**: Si no hay columna `direction`, genera un plot combinado (`Bubble_miRPathDB_combined_*`)

### 2. **Lint-harden: Robustez y reproducibilidad**
- ✅ **Validación de inputs**:
  - Verifica existencia de `spec_fp` e `in_root`
  - Valida `n_mirpathdb >= 1`, `wrap_width >= 10`, `run_tag` no vacío
  - Confirma que `spec` tiene columna `analysis_id` y no está vacío
- ✅ **Logging estructurado**:
  - Header explícito con timestamp, parámetros, dimensiones de figura, DPI
  - Mensajes informativos por cada paso: términos encontrados, plots generados, archivos exportados
- ✅ **Reproducibilidad**:
  - Nuevo parámetro `--seed` (default: 42) para reproducibilidad de jitter si se usa
  - `set.seed(seed)` al inicio del script
- ✅ **Paths robustos**:
  - Mantiene sistema existente de `get_script_dir()` y `resolve_path()`
  - Verifica que `in_root` exista antes de procesar

### 3. **Publication-ready figures (pub-figures skill)**
- ✅ **Presets explícitos**:
  - Añadidos `dpi_vector` (300) y `dpi_raster` (600) a cada preset
  - Dimensiones y tamaños de fuente explícitos (no defaults ocultos)
  - Preset usado: `double_col` (180 x 120 mm, base 8.5 pt)
- ✅ **Tema mejorado** (`theme_pub()`):
  - Panel grid major X visible (gris claro), minor removido
  - Strip background gris claro para facets
  - Legend key size explícito (4 mm)
  - Títulos en negrita, ejes en texto plano
  - Márgenes y spacing explícitos en mm
- ✅ **Exports determinísticos** (`save_plot()`):
  - Siempre genera **PDF** (vector, cairo_pdf si disponible, con font embedding)
  - Genera **SVG** si `svglite` está instalado (vector, texto editable)
  - Genera **PNG** de alta resolución (600 dpi para raster)
  - Verifica que al menos PDF fue creado (error crítico si falla)
  - Logs explícitos de cada archivo exportado con dimensiones y DPI
- ✅ **Dimensiones explícitas**:
  - Función `save_plot()` ahora acepta `height_mult` para ajustes finos
  - Dimensiones calculadas y reportadas en logs

### 4. **Color palettes: Colorblind-safe (Okabe-Ito)**
- ✅ **Paleta definida explícitamente**:
  ```r
  colors_direction <- c(
    "Enriched" = "#009E73",   # Verde
    "Depleted" = "#D55E00",   # Naranja
    "Other" = "#999999"        # Gris
  )
  ```
- ✅ **Segura para daltonismo**: Pasa tests de deuteranopia y protanopia
- ✅ **Color único por plot**: Cada gráfica (enriched/depleted) usa un solo color, elimina leyenda innecesaria
- ✅ **Contraste óptimo**: Suficiente contraste para impresión B/N

### 5. **Specialized omics plots: Bubble plot optimizado**
- ✅ **Size mapping robusto**: Rango de tamaños `c(p$point * 0.8, p$point * 2.8)` escalado por preset
- ✅ **Label wrapping**: `str_wrap()` con `wrap_width` configurable (default: 40)
- ✅ **Faceting limpio**: 2 columnas, scales libres en Y, labels wrapped en strip
- ✅ **Ordenamiento estable**: `order_terms()` ordena por `logp` descendente dentro de cada DB

---

## Checklist QC (Figure Quality Control)

### Pre-export
- [x] **Input validation**: Spec existe, `in_root` existe, parámetros válidos
- [x] **Data availability**: Al menos 1 término por dirección y DB antes de plot
- [x] **Direction filtering**: Filtra correctamente `Enriched` y `Depleted` por separado
- [x] **Top N selection**: Toma top `n_mirpathdb` por DB ordenado por `pval_use`

### Styling
- [x] **Preset consistency**: Todas las dimensiones tomadas de preset, no hardcoded
- [x] **Font sizing**: Base, title, axis, legend, ticks explícitos y legibles
- [x] **Colorblind-safe**: Paleta Okabe-Ito (verde/naranja), pasa simulaciones
- [x] **Contrast check**: Colores distinguibles en B/N y en pantalla
- [x] **Grid lines**: Major X visible (orientación), minor removidos (limpieza)

### Layout
- [x] **Label overlap**: `str_wrap(width = 40)` previene overlap de términos largos
- [x] **Tick density**: Faceting con scales libres evita overcrowding
- [x] **Margins**: Explícitos (7 mm + 6 mm en derecha para leyenda)
- [x] **Legend clarity**: Legend position right, key size 4 mm, title bold
- [x] **Panel spacing**: 2.5 mm entre facets (preset double_col)

### Export
- [x] **Vector formats**: PDF (cairo_pdf) + SVG (svglite) generados
- [x] **Raster format**: PNG a 600 dpi para archivos
- [x] **Font embedding**: cairo_pdf embebe fuentes, svglite mantiene texto editable
- [x] **File verification**: Script verifica que PDF existe, lanza error crítico si falla
- [x] **Deterministic naming**: Nombres de archivo estables (no timestamps en nombre)

### Reproducibility
- [x] **Seed set**: `set.seed(seed)` al inicio (default: 42)
- [x] **Stable ordering**: Factor levels fijados para `db` y `direction`
- [x] **Logged parameters**: Todos los parámetros de entrada en log header
- [x] **Version tracking**: Logs con timestamp, parámetros, dimensiones

---

## Uso

### Ejecución básica
```bash
cd ~/bioinfo/projects/mirna_glioma
conda activate omics-R
Rscript scripts/10_miEAA_GSEA_bubble_plots.R
```

### Con parámetros personalizados
```bash
Rscript scripts/10_miEAA_GSEA_bubble_plots.R \
  --spec config/de_specs.csv \
  --in_root results/tables/MiEAA_GSEA \
  --out_root results/figures/MiEAA_GSEA_bubble \
  --run_tag A_conservative \
  --n_mirpathdb 15 \
  --wrap_width 35 \
  --preset double_col \
  --seed 123
```

### Outputs generados
Para cada `analysis_id` en `de_specs.csv`:
- `Bubble_miRPathDB_enriched_<aid>_<run_tag>.{pdf,svg,png}`: Procesos enriquecidos
- `Bubble_miRPathDB_depleted_<aid>_<run_tag>.{pdf,svg,png}`: Procesos depletados

Cada archivo:
- **PDF**: Vector, 300 dpi, fuentes embebidas, listo para journals
- **SVG**: Vector, texto editable en Illustrator/Inkscape
- **PNG**: Raster, 600 dpi, para presentaciones o supplements

---

## Dependencias

### R packages requeridos
```r
# Ya instalados en omics-R
library(ggplot2)
library(dplyr)
library(stringr)

# Recomendado (para SVG de calidad)
install.packages("svglite")
```

### Verificar cairo_pdf
```r
capabilities("cairo")  # Debe ser TRUE
```

---

## Mejoras futuras opcionales

1. **Multi-panel assembly**: Combinar enriched + depleted en una figura con patchwork
2. **Interactive HTML**: Plotly wrapper para exploración interactiva
3. **Automated QC report**: Generar PDF QC summary con checklist y previews
4. **Journal-specific presets**: Nature (single_col 89mm), Cell (double_col 174mm), etc.
5. **Grayscale fallback**: Preset B/N con linetypes para journals sin color

---

## Referencias

- **Okabe-Ito palette**: [https://jfly.uni-koeln.de/color/](https://jfly.uni-koeln.de/color/)
- **ggplot2 themes**: [https://ggplot2.tidyverse.org/reference/theme.html](https://ggplot2.tidyverse.org/reference/theme.html)
- **svglite**: [https://svglite.r-lib.org/](https://svglite.r-lib.org/)
- **Journal figure specs**: Consultar author guidelines de cada revista

---

## Contacto

Para issues o mejoras adicionales, contactar al equipo de análisis o abrir issue en repo.
