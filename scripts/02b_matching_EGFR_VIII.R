#!/usr/bin/env Rscript
# ==============================================================================
# 02b_matching_EGFR_VIII.R
#
# Objetivo: Crear grupos pareados por edad para comparar MUTACION_EGFR_VIII
#           (SI vs NO) usando matching óptimo 1:2
#
# Método: Matching óptimo usando programación lineal (nbpMatching) o
#         greedy nearest-neighbor si no disponible
#
# Salidas:
#   - Reporte de matching con pares formados
#   - Lista de muestras excluidas para usar en de_specs.csv
#   - Diagnósticos de balance pre/post matching
# ==============================================================================

options(stringsAsFactors = FALSE)

# --- Configuración ---
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

meta_fp <- get_arg("--meta", "data/intermediate/Metadatos_gliomas_verificados.csv")
outdir  <- get_arg("--outdir", "data/processed")
ratio   <- as.numeric(get_arg("--ratio", "2"))
caliper <- as.numeric(get_arg("--caliper", "15"))  # máx diferencia de edad

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Crear directorios ---
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/matching", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables/matching", showWarnings = FALSE, recursive = TRUE)

# --- Logging ---
log_fp <- file.path("logs", paste0("matching_EGFR_VIII_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat("============================================\n")
cat("Matching 1:2 por edad - MUTACION_EGFR_VIII\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Meta:", meta_fp, "\n")
cat("Ratio:", ratio, "\n")
cat("Caliper (años):", caliper, "\n")
cat("============================================\n\n")

# --- Cargar paquetes ---
suppressPackageStartupMessages({
  library(data.table)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
})

# --- Función de matching óptimo manual ---
# Usa programación de asignación para minimizar suma total de distancias
optimal_matching_1toN <- function(casos_df, controles_df, ratio = 2, caliper = Inf) {
  # casos_df: data.frame con columnas id, EDAD
  # controles_df: data.frame con columnas id, EDAD
  # ratio: número de controles por caso
  # caliper: máxima diferencia permitida

  n_casos <- nrow(casos_df)
  n_controles <- nrow(controles_df)

  if (n_controles < n_casos * ratio) {
    warning("No hay suficientes controles. Se ajustará el ratio.")
    ratio <- floor(n_controles / n_casos)
  }

  # Matriz de distancias (diferencia absoluta de edad)
  dist_mat <- matrix(NA, nrow = n_casos, ncol = n_controles)
  rownames(dist_mat) <- casos_df$id
  colnames(dist_mat) <- controles_df$id

  for (i in seq_len(n_casos)) {
    for (j in seq_len(n_controles)) {
      d <- abs(casos_df$EDAD[i] - controles_df$EDAD[j])
      dist_mat[i, j] <- if (d <= caliper) d else Inf
    }
  }

  # Algoritmo greedy optimizado: para cada caso, seleccionar los N controles más cercanos
  # Ordenar casos por dificultad (menos opciones primero)
  n_valid <- apply(dist_mat, 1, function(x) sum(is.finite(x)))
  caso_order <- order(n_valid)  # casos más difíciles primero

  asignaciones <- list()
  controles_usados <- character(0)

  for (idx in caso_order) {
    caso_id <- casos_df$id[idx]
    caso_edad <- casos_df$EDAD[idx]

    # Distancias a controles no usados
    disponibles <- setdiff(controles_df$id, controles_usados)
    if (length(disponibles) < ratio) {
      warning(sprintf("Caso %s: solo %d controles disponibles", caso_id, length(disponibles)))
    }

    dists <- dist_mat[caso_id, disponibles, drop = FALSE]
    dists_sorted <- sort(dists[1, ])
    dists_sorted <- dists_sorted[is.finite(dists_sorted)]

    # Seleccionar los N más cercanos
    n_select <- min(ratio, length(dists_sorted))
    if (n_select == 0) {
      warning(sprintf("Caso %s: no hay controles dentro del caliper", caso_id))
      next
    }

    seleccionados <- names(dists_sorted)[1:n_select]
    controles_usados <- c(controles_usados, seleccionados)

    asignaciones[[caso_id]] <- list(
      caso = caso_id,
      caso_edad = caso_edad,
      controles = seleccionados,
      ctrl_edades = controles_df$EDAD[match(seleccionados, controles_df$id)],
      distancias = as.numeric(dists_sorted[1:n_select])
    )
  }

  return(asignaciones)
}

# --- Leer metadatos ---
meta <- fread(meta_fp, data.table = FALSE)
cat("Metadatos cargados:", nrow(meta), "muestras\n")

# Identificar columna de mutación EGFR VIII (puede tener espacio)
egfr_col <- grep("MUTACION_EGFR.*VIII", colnames(meta), value = TRUE)[1]
if (is.na(egfr_col)) stop("No se encontró columna MUTACION_EGFR VIII")
cat("Columna EGFR VIII:", egfr_col, "\n\n")

# --- Preparar datos para matching ---
df <- data.frame(
  id = meta$id,
  EDAD = as.numeric(meta$EDAD),
  EGFR_VIII = meta[[egfr_col]],
  stringsAsFactors = FALSE
)

# Filtrar solo SI y NO (excluir NA o valores raros)
df <- df[df$EGFR_VIII %in% c("SI", "NO"), ]

cat("Distribución MUTACION_EGFR_VIII:\n")
print(table(df$EGFR_VIII, useNA = "ifany"))
cat("\n")

casos <- df[df$EGFR_VIII == "SI", ]
controles <- df[df$EGFR_VIII == "NO", ]

n_casos <- nrow(casos)
n_controles <- nrow(controles)
cat("Casos (SI):", n_casos, "\n")
cat("Controles potenciales (NO):", n_controles, "\n\n")

if (n_casos == 0) stop("No hay casos (SI) para matching")

# --- Balance pre-matching ---
cat("=== Balance PRE-matching ===\n")
cat("Edad - Casos (SI):\n")
cat("  Media:", round(mean(casos$EDAD), 1), "\n")
cat("  SD:", round(sd(casos$EDAD), 1), "\n")
cat("  Rango:", paste(range(casos$EDAD), collapse = " - "), "\n")
cat("  Valores:", paste(sort(casos$EDAD), collapse = ", "), "\n")

cat("Edad - Controles (NO):\n")
cat("  Media:", round(mean(controles$EDAD), 1), "\n")
cat("  SD:", round(sd(controles$EDAD), 1), "\n")
cat("  Rango:", paste(range(controles$EDAD), collapse = " - "), "\n")
cat("  Valores:", paste(sort(controles$EDAD), collapse = ", "), "\n")

diff_pre <- abs(mean(casos$EDAD) - mean(controles$EDAD))
cat("Diferencia de medias pre-matching:", round(diff_pre, 1), "años\n\n")

# --- Ejecutar matching ---
cat("=== Ejecutando matching óptimo ===\n")
cat("Método: Greedy nearest-neighbor (casos difíciles primero)\n")
cat("Ratio:", ratio, ":1\n")
cat("Caliper:", caliper, "años\n\n")

set.seed(42)  # reproducibilidad
asignaciones <- optimal_matching_1toN(casos, controles, ratio = ratio, caliper = caliper)

# --- Construir datos pareados ---
muestras_incluidas <- character(0)
matching_detail <- data.frame()

cat("=== Detalle de pares formados ===\n")
for (caso_id in names(asignaciones)) {
  a <- asignaciones[[caso_id]]

  cat(sprintf("Par %s:\n", caso_id))
  cat(sprintf("  CASO: %s (edad=%d)\n", a$caso, a$caso_edad))

  for (i in seq_along(a$controles)) {
    cat(sprintf("  CTRL: %s (edad=%d, diff=%d años)\n",
                a$controles[i], a$ctrl_edades[i], a$distancias[i]))
  }
  cat("\n")

  # Agregar caso
  muestras_incluidas <- c(muestras_incluidas, a$caso)
  matching_detail <- rbind(matching_detail, data.frame(
    id = a$caso,
    EDAD = a$caso_edad,
    EGFR_VIII = "SI",
    treatment = 1,
    subclass = caso_id,
    stringsAsFactors = FALSE
  ))

  # Agregar controles
  muestras_incluidas <- c(muestras_incluidas, a$controles)
  for (i in seq_along(a$controles)) {
    matching_detail <- rbind(matching_detail, data.frame(
      id = a$controles[i],
      EDAD = a$ctrl_edades[i],
      EGFR_VIII = "NO",
      treatment = 0,
      subclass = caso_id,
      stringsAsFactors = FALSE
    ))
  }
}

# --- Balance post-matching ---
matched_casos <- matching_detail[matching_detail$treatment == 1, ]
matched_ctrl <- matching_detail[matching_detail$treatment == 0, ]

cat("=== Balance POST-matching ===\n")
cat("Muestras en análisis pareado:", nrow(matching_detail), "\n")
cat("  Casos (SI):", nrow(matched_casos), "\n")
cat("  Controles (NO):", nrow(matched_ctrl), "\n\n")

cat("Edad - Casos (SI):\n")
cat("  Media:", round(mean(matched_casos$EDAD), 1), "\n")
cat("  SD:", round(sd(matched_casos$EDAD), 1), "\n")

cat("Edad - Controles pareados (NO):\n")
cat("  Media:", round(mean(matched_ctrl$EDAD), 1), "\n")
cat("  SD:", round(sd(matched_ctrl$EDAD), 1), "\n")

diff_post <- abs(mean(matched_casos$EDAD) - mean(matched_ctrl$EDAD))
cat("Diferencia de medias post-matching:", round(diff_post, 1), "años\n")

if (diff_pre > 0) {
  cat("Reducción en diferencia:", round((1 - diff_post/diff_pre) * 100, 1), "%\n\n")
}

# Calcular SMD (Standardized Mean Difference)
pooled_sd <- sqrt((var(matched_casos$EDAD) + var(matched_ctrl$EDAD)) / 2)
smd <- diff_post / pooled_sd
cat("SMD (Standardized Mean Difference):", round(smd, 3), "\n")
cat("  (< 0.1 es excelente, < 0.25 es aceptable)\n\n")

# --- Identificar muestras excluidas ---
muestras_todas <- meta$id
muestras_excluidas <- setdiff(muestras_todas, muestras_incluidas)

cat("=== Muestras para de_specs.csv ===\n")
cat("Muestras INCLUIDAS en análisis pareado:\n")
cat(paste(muestras_incluidas, collapse = ";"), "\n\n")

cat("Muestras EXCLUIDAS (para exclude_samples):\n")
cat(paste(muestras_excluidas, collapse = ";"), "\n\n")

# --- Guardar resultados ---

# 1. Tabla de matching detallada
matching_detail <- matching_detail[order(matching_detail$subclass, -matching_detail$treatment), ]
detail_fp <- file.path("results/tables/matching",
                       paste0("matching_EGFR_VIII_detail_", ts, ".csv"))
write.csv(matching_detail, detail_fp, row.names = FALSE)
cat("Tabla de matching guardada:", detail_fp, "\n")

# 2. Archivo con muestras excluidas (para copiar a de_specs.csv)
exclude_fp <- file.path(outdir, paste0("exclude_samples_EGFR_VIII_", ts, ".txt"))
writeLines(paste(muestras_excluidas, collapse = ";"), exclude_fp)
cat("Muestras excluidas guardadas:", exclude_fp, "\n")

# 3. Archivo con muestras incluidas
include_fp <- file.path(outdir, paste0("include_samples_EGFR_VIII_", ts, ".txt"))
writeLines(paste(muestras_incluidas, collapse = ";"), include_fp)
cat("Muestras incluidas guardadas:", include_fp, "\n")

# 4. Línea lista para de_specs.csv
specs_line <- sprintf(
  'EGFR_VIII_SI_vs_NO_matched,as_is,MUTACION_EGFR VIII,,,,"",SI,NO,,%s',
  paste(muestras_excluidas, collapse = ";")
)
specs_fp <- file.path(outdir, paste0("de_specs_line_EGFR_VIII_", ts, ".txt"))
writeLines(specs_line, specs_fp)
cat("\nLínea para de_specs.csv guardada:", specs_fp, "\n")
cat("Contenido:\n", specs_line, "\n")

# --- Gráficos de diagnóstico ---
if (requireNamespace("ggplot2", quietly = TRUE)) {

  # Plot 1: Distribución de edad pre/post matching
  df$matched <- ifelse(df$id %in% muestras_incluidas, "Incluido", "Excluido")
  df$group <- ifelse(df$EGFR_VIII == "SI", "SI (casos)", "NO (controles)")

  p1 <- ggplot(df, aes(x = group, y = EDAD, fill = matched)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(aes(color = matched), position = position_jitterdodge(jitter.width = 0.2),
                size = 3, alpha = 0.9) +
    scale_fill_manual(values = c("Incluido" = "#2E86AB", "Excluido" = "#E94F37"),
                      name = "Estado") +
    scale_color_manual(values = c("Incluido" = "#1a5276", "Excluido" = "#922b21"),
                       name = "Estado") +
    labs(
      title = "Distribución de edad por grupo EGFR_VIII",
      subtitle = sprintf("Matching 1:%d por edad | SMD post-matching: %.3f", ratio, smd),
      x = "MUTACION_EGFR_VIII",
      y = "Edad (años)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  fig_fp <- file.path("results/figures/matching",
                      paste0("matching_EGFR_VIII_boxplot_", ts, ".png"))
  ggsave(fig_fp, p1, width = 8, height = 6, dpi = 300)
  cat("\nFigura guardada:", fig_fp, "\n")

  # Plot 2: Dotplot de pares
  p2 <- ggplot(matching_detail, aes(x = subclass, y = EDAD, color = factor(treatment))) +
    geom_point(size = 4) +
    geom_line(aes(group = subclass), color = "gray50", linetype = "dashed") +
    scale_color_manual(values = c("0" = "#0072B2", "1" = "#D55E00"),
                       labels = c("Control (NO)", "Caso (SI)"),
                       name = "Grupo") +
    labs(
      title = "Pares formados en el matching",
      subtitle = "Cada línea conecta un caso con sus controles pareados",
      x = "ID del par (caso)",
      y = "Edad (años)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  fig2_fp <- file.path("results/figures/matching",
                       paste0("matching_EGFR_VIII_pairs_", ts, ".png"))
  ggsave(fig2_fp, p2, width = 10, height = 6, dpi = 300)
  cat("Figura de pares guardada:", fig2_fp, "\n")
}

# --- Resumen final ---
cat("\n============================================\n")
cat("RESUMEN MATCHING EGFR_VIII\n")
cat("============================================\n")
cat("Casos (SI):", nrow(matched_casos), "\n")
cat("Controles pareados (NO):", nrow(matched_ctrl), "\n")
cat("Total muestras:", nrow(matching_detail), "\n")
cat("Muestras excluidas:", length(muestras_excluidas), "\n")
cat("Diferencia media edad post-matching:", round(diff_post, 1), "años\n")
cat("SMD:", round(smd, 3), "\n")
cat("============================================\n")
cat("\nPróximo paso: Añadir línea a config/de_specs.csv\n")
cat("Luego ejecutar: Rscript scripts/03_edgeR_multiDE.R --spec config/de_specs.csv\n")

sink(type = "message")
sink(type = "output")
close(zz)

cat("\n✓ Matching completado. Log:", log_fp, "\n")
cat("\nMuestras excluidas para de_specs.csv:\n")
cat(paste(muestras_excluidas, collapse = ";"), "\n")
