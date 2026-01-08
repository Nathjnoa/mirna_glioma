#!/usr/bin/env Rscript

# ============================================================
# 02_logcpm_mds.R
# - logCPM (edgeR + TMM)
# - MDS (sobre logCPM)
# - RLE pre/post normalización (sin/con TMM)
# - Excluye muestras ^A (post mortem) y columnas de anotación
# - Log -> logs/
# - Figuras -> results/figures/qc_mds/
# ============================================================

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

in_fp <- if (length(args) >= 1) args[1] else file.path("data", "intermediate", "Gliomas_all_counts_merged.csv")

# Puedes cambiar esto si algún día quieres NO filtrar A*
drop_prefix <- "^A"   # regex para eliminar muestras post mortem

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path("results", "figures", "qc_mds")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
tab_dir <- file.path("results", "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_fp <- file.path("logs", paste0("qc_logcpm_mds_", ts, ".txt"))

zz <- file(log_fp, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat("============================================\n")
cat("QC: logCPM + MDS (edgeR)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Input:", in_fp, "\n")
cat("Drop prefix regex:", drop_prefix, "\n")
cat("Figuras:", fig_dir, "\n")
cat("Tablas:", tab_dir, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

# ---- Dependencias ----
if (!requireNamespace("edgeR", quietly = TRUE)) stop("edgeR no disponible.")
if (!requireNamespace("limma", quietly = TRUE)) stop("limma no disponible (necesario para MDS).")

cat("edgeR version:", as.character(utils::packageVersion("edgeR")), "\n")
cat("limma version:", as.character(utils::packageVersion("limma")), "\n\n")

# ---- Lectura ----
if (!file.exists(in_fp)) stop("Archivo de entrada no encontrado: ", in_fp)

read_counts_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE)
  }
}

df <- read_counts_table(in_fp)
cat("Dimensiones (df):", nrow(df), "filas x", ncol(df), "columnas\n\n")

# ---- Columnas de anotación ----
anno_candidates <- c("id", "type", "entrez_id", "HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
cat("Columnas de anotación detectadas:\n")
print(anno_present); cat("\n")

# ---- Columnas de conteo ----
count_cols <- setdiff(colnames(df), anno_present)
if (length(count_cols) == 0) stop("No quedaron columnas de conteo tras excluir anotación.")

# ---- Filtrar muestras post mortem (prefijo A) ----
before_n <- length(count_cols)
drop_samples <- grep(drop_prefix, count_cols, value = TRUE)

cat("Filtrado de muestras (drop '", drop_prefix, "'):\n", sep = "")
cat("- Muestras antes:", before_n, "\n")
cat("- A eliminar:", length(drop_samples), "\n")
if (length(drop_samples) > 0) {
  cat("- Nombres eliminados:\n")
  print(drop_samples)
}
count_cols <- setdiff(count_cols, drop_samples)
cat("- Muestras después:", length(count_cols), "\n\n")

if (length(count_cols) == 0) stop("Tras filtrar no quedan muestras para QC.")

# ---- Matriz de conteos ----
feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))
feature_type_df <- data.frame(
  id = make.unique(as.character(feature_id)),
  type = if ("type" %in% colnames(df)) df$type else NA,
  stringsAsFactors = FALSE
)
counts_df <- df[, count_cols, drop = FALSE]

# coerción numérica segura
to_numeric <- function(x) suppressWarnings(as.numeric(x))
counts_df[] <- lapply(counts_df, to_numeric)

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

cat("Dimensiones (counts_mat):", nrow(counts_mat), "features x", ncol(counts_mat), "muestras\n\n")
cat("Chequeos básicos:\n")
cat("- NA totales:", sum(is.na(counts_mat)), "\n")
cat("- Negativos:", sum(counts_mat < 0, na.rm = TRUE), "\n")
cat("- No enteros (tol 1e-6):", sum(abs(counts_mat - round(counts_mat)) > 1e-6, na.rm = TRUE), "\n\n")

# ---- Filtrar muestras por library size (pre normalización) ----
min_lib_size <- 100000
lib_size_pre_raw <- colSums(counts_mat, na.rm = TRUE)
before_lib <- length(lib_size_pre_raw)
drop_lib <- names(lib_size_pre_raw)[lib_size_pre_raw < min_lib_size]

cat("Filtrado de muestras por library size (<", min_lib_size, "):\n", sep = "")
cat("- Muestras antes:", before_lib, "\n")
cat("- A eliminar:", length(drop_lib), "\n")
if (length(drop_lib) > 0) {
  cat("- Nombres eliminados:\n")
  print(drop_lib)
}

keep_lib <- lib_size_pre_raw >= min_lib_size
counts_mat <- counts_mat[, keep_lib, drop = FALSE]
lib_size_pre_raw <- lib_size_pre_raw[keep_lib]

cat("- Muestras después:", ncol(counts_mat), "\n\n")

if (ncol(counts_mat) == 0) {
  stop("No quedan muestras tras filtrar por library size (<100000).")
}

# ---- Limpieza mínima de features (antes de edgeR/logCPM/MDS) ----
min_expr <- 10
n_feat_before <- nrow(counts_mat)

cat("Limpieza mínima de features:\n")
cat("- Features iniciales:", n_feat_before, "\n")

# 1) remover all-zero
before_all_zero <- nrow(counts_mat)
keep_nz <- rowSums(counts_mat, na.rm = TRUE) > 0
removed_all_zero <- sum(!keep_nz)
counts_mat <- counts_mat[keep_nz, , drop = FALSE]
cat("  * All-zero: antes:", before_all_zero, "eliminadas:", removed_all_zero, "después:", nrow(counts_mat), "\n")

# 2) remover varianza 0 en conteos (robusto a NA)
before_var0 <- nrow(counts_mat)
row_var <- apply(counts_mat, 1, function(x) stats::var(x, na.rm = TRUE))
keep_var <- is.finite(row_var) & (row_var > 0)
removed_var0 <- sum(!keep_var)
counts_mat <- counts_mat[keep_var, , drop = FALSE]
cat("  * Varianza 0: antes:", before_var0, "eliminadas:", removed_var0, "después:", nrow(counts_mat), "\n")

# 3) filtrar baja expresión (suma de conteos < min_expr)
before_lowexpr <- nrow(counts_mat)
keep_lowexpr <- rowSums(counts_mat, na.rm = TRUE) >= min_expr
removed_lowexpr <- sum(!keep_lowexpr)
counts_mat <- counts_mat[keep_lowexpr, , drop = FALSE]
cat("  * Expresión <", min_expr, ": antes:", before_lowexpr, "eliminadas:", removed_lowexpr, "después:", nrow(counts_mat), "\n")

cat("- Features finales:", nrow(counts_mat), "\n\n")

if (nrow(counts_mat) < 1000) {
  stop("Tras limpieza mínima quedan <1000 features; revisar input o criterios.")
}

# ---- edgeR: DGEList + TMM ----
y <- edgeR::DGEList(counts = counts_mat)

# prior.count estabiliza log para conteos bajos/cero
prior_count <- 2

# ---- logCPM PRE normalización (sin TMM) ----
y_pre <- y
y_pre$samples$norm.factors <- rep(1, ncol(y$counts))

logcpm_pre <- edgeR::cpm(y_pre, log = TRUE, prior.count = prior_count)

# RLE pre: centrar por feature (mediana)
if (requireNamespace("matrixStats", quietly = TRUE)) {
  med_pre <- matrixStats::rowMedians(logcpm_pre, na.rm = TRUE)
} else {
  med_pre <- apply(logcpm_pre, 1, stats::median, na.rm = TRUE)
}
rle_pre <- sweep(logcpm_pre, 1, med_pre, FUN = "-")

png(file.path(fig_dir, paste0("01_RLE_preTMM_", ts, ".png")), width = 2200, height = 900, res = 150)
boxplot(rle_pre, outline = FALSE,
        main = "RLE pre-normalización (logCPM sin TMM)",
        ylab = "logCPM - mediana(feature)",
        las = 2, cex.axis = 0.6)
abline(h = 0)
dev.off()

# ---- logCPM POST normalización (TMM) ----
y <- edgeR::calcNormFactors(y, method = "TMM")

lib_size <- y$samples$lib.size
norm_factors <- y$samples$norm.factors
eff_lib  <- lib_size * norm_factors

cat("Resumen library size (raw):\n"); print(summary(lib_size)); cat("\n")
cat("Resumen norm.factors (TMM):\n"); print(summary(norm_factors)); cat("\n")
cat("Resumen effective library size (TMM):\n"); print(summary(eff_lib)); cat("\n")

logcpm_post <- edgeR::cpm(y, log = TRUE, prior.count = prior_count)

if (requireNamespace("matrixStats", quietly = TRUE)) {
  med_post <- matrixStats::rowMedians(logcpm_post, na.rm = TRUE)
} else {
  med_post <- apply(logcpm_post, 1, stats::median, na.rm = TRUE)
}
rle_post <- sweep(logcpm_post, 1, med_post, FUN = "-")

png(file.path(fig_dir, paste0("02_RLE_postTMM_", ts, ".png")), width = 2200, height = 900, res = 150)
boxplot(rle_post, outline = FALSE,
        main = "RLE post-normalización (logCPM con TMM)",
        ylab = "logCPM - mediana(feature)",
        las = 2, cex.axis = 0.6)
abline(h = 0)
dev.off()

cat("RLE pre y post generados (logCPM pre/post TMM).\n")
cat("Dimensiones (logcpm_post):", nrow(logcpm_post), "features x", ncol(logcpm_post), "muestras\n\n")

# ---- Exportar matriz normalizada (logCPM post TMM) ----
logcpm_fp <- file.path(tab_dir, paste0("logCPM_TMM_", ts, ".csv"))
utils::write.csv(
  cbind(feature_id = rownames(logcpm_post), as.data.frame(logcpm_post, check.names = FALSE)),
  logcpm_fp,
  row.names = FALSE
)
cat("Matriz logCPM (post TMM) exportada en:", logcpm_fp, "\n\n")

# ---- Density plots (logCPM) pre vs post TMM ----
xlim_all <- range(c(logcpm_pre, logcpm_post), finite = TRUE)

png(file.path(fig_dir, paste0("03_density_logCPM_preTMM_", ts, ".png")),
    width = 1600, height = 1200, res = 150)
plot(density(logcpm_pre[, 1], na.rm = TRUE),
     main = "Densidades por muestra (logCPM) - PRE TMM",
     xlab = "logCPM", ylab = "Densidad",
     xlim = xlim_all)
if (ncol(logcpm_pre) > 1) for (j in 2:ncol(logcpm_pre)) lines(density(logcpm_pre[, j], na.rm = TRUE))
dev.off()

png(file.path(fig_dir, paste0("04_density_logCPM_postTMM_", ts, ".png")),
    width = 1600, height = 1200, res = 150)
plot(density(logcpm_post[, 1], na.rm = TRUE),
     main = "Densidades por muestra (logCPM) - POST TMM",
     xlab = "logCPM", ylab = "Densidad",
     xlim = xlim_all)
if (ncol(logcpm_post) > 1) for (j in 2:ncol(logcpm_post)) lines(density(logcpm_post[, j], na.rm = TRUE))
dev.off()

cat("Density plots generados: PRE y POST TMM\n\n")

# ---- Tamaños de librería por muestra (pre/post TMM) ----
lib_size_pre <- colSums(counts_mat, na.rm = TRUE)

png(file.path(fig_dir, paste0("05_library_size_preTMM_", ts, ".png")),
    width = 2000, height = 900, res = 150)
barplot(lib_size_pre,
        main = "Library size por muestra - PRE TMM",
        ylab = "Total counts",
        names.arg = colnames(counts_mat),
        las = 2, cex.names = 0.6)
abline(h = median(lib_size_pre), col = "red", lty = 2)
dev.off()

png(file.path(fig_dir, paste0("06_library_size_postTMM_", ts, ".png")),
    width = 2000, height = 900, res = 150)
barplot(eff_lib,
        main = "Library size efectivo por muestra - POST TMM",
        ylab = "Effective library size (TMM)",
        names.arg = colnames(logcpm_post),
        las = 2, cex.names = 0.6)
abline(h = median(eff_lib), col = "red", lty = 2)
dev.off()

cat("Figuras de library size generadas: PRE y POST TMM\n\n")

# ---- Tablas de library size (pre/post TMM) ----
libsize_pre_df <- data.frame(
  sample = colnames(counts_mat),
  library_size = as.numeric(lib_size_pre),
  stringsAsFactors = FALSE
)

libsize_post_df <- data.frame(
  sample = colnames(logcpm_post),
  library_size_eff = as.numeric(eff_lib),
  stringsAsFactors = FALSE
)

lib_pre_fp <- file.path(tab_dir, paste0("library_size_preTMM_", ts, ".csv"))
lib_post_fp <- file.path(tab_dir, paste0("library_size_postTMM_", ts, ".csv"))

utils::write.csv(libsize_pre_df, lib_pre_fp, row.names = FALSE)
utils::write.csv(libsize_post_df, lib_post_fp, row.names = FALSE)

cat("Tablas de library size guardadas:\n")
cat("- Pre TMM:", lib_pre_fp, "\n")
cat("- Post TMM:", lib_post_fp, "\n\n")

# ---- MDS sobre logCPM ----
# Usamos limma::plotMDS para obtener coordenadas con plot=FALSE
# (top=500 usa genes/features más variables para distancias)
top_n <- 500
mds <- limma::plotMDS(logcpm_post, top = top_n, plot = FALSE)

mds_df <- data.frame(
  sample = colnames(logcpm_post),
  Dim1 = mds$x,
  Dim2 = mds$y,
  lib_size = as.numeric(lib_size),
  eff_lib = as.numeric(eff_lib),
  stringsAsFactors = FALSE
)

mds_fp <- file.path(tab_dir, paste0("mds_coordinates_logCPM_", ts, ".csv"))
utils::write.csv(mds_df, mds_fp, row.names = FALSE)

cat("MDS calculado con top =", top_n, "features.\n")
cat("Coordenadas guardadas en:", mds_fp, "\n\n")

# Plot MDS (con etiquetas)
png(file.path(fig_dir, paste0("02_MDS_logCPM_labels_", ts, ".png")), width = 1400, height = 1000, res = 150)
plot(mds$x, mds$y,
     xlab = "MDS1", ylab = "MDS2",
     main = paste0("MDS (logCPM, TMM) | top=", top_n))
text(mds$x, mds$y, labels = colnames(logcpm_post), pos = 3, cex = 0.7)
dev.off()

# Plot MDS (sin etiquetas, por si queda muy cargado)
png(file.path(fig_dir, paste0("03_MDS_logCPM_points_", ts, ".png")), width = 1400, height = 1000, res = 150)
plot(mds$x, mds$y,
     xlab = "MDS1", ylab = "MDS2",
     main = paste0("MDS (logCPM, TMM) | top=", top_n))
dev.off()

cat("============================================\n")
cat("FIN\n")
cat("- Log:", log_fp, "\n")
cat("- Figuras:", fig_dir, "\n")
cat("- Tabla MDS:", mds_fp, "\n")
type_info <- NA
try({
  ft_map <- feature_type_df
  ft_map <- ft_map[match(rownames(counts_mat), ft_map$id), , drop = FALSE]
  type_tab <- sort(table(ft_map$type), decreasing = TRUE)
  type_info <- type_tab
  cat("- Conteo por tipo (features retenidas):\n")
  print(type_tab)
}, silent = TRUE)
cat("============================================\n")

sink(type = "message")
sink(type = "output")
close(zz)
