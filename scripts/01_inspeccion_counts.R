#!/usr/bin/env Rscript

# ============================================================
# 01_inspeccion_counts.R (corregido)
# - Excluye columnas de anotación (incluye entrez_id) del bloque de conteos
# - Library size "oficial" con edgeR::DGEList
# - Log -> logs/
# - Figuras -> results/figures/inspeccion/
# ============================================================

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

in_fp <- if (length(args) >= 1) args[1] else file.path("data", "intermediate", "Gliomas_all_counts_merged.csv")

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path("results", "figures", "inspeccion")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_fp <- file.path("logs", paste0("inspeccion_counts_", ts, ".txt"))

zz <- file(log_fp, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat("============================================\n")
cat("Inspección de matriz de conteos (edgeR)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Input:", in_fp, "\n")
cat("Figuras:", fig_dir, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

if (!file.exists(in_fp)) {
  cat("ERROR: No existe el archivo de entrada:\n", in_fp, "\n", sep = "")
  sink(type = "message"); sink(type = "output")
  close(zz)
  stop("Archivo de entrada no encontrado.")
}

cat("Primeras 5 líneas (crudo):\n")
first_lines <- readLines(in_fp, n = 5, warn = FALSE)
cat(paste0(first_lines, collapse = "\n"), "\n\n")

read_counts_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    cat("Leyendo con data.table::fread\n\n")
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    cat("Leyendo con utils::read.csv\n\n")
    utils::read.csv(path, check.names = FALSE)
  }
}

df <- read_counts_table(in_fp)

cat("Dimensiones (df):", nrow(df), "filas x", ncol(df), "columnas\n\n")
cat("Columnas (primeras 30):\n")
print(utils::head(colnames(df), 30))
cat("\nVista previa (3 filas x 8 cols):\n")
print(df[seq_len(min(3, nrow(df))), seq_len(min(8, ncol(df))), drop = FALSE])
cat("\n")

# ---- Identificar columnas de anotación por nombre (robusto) ----
# Ajusta aquí si tu archivo cambia.
anno_candidates <- c("id", "type", "entrez_id", "HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))

cat("Columnas de anotación detectadas:\n")
print(anno_present)
cat("\n")

if (!("id" %in% colnames(df))) {
  cat("WARNING: No veo columna 'id'. Usaré rownames numéricos.\n\n")
}

# ---- Separar IDs de feature y matriz de conteos ----
feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))

# Conteos = todas las columnas EXCEPTO anotación detectada
count_cols <- setdiff(colnames(df), anno_present)

# ---- Filtrar muestras: eliminar sanos post mortem (columnas que empiezan con "A") ----
before_n <- length(count_cols)
drop_samples <- grep("^A", count_cols, value = TRUE)

cat("Filtrado de muestras (drop '^A'):\n")
cat("- Muestras antes:", before_n, "\n")
cat("- A eliminar:", length(drop_samples), "\n")
if (length(drop_samples) > 0) {
  cat("- Nombres eliminados:\n")
  print(drop_samples)
}

count_cols <- setdiff(count_cols, drop_samples)
cat("- Muestras después:", length(count_cols), "\n\n")

if (length(count_cols) == 0) {
  sink(type = "message"); sink(type = "output")
  close(zz)
  stop("Tras eliminar columnas que empiezan con 'A' no quedan muestras para inspección.")
}

counts_df <- df[, count_cols, drop = FALSE]

# ---- Coerción segura a numérico (por si alguna columna llega como character) ----
to_numeric <- function(x) suppressWarnings(as.numeric(x))
counts_df[] <- lapply(counts_df, to_numeric)

# Reportar columnas que quedaron completamente NA tras coerción (típico de errores de parseo)
all_na_cols <- names(counts_df)[vapply(counts_df, function(x) all(is.na(x)), logical(1))]
if (length(all_na_cols) > 0) {
  cat("WARNING: Estas columnas quedaron 100% NA tras coerción a numérico:\n")
  print(all_na_cols)
  cat("\n")
}

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"

rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

cat("Dimensiones (counts_mat):", nrow(counts_mat), "features x", ncol(counts_mat), "muestras\n\n")

# ---- Chequeos de integridad ----
na_total <- sum(is.na(counts_mat))
neg_total <- sum(counts_mat < 0, na.rm = TRUE)
nonint_total <- sum(abs(counts_mat - round(counts_mat)) > 1e-6, na.rm = TRUE)

cat("Chequeos de integridad:\n")
cat("- NA totales:", na_total, "\n")
cat("- Negativos:", neg_total, "\n")
cat("- No enteros (tol 1e-6):", nonint_total, "\n\n")

# ---- Library size "oficial" con edgeR ----
if (!requireNamespace("edgeR", quietly = TRUE)) {
  cat("ERROR: edgeR no está instalado/cargable en este entorno.\n")
  cat("Instálalo en omics-R (BiocManager::install('edgeR')).\n")
  sink(type = "message"); sink(type = "output")
  close(zz)
  stop("edgeR no disponible.")
}

cat("edgeR version:", as.character(utils::packageVersion("edgeR")), "\n\n")

y <- edgeR::DGEList(counts = counts_mat)
lib_size <- y$samples$lib.size  # == colSums(counts_mat), pero “oficial” edgeR

cat("Resumen de library sizes (edgeR::DGEList -> y$samples$lib.size):\n")
print(summary(lib_size))
cat("\nTop 10 menor library size:\n")
print(head(sort(lib_size), 10))
cat("\nTop 10 mayor library size:\n")
print(head(sort(lib_size, decreasing = TRUE), 10))
cat("\n")

# (Opcional de inspección) Effective library size si luego se calculan norm factors
y_norm <- edgeR::calcNormFactors(y)
eff_lib <- y_norm$samples$lib.size * y_norm$samples$norm.factors
cat("Resumen de effective library size (lib.size * norm.factors; edgeR TMM):\n")
print(summary(eff_lib))
cat("\n")

# ---- Fracción de ceros por muestra (solo sobre conteos reales) ----
zero_frac <- colMeans(counts_mat == 0, na.rm = TRUE)
cat("Resumen fracción de ceros por muestra:\n")
print(summary(zero_frac))
cat("\nTop 10 mayor fracción de ceros:\n")
print(head(sort(zero_frac, decreasing = TRUE), 10))
cat("\n")

# ---- Estadísticas por feature ----
feature_total <- rowSums(counts_mat, na.rm = TRUE)
cat("Resumen de conteo total por feature (rowSums):\n")
print(summary(feature_total))
cat("\nTop 10 features más abundantes:\n")
print(head(sort(feature_total, decreasing = TRUE), 10))
cat("\n")

# ---- Figuras ----
png(file.path(fig_dir, paste0("01_library_sizes_edgeR_", ts, ".png")), width = 1400, height = 900, res = 150)
hist(lib_size, main = "Library sizes (edgeR DGEList)", xlab = "Total counts por muestra", las = 1)
dev.off()

png(file.path(fig_dir, paste0("02_zero_fraction_", ts, ".png")), width = 1400, height = 900, res = 150)
hist(zero_frac, main = "Fracción de ceros por muestra", xlab = "mean(count==0)", las = 1)
dev.off()

box_fp <- file.path(fig_dir, paste0("03_boxplot_log1p_counts_", ts, ".png"))
ok_box <- TRUE
tryCatch({
  png(box_fp, width = 2000, height = 900, res = 150)
  boxplot(log1p(counts_mat),
          main = "log1p(counts) por muestra",
          ylab = "log1p(count)",
          las = 2, cex.axis = 0.5)
  dev.off()
}, error = function(e) {
  ok_box <<- FALSE
  cat("WARNING: No pude generar el boxplot:\n")
  cat("  ", conditionMessage(e), "\n\n")
  try(dev.off(), silent = TRUE)
})

# ---- Guardar objetos clave ----
rds_fp <- file.path("data", "intermediate", paste0("inspect_counts_objects_", ts, ".rds"))
saveRDS(
  list(
    input = in_fp,
    anno_present = anno_present,
    counts_dim = dim(counts_mat),
    lib_size = lib_size,
    eff_lib = eff_lib,
    zero_frac = zero_frac,
    feature_total = feature_total
  ),
  file = rds_fp
)

cat("============================================\n")
cat("FIN\n")
cat("- Log:", log_fp, "\n")
cat("- Figuras:", fig_dir, "\n")
cat("- RDS:", rds_fp, "\n")
if (!ok_box) cat("- Nota: boxplot no se generó.\n")
cat("============================================\n")

sink(type = "message")
sink(type = "output")
close(zz)
