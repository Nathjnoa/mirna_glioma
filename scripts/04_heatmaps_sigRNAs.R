#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

read_table <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.delim(path, check.names = FALSE)
  }
}

# ---- CLI ----
counts_fp <- get_arg("--counts", file.path("data","intermediate","Gliomas_all_counts_merged.csv"))
meta_fp   <- get_arg("--meta",   file.path("data","intermediate","Metadatos_gliomas_verificados.csv"))
spec_fp   <- get_arg("--spec",   file.path("config","de_specs.csv"))
de_dir    <- get_arg("--de_dir", file.path("results","DE"))

drop_re   <- get_arg("--drop_regex", "^A")
min_lib   <- as.numeric(get_arg("--min_libsize", "100000"))
prior_cnt <- as.numeric(get_arg("--prior_count", "2"))

fdrs_str  <- get_arg("--fdrs", "0.05,0.1")
fdrs <- as.numeric(trimws(strsplit(fdrs_str, ",")[[1]]))
fdrs <- fdrs[is.finite(fdrs) & fdrs > 0 & fdrs < 1]
if (length(fdrs) == 0) stop("No hay FDRs válidos en --fdrs")

max_features <- as.integer(get_arg("--max_features", "200"))  # si hay demasiados, toma top por FDR
width_px  <- as.integer(get_arg("--width",  "2000"))
height_px <- as.integer(get_arg("--height", "1600"))
res_dpi   <- as.integer(get_arg("--res",    "150"))

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
fig_root <- file.path("results","figures","heatmaps")
dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("heatmaps_sigRNAs_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output"); sink(zz, type = "message")

cat("============================================\n")
cat("Heatmaps de RNAs significativos (Z-score)\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Counts:", counts_fp, "\n")
cat("Meta:", meta_fp, "\n")
cat("Spec:", spec_fp, "\n")
cat("DE dir:", de_dir, "\n")
cat("drop_regex:", drop_re, "\n")
cat("min_libsize:", min_lib, "\n")
cat("prior_count:", prior_cnt, "\n")
cat("FDRs:", paste(fdrs, collapse=", "), "\n")
cat("max_features:", max_features, "\n")
cat("Fig root:", fig_root, "\n")
cat("============================================\n\n")

# ---- deps ----
has_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
has_circ <- requireNamespace("circlize", quietly = TRUE)
if (!has_ch || !has_circ) {
  stop("Falta instalar ComplexHeatmap (y circlize) en este entorno.")
}
suppressPackageStartupMessages({
  library(edgeR)
  library(ComplexHeatmap)
  library(circlize)
})

# ---- Leer spec ----
if (!file.exists(spec_fp)) stop("No existe spec: ", spec_fp)
spec <- read_table(spec_fp)

required_cols <- c("analysis_id","mode","var","event_var","cut","cut_type","case","ref","collapse_map","exclude_rule")
missing_cols <- setdiff(required_cols, colnames(spec))
if (length(missing_cols) > 0) stop("Faltan columnas en spec: ", paste(missing_cols, collapse=", "))

spec$analysis_id <- as.character(spec$analysis_id)
spec$mode <- as.character(spec$mode)

if (anyDuplicated(spec$analysis_id)) {
  dups <- unique(spec$analysis_id[duplicated(spec$analysis_id)])
  stop("analysis_id duplicados en spec: ", paste(dups, collapse=", "))
}

cat("Analyses en spec:", nrow(spec), "\n\n")

# ---- Leer counts y meta base ----
if (!file.exists(counts_fp)) stop("No existe counts: ", counts_fp)
if (!file.exists(meta_fp)) stop("No existe meta: ", meta_fp)

df <- if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fread(counts_fp, data.table = FALSE, check.names = FALSE)
} else read.csv(counts_fp, check.names = FALSE)

meta <- if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fread(meta_fp, data.table = FALSE, check.names = FALSE)
} else read.csv(meta_fp, check.names = FALSE)

# limpiar filas meta con id vacío por si acaso
if ("id" %in% colnames(meta)) {
  meta <- meta[!(is.na(meta$id) | trimws(as.character(meta$id)) == ""), , drop = FALSE]
} else {
  stop("El metadata debe tener columna 'id'.")
}

# columnas de anotación en counts
anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop post mortem por regex
drop_samples <- grep(drop_re, count_cols, value = TRUE)
count_cols <- setdiff(count_cols, drop_samples)

feature_id <- if ("id" %in% colnames(df)) df[["id"]] else seq_len(nrow(df))
counts_df <- df[, count_cols, drop = FALSE]
counts_df[] <- lapply(counts_df, function(x) suppressWarnings(as.numeric(x)))

counts_mat <- as.matrix(counts_df)
storage.mode(counts_mat) <- "numeric"
rownames(counts_mat) <- make.unique(as.character(feature_id))
colnames(counts_mat) <- make.unique(colnames(counts_mat))

# limpieza mínima features
counts_mat <- counts_mat[rowSums(counts_mat, na.rm = TRUE) > 0, , drop = FALSE]
rv <- apply(counts_mat, 1, function(x) stats::var(x, na.rm = TRUE))
counts_mat <- counts_mat[is.finite(rv) & rv > 0, , drop = FALSE]

cat("Base counts_mat (tras drop_regex y limpieza mínima): ",
    nrow(counts_mat), "features x ", ncol(counts_mat), "samples\n\n", sep="")

# filtro global por libsize
y0 <- DGEList(counts_mat)
lib0 <- y0$samples$lib.size
keep_lib <- lib0 >= min_lib
cat("Filtro global lib_size >= ", min_lib, ":\n", sep="")
cat("- samples antes:", length(lib0), "\n")
cat("- samples removidas:", sum(!keep_lib), "\n")
cat("- samples después:", sum(keep_lib), "\n\n")
counts_mat <- counts_mat[, keep_lib, drop = FALSE]

samples_final <- colnames(counts_mat)
meta_ids <- as.character(meta$id)

missing_in_meta <- setdiff(samples_final, meta_ids)
if (length(missing_in_meta) > 0) {
  cat("WARNING: muestras en counts pero no en meta$id. Se descartarán de counts:\n",
      paste(missing_in_meta, collapse = ", "), "\n")
  keep_cols <- setdiff(samples_final, missing_in_meta)
  if (length(keep_cols) == 0) {
    stop("No quedan muestras tras descartar las ausentes en metadata.")
  }
  counts_mat <- counts_mat[, keep_cols, drop = FALSE]
  samples_final <- colnames(counts_mat)
}
meta_base <- meta[match(samples_final, meta_ids), , drop = FALSE]
stopifnot(all(as.character(meta_base$id) == samples_final))

# ---- helpers para derivar grupo según spec ----
as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

parse_collapse_map_levels <- function(map_str) {
  # "0|1=01;2=2" -> list(c("0","1")->"01", c("2")->"2")
  out <- list()
  parts <- strsplit(map_str, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  for (p in parts) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    lhs <- trimws(kv[1]); rhs <- trimws(kv[2])
    lhs_levels <- trimws(strsplit(lhs, "\\|")[[1]])
    out[[rhs]] <- lhs_levels
  }
  out
}

parse_range_map <- function(map_str) {
  # "40-60=age40_60;71-200=gt70" -> list(list(min=40,max=60,label="age40_60"), ...)
  out <- list()
  parts <- strsplit(map_str, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  for (p in parts) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    rng <- trimws(kv[1]); lab <- trimws(kv[2])
    mm <- strsplit(rng, "-", fixed = TRUE)[[1]]
    if (length(mm) != 2) next
    mn <- as_num_safe(mm[1]); mx <- as_num_safe(mm[2])
    if (!is.finite(mn) || !is.finite(mx)) next
    out[[length(out) + 1]] <- list(min = mn, max = mx, label = lab)
  }
  out
}

derive_group <- function(spec_row, meta_sub) {
  mode <- spec_row$mode
  var <- spec_row$var
  event_var <- spec_row$event_var
  cut <- spec_row$cut
  cut_type <- spec_row$cut_type
  case <- spec_row$case
  ref <- spec_row$ref
  collapse_map <- spec_row$collapse_map
  exclude_rule <- spec_row$exclude_rule

  if (!(var %in% colnames(meta_sub))) stop("No existe columna en meta: ", var)

  if (mode == "as_is") {
    g <- as.character(meta_sub[[var]])
    g[trimws(g) == ""] <- NA
    g <- as.factor(g)
    # restringir a case/ref si están definidos
    if (nzchar(case) && nzchar(ref)) {
      keep <- as.character(g) %in% c(case, ref)
      g[!keep] <- NA
      g <- droplevels(g)
    }
    return(list(group = g, x_cont = NULL))

  } else if (mode == "binary_cut") {
    x <- as_num_safe(meta_sub[[var]])
    if (!is.finite(as.numeric(cut))) stop("cut no numérico en binary_cut: ", cut)
    thr <- as.numeric(cut)
    if (cut_type == "gt") {
      g <- ifelse(x > thr, case, ref)
    } else if (cut_type == "ge") {
      g <- ifelse(x >= thr, case, ref)
    } else {
      stop("cut_type inválido en binary_cut: ", cut_type)
    }
    g[!is.finite(x)] <- NA
    return(list(group = as.factor(g), x_cont = NULL))

  } else if (mode == "collapse_levels") {
    if (!nzchar(collapse_map)) stop("collapse_map vacío en collapse_levels")
    mapping <- parse_collapse_map_levels(collapse_map)
    raw <- as.character(meta_sub[[var]])
    raw[trimws(raw) == ""] <- NA
    g <- rep(NA_character_, length(raw))
    for (lab in names(mapping)) {
      g[raw %in% mapping[[lab]]] <- lab
    }
    return(list(group = as.factor(g), x_cont = NULL))

  } else if (mode == "survival_cut") {
    if (!(event_var %in% colnames(meta_sub))) stop("No existe event_var en meta: ", event_var)
    tmo <- as_num_safe(meta_sub[[var]])
    ev  <- as_num_safe(meta_sub[[event_var]])
    if (!is.finite(as.numeric(cut))) stop("cut no numérico en survival_cut: ", cut)
    thr <- as.numeric(cut)

    # grupos base: ge vs lt (definidos por case/ref en spec)
    g <- rep(NA_character_, length(tmo))

    # lt: evento==1 y t < cut
    g[is.finite(tmo) & is.finite(ev) & ev == 1 & tmo < thr] <- ref  # ref = lt24 (por spec)
    # ge: t >= cut (evento 0 o 1)
    g[is.finite(tmo) & tmo >= thr] <- case  # case = ge24

    # exclusión censurados < cut
    if (exclude_rule == "exclude_if_censored_lt_cut") {
      cens_lt <- is.finite(tmo) & is.finite(ev) & ev == 0 & tmo < thr
      g[cens_lt] <- NA
    }

    return(list(group = as.factor(g), x_cont = NULL))

  } else if (mode == "range_vs") {
    if (!nzchar(collapse_map)) stop("collapse_map vacío en range_vs")
    ranges <- parse_range_map(collapse_map)
    x <- as_num_safe(meta_sub[[var]])
    g <- rep(NA_character_, length(x))
    for (rr in ranges) {
      idx <- is.finite(x) & x >= rr$min & x <= rr$max
      g[idx] <- rr$label
    }
    # restringir a case/ref
    if (nzchar(case) && nzchar(ref)) {
      keep <- g %in% c(case, ref)
      g[!keep] <- NA
    }
    return(list(group = as.factor(g), x_cont = NULL))

  } else if (mode == "continuous") {
    x <- as_num_safe(meta_sub[[var]])
    if (!is.finite(as.numeric(cut))) stop("cut no numérico en continuous (escala): ", cut)
    scale_by <- as.numeric(cut)
    x_scaled <- x / scale_by
    x_scaled[!is.finite(x_scaled)] <- NA
    return(list(group = NULL, x_cont = x_scaled))

  } else {
    stop("mode no soportado: ", mode)
  }
}

# ---- compute base logCPM post-TMM (para toda la base) ----
y_base <- DGEList(counts_mat)
y_base <- calcNormFactors(y_base, method = "TMM")
logcpm_base <- cpm(y_base, log = TRUE, prior.count = prior_cnt)

# ---- Heatmap function ----
make_z <- function(mat) {
  # Z-score por fila; elimina filas con SD==0
  z <- t(scale(t(mat)))
  z <- z[apply(z, 1, function(x) all(is.finite(x))), , drop = FALSE]
  z
}

plot_heatmap <- function(z, annotation_col, out_png, main_title) {
  cluster_rows <- nrow(z) > 1
  hc_rows <- NULL
  if (cluster_rows) {
    d_rows <- stats::as.dist(1 - stats::cor(t(z), use = "pairwise.complete.obs", method = "pearson"))
    hc_rows <- stats::hclust(d_rows, method = "complete")
  }

  top_ha <- NULL
  if (!is.null(annotation_col)) {
    top_ha <- ComplexHeatmap::HeatmapAnnotation(df = annotation_col)
  }

  png(out_png, width = width_px, height = height_px, res = res_dpi)
  ht <- ComplexHeatmap::Heatmap(
    z,
    name = "Z",
    cluster_rows = if (cluster_rows) hc_rows else FALSE,
    cluster_columns = FALSE,
    show_row_names = (nrow(z) <= 80),
    show_column_names = TRUE,
    top_annotation = top_ha,
    column_title = main_title
  )
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# ---- Main loop ----
for (i in seq_len(nrow(spec))) {
  s <- spec[i, , drop = FALSE]
  analysis_id <- as.character(s$analysis_id)
  analysis_fs <- safe_name(analysis_id)

  cat("--------------------------------------------------\n")
  cat("Analysis:", analysis_id, " (mode=", s$mode, ")\n", sep="")
  cat("var:", s$var, "\n")
  cat("case/ref:", s$case, "/", s$ref, "\n")
  cat("cut:", s$cut, " cut_type:", s$cut_type, " exclude_rule:", s$exclude_rule, "\n")
  cat("--------------------------------------------------\n")

  # localizar tabla DE más reciente para este análisis
  analysis_dir <- file.path(de_dir, analysis_fs)
  if (!dir.exists(analysis_dir)) {
    cat("SKIP: no existe carpeta de DE:", analysis_dir, "\n\n")
    next
  }

  # buscar tsv "all"
  all_files <- list.files(analysis_dir, pattern = "DE_.*_all_.*\\.tsv$", full.names = TRUE)
  if (length(all_files) == 0) {
    cat("SKIP: no se encontró tabla *_all_*.tsv en:", analysis_dir, "\n\n")
    next
  }
  # elegir el más reciente por mtime
  all_fp <- all_files[order(file.info(all_files)$mtime, decreasing = TRUE)][1]
  cat("Usando DE table:", all_fp, "\n")

  de <- read.delim(all_fp, check.names = FALSE)
  if (!("feature_id" %in% colnames(de))) {
    cat("SKIP: la tabla no tiene columna feature_id\n\n")
    next
  }
  if (!("FDR" %in% colnames(de))) {
    cat("SKIP: la tabla no tiene columna FDR\n\n")
    next
  }

  # derivar subset de muestras y anotación
  der <- derive_group(s, meta_base)
  mode <- as.character(s$mode)

  if (mode == "continuous") {
    x_cont <- der$x_cont
    keep_samp <- is.finite(x_cont)
    if (sum(keep_samp) < 6) {
      cat("SKIP: muy pocas muestras con valor continuo válido (n=", sum(keep_samp), ")\n\n", sep="")
      next
    }
    samp_ids <- meta_base$id[keep_samp]
    # ordenar por edad
    ord <- order(x_cont[keep_samp], decreasing = FALSE)
    samp_ids <- samp_ids[ord]
    ann <- data.frame(value = x_cont[keep_samp][ord], row.names = samp_ids)
    ann_col <- ann

  } else {
    g <- der$group
    keep_samp <- !is.na(g) & trimws(as.character(g)) != ""
    g <- droplevels(g[keep_samp])
    if (nlevels(g) < 2) {
      cat("SKIP: <2 niveles en group tras aplicar reglas\n\n")
      next
    }
    tab <- table(g)
    cat("Muestras por grupo:\n")
    print(tab)
    if (any(tab < 2)) {
      cat("SKIP: algún grupo tiene <2 muestras\n\n")
      next
    }
    samp_ids <- meta_base$id[keep_samp]
    # ordenar por grupo para visualizar separación
    ord <- order(as.character(g))
    samp_ids <- samp_ids[ord]
    g <- g[ord]
    ann_col <- data.frame(group = g, row.names = samp_ids)
  }

  if (!is.null(ann_col) && "group" %in% colnames(ann_col)) {
    ann_col$group <- as.factor(ann_col$group)
  }

  # expresión logCPM para esas muestras
  present_samp <- intersect(samp_ids, colnames(logcpm_base))
  if (length(present_samp) != length(samp_ids)) {
    cat("WARNING: algunas muestras no están en logcpm_base (revisar IDs)\n")
    samp_ids <- present_samp
    if (!is.null(ann_col)) ann_col <- ann_col[samp_ids, , drop = FALSE]
  }

  # iterar umbrales FDR
  for (thr in sort(fdrs)) {
    sig <- de[is.finite(as.numeric(de$FDR)) & de$FDR < thr, , drop = FALSE]
    n_sig <- nrow(sig)
    cat(sprintf("FDR<%.2g: n_sig=%d\n", thr, n_sig))

    if (n_sig == 0) next

    # si demasiados, reducir a top por FDR
    sig$FDR_num <- as.numeric(sig$FDR)
    sig <- sig[order(sig$FDR_num, na.last = TRUE), , drop = FALSE]
    if (n_sig > max_features) {
      sig <- sig[seq_len(max_features), , drop = FALSE]
      cat("  -> reducido a top ", max_features, " por FDR\n", sep="")
    }

    feats <- unique(as.character(sig$feature_id))
    feats <- feats[feats %in% rownames(logcpm_base)]
    mat <- logcpm_base[feats, samp_ids, drop = FALSE]
    z <- make_z(mat)
    if (nrow(z) < 1) {
      cat("  -> SKIP heatmap: sin features tras Z-score (SD=0/NA)\n")
      next
    }

    out_dir <- file.path(fig_root, analysis_fs)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    out_png <- file.path(out_dir, paste0("heatmap_Z_FDR", gsub("\\.", "", sprintf("%.2g", thr)), "_", ts, ".png"))
    main_title <- if (mode == "continuous") {
      paste0(analysis_id, " | Z-score (logCPM) | FDR<", thr, " | ordered by ", s$var)
    } else {
      paste0(analysis_id, " | Z-score (logCPM) | FDR<", thr)
    }

    plot_heatmap(z, ann_col, out_png, main_title)

    # guardar lista de features
    out_list <- file.path(out_dir, paste0("features_FDR", gsub("\\.", "", sprintf("%.2g", thr)), "_", ts, ".txt"))
    writeLines(rownames(z), con = out_list)

    cat("  -> saved:", out_png, "\n")
  }

  cat("\n")
}

cat("FIN. Log:", log_fp, "\n")
sink(type="message"); sink(type="output"); close(zz)
