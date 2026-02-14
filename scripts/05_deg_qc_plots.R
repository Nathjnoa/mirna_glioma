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

as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

# ---- CLI defaults ----
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

max_features <- as.integer(get_arg("--max_features", "30"))  # top N por FDR si hay demasiados
make_pdf <- tolower(get_arg("--pdf", "true")) %in% c("true","t","1","yes","y")
label_outliers <- tolower(get_arg("--label_outliers", "false")) %in% c("true","t","1","yes","y")
n_label <- as.integer(get_arg("--n_label", "2"))

width_px  <- as.integer(get_arg("--width",  "3000"))
height_px <- as.integer(get_arg("--height", "2400"))
res_dpi   <- as.integer(get_arg("--res",    "300"))

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

out_root <- get_arg("--out_root", file.path("results","figures","deg_qc"))
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("deg_qc_plots_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output"); sink(zz, type = "message")

cat("============================================\n")
cat("QC plots por DEG (box/jitter y scatter continuo)\n")
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
cat("make_pdf:", make_pdf, "\n")
cat("label_outliers:", label_outliers, " n_label:", n_label, "\n")
cat("Out root:", out_root, "\n")
cat("Log:", log_fp, "\n")
cat("============================================\n\n")

suppressPackageStartupMessages({
  library(edgeR)
})

# ---- Spec ----
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

# ---- Read counts/meta ----
if (!file.exists(counts_fp)) stop("No existe counts: ", counts_fp)
if (!file.exists(meta_fp)) stop("No existe meta: ", meta_fp)

df <- if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fread(counts_fp, data.table = FALSE, check.names = FALSE)
} else read.csv(counts_fp, check.names = FALSE)

meta <- if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fread(meta_fp, data.table = FALSE, check.names = FALSE)
} else read.csv(meta_fp, check.names = FALSE)

if (!("id" %in% colnames(meta))) stop("El metadata debe tener columna 'id'.")
meta <- meta[!(is.na(meta$id) | trimws(as.character(meta$id)) == ""), , drop = FALSE]

anno_candidates <- c("id","type","entrez_id","HGNC_symbol")
anno_present <- intersect(anno_candidates, colnames(df))
count_cols <- setdiff(colnames(df), anno_present)

# drop postmortem
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

# filtro global libsize
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

# ---- base logCPM post-TMM ----
y_base <- DGEList(counts_mat)
y_base <- calcNormFactors(y_base, method = "TMM")
logcpm_base <- cpm(y_base, log = TRUE, prior.count = prior_cnt)

# ---- helpers derive group/x ----
parse_collapse_map_levels <- function(map_str) {
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

derive_samples_and_annotation <- function(spec_row, meta_sub) {
  mode <- as.character(spec_row$mode)
  var <- as.character(spec_row$var)
  event_var <- as.character(spec_row$event_var)
  cut <- as_num_safe(spec_row$cut)
  cut_type <- as.character(spec_row$cut_type)
  case <- as.character(spec_row$case)
  ref <- as.character(spec_row$ref)
  collapse_map <- as.character(spec_row$collapse_map)
  exclude_rule <- as.character(spec_row$exclude_rule)

  if (!(var %in% colnames(meta_sub))) stop("No existe columna en meta: ", var)

  if (mode == "continuous") {
    x <- as_num_safe(meta_sub[[var]])
    keep <- is.finite(x)
    if (exclude_rule == "drop_na") keep <- keep & is.finite(x)
    if (!is.finite(cut) || cut <= 0) stop("cut inválido para continuous (escala): ", spec_row$cut)
    x_scaled <- x / cut
    keep <- keep & is.finite(x_scaled)

    samp_ids <- meta_sub$id[keep]
    ord <- order(x[keep], decreasing = FALSE)
    samp_ids <- samp_ids[ord]

    ann <- data.frame(value = x[keep][ord], row.names = samp_ids)
    return(list(samp_ids = samp_ids, mode = mode, ann = ann, group = NULL, x = x[keep][ord], var = var))

  } else {
    g <- NULL

    if (mode == "as_is") {
      g <- as.character(meta_sub[[var]])
      g[trimws(g) == ""] <- NA
      # restringir a case/ref si definidos
      if (nzchar(case) && nzchar(ref)) {
        keep2 <- g %in% c(case, ref)
        g[!keep2] <- NA
      }

    } else if (mode == "binary_cut") {
      x <- as_num_safe(meta_sub[[var]])
      if (!is.finite(cut)) stop("cut inválido en binary_cut: ", spec_row$cut)
      if (cut_type == "gt") {
        g <- ifelse(x > cut, case, ref)
      } else if (cut_type == "ge") {
        g <- ifelse(x >= cut, case, ref)
      } else stop("cut_type inválido en binary_cut: ", cut_type)
      g[!is.finite(x)] <- NA

    } else if (mode == "collapse_levels") {
      if (!nzchar(collapse_map)) stop("collapse_map vacío en collapse_levels")
      mapping <- parse_collapse_map_levels(collapse_map)
      raw <- as.character(meta_sub[[var]])
      raw[trimws(raw) == ""] <- NA
      g <- rep(NA_character_, length(raw))
      for (lab in names(mapping)) g[raw %in% mapping[[lab]]] <- lab

    } else if (mode == "survival_cut") {
      if (!(event_var %in% colnames(meta_sub))) stop("No existe event_var en meta: ", event_var)
      tmo <- as_num_safe(meta_sub[[var]])
      ev  <- as_num_safe(meta_sub[[event_var]])
      if (!is.finite(cut)) stop("cut inválido en survival_cut: ", spec_row$cut)

      g <- rep(NA_character_, length(tmo))
      # ref = lt (evento antes del corte)
      g[is.finite(tmo) & is.finite(ev) & ev == 1 & tmo < cut] <- ref
      # case = ge (>= corte)
      g[is.finite(tmo) & tmo >= cut] <- case

      if (exclude_rule == "exclude_if_censored_lt_cut") {
        cens_lt <- is.finite(tmo) & is.finite(ev) & ev == 0 & tmo < cut
        g[cens_lt] <- NA
      }

    } else if (mode == "range_vs") {
      if (!nzchar(collapse_map)) stop("collapse_map vacío en range_vs")
      x <- as_num_safe(meta_sub[[var]])
      ranges <- parse_range_map(collapse_map)
      g <- rep(NA_character_, length(x))
      for (rr in ranges) {
        idx <- is.finite(x) & x >= rr$min & x <= rr$max
        g[idx] <- rr$label
      }
      if (nzchar(case) && nzchar(ref)) {
        keep2 <- g %in% c(case, ref)
        g[!keep2] <- NA
      }

    } else {
      stop("mode no soportado aquí: ", mode)
    }

    if (exclude_rule == "drop_na" || exclude_rule == "drop_unmapped") {
      # ya quedaron NAs; solo mantenemos no-NA
    }

    keep <- !is.na(g) & trimws(g) != ""
    g <- as.factor(g[keep])

    samp_ids <- meta_sub$id[keep]
    ord <- order(as.character(g))
    samp_ids <- samp_ids[ord]
    g <- g[ord]
    ann <- data.frame(group = g, row.names = samp_ids)

    return(list(samp_ids = samp_ids, mode = mode, ann = ann, group = g, x = NULL, var = var))
  }
}

# ---- Publication-quality theme ----
suppressPackageStartupMessages({
  library(ggplot2)
})

theme_pub_qc <- function() {
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 12, color = "#555555", hjust = 0),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
      plot.margin = margin(10, 15, 10, 10, unit = "mm")
    )
}

# Colorblind-safe palette
colors_qc <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9")

# ---- plotting ----
plot_box_jitter <- function(y, group, title, subtitle = NULL, out_png) {
  samp_ids <- names(y)
  if (is.null(samp_ids)) samp_ids <- paste0("s", seq_along(y))

  df_plot <- data.frame(
    sample_id = samp_ids,
    expression = as.numeric(y),
    group = as.factor(group),
    stringsAsFactors = FALSE
  )

  # Count per group
  n_per_group <- table(df_plot$group)
  group_labels <- paste0(names(n_per_group), "\n(n=", as.integer(n_per_group), ")")
  names(group_labels) <- names(n_per_group)

  p <- ggplot(df_plot, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, linewidth = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2, color = "#333333") +
    scale_fill_manual(values = colors_qc, guide = "none") +
    scale_x_discrete(labels = group_labels) +
    labs(
      x = NULL,
      y = "logCPM (TMM normalized)",
      title = title,
      subtitle = if (!is.null(subtitle)) subtitle else NULL
    ) +
    coord_cartesian(ylim = c(4, 16)) +
    theme_pub_qc()

  # Add sample labels if few samples
  if (nrow(df_plot) <= 30) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = sample_id),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.2
      )
    }
  }

  ggsave(out_png, p, width = 10, height = 8, dpi = 300, bg = "white")
}

plot_scatter_cont <- function(x, y, feat, title, out_png, label_outliers = FALSE, n_label = 2) {
  samp_ids <- names(y)
  if (is.null(samp_ids)) samp_ids <- paste0("s", seq_along(y))

  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(invisible(NULL))

  df_plot <- data.frame(
    sample_id = samp_ids[ok],
    x = as.numeric(x[ok]),
    y = as.numeric(y[ok]),
    stringsAsFactors = FALSE
  )

  # Correlation
  sp <- suppressWarnings(stats::cor.test(df_plot$x, df_plot$y, method = "spearman", exact = FALSE))
  rho <- unname(sp$estimate)
  pval <- sp$p.value

  # Format p-value
  p_txt <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)

  p <- ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(alpha = 0.8, size = 2.5, color = "#0072B2") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00", fill = "#D55E0030", linewidth = 1) +
    labs(
      x = title,
      y = "logCPM (TMM normalized)",
      title = feat,
      subtitle = sprintf("Spearman rho = %.3f, %s, n = %d", rho, p_txt, nrow(df_plot))
    ) +
    coord_cartesian(ylim = c(4, 16)) +
    theme_pub_qc()

  # Add sample labels if few samples
  if (nrow(df_plot) <= 30) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = sample_id),
        size = 3,
        max.overlaps = 15,
        box.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.2
      )
    }
  }

  ggsave(out_png, p, width = 10, height = 8, dpi = 300, bg = "white")
  invisible(list(rho = rho, p = pval, n = nrow(df_plot)))
}

# ---- loop analyses ----
for (i in seq_len(nrow(spec))) {
  s <- spec[i, , drop = FALSE]
  analysis_id <- as.character(s$analysis_id)
  analysis_fs <- safe_name(analysis_id)
  mode <- as.character(s$mode)

  cat("--------------------------------------------------\n")
  cat("Analysis:", analysis_id, " mode=", mode, "\n", sep="")
  cat("var:", s$var, "\n")
  cat("--------------------------------------------------\n")

  analysis_dir <- file.path(de_dir, analysis_fs)
  if (!dir.exists(analysis_dir)) {
    cat("SKIP: no existe carpeta de DE:", analysis_dir, "\n\n")
    next
  }

  all_files <- list.files(analysis_dir, pattern = "DE_.*_all_.*\\.tsv$", full.names = TRUE)
  if (length(all_files) == 0) {
    cat("SKIP: no se encontró tabla *_all_*.tsv en:", analysis_dir, "\n\n")
    next
  }
  all_fp <- all_files[order(file.info(all_files)$mtime, decreasing = TRUE)][1]
  cat("Usando DE table:", all_fp, "\n")
  de <- read.delim(all_fp, check.names = FALSE)

  if (!all(c("feature_id","FDR") %in% colnames(de))) {
    cat("SKIP: la tabla DE no tiene feature_id/FDR\n\n")
    next
  }

  # Derivar muestras y anotación
  der <- derive_samples_and_annotation(s, meta_base)
  samp_ids <- der$samp_ids
  ann <- der$ann

  # asegurar que existan en logcpm_base
  samp_ids <- intersect(samp_ids, colnames(logcpm_base))
  if (length(samp_ids) < 6) {
    cat("SKIP: muy pocas muestras tras intersección con logcpm_base: n=", length(samp_ids), "\n\n", sep="")
    next
  }
  ann <- ann[samp_ids, , drop = FALSE]

  # output dirs
  out_dir_analysis <- file.path(out_root, analysis_fs)
  dir.create(out_dir_analysis, showWarnings = FALSE, recursive = TRUE)

  for (thr in sort(fdrs)) {
    sig <- de[is.finite(as_num_safe(de$FDR)) & as_num_safe(de$FDR) < thr, , drop = FALSE]
    n_sig <- nrow(sig)
    cat(sprintf("FDR<%.2g: n_sig=%d\n", thr, n_sig))
    if (n_sig == 0) next

    sig$FDR_num <- as_num_safe(sig$FDR)
    sig <- sig[order(sig$FDR_num, na.last = TRUE), , drop = FALSE]
    if (n_sig > max_features) {
      sig <- sig[seq_len(max_features), , drop = FALSE]
      cat("  -> reducido a top ", max_features, " por FDR\n", sep="")
    }

    feat_list <- unique(as.character(sig$feature_id))
    feat_list <- feat_list[feat_list %in% rownames(logcpm_base)]
    if (length(feat_list) == 0) {
      cat("  -> 0 features presentes en logcpm_base\n")
      next
    }

    thr_tag <- paste0("FDR", gsub("\\.", "", sprintf("%.2g", thr)))
    out_dir_thr <- file.path(out_dir_analysis, thr_tag)
    dir.create(out_dir_thr, showWarnings = FALSE, recursive = TRUE)

    # PDF multipágina
    pdf_fp <- NULL
    if (make_pdf) {
      pdf_fp <- file.path(out_dir_thr, paste0("QCplots_", analysis_fs, "_", thr_tag, "_", ts, ".pdf"))
      pdf(pdf_fp, width = 10, height = 7)
    }

    # Para continuous: vector x (original) para el orden de muestras
    x_cont <- NULL
    if (mode == "continuous") {
      x_cont <- ann$value
      names(x_cont) <- rownames(ann)
    }

    # generar plots
    for (feat in feat_list) {
      y <- logcpm_base[feat, samp_ids]
      names(y) <- samp_ids

      if (mode == "continuous") {
        x_vals <- x_cont[samp_ids]
        ok_expr <- is.finite(y)
        ok_x <- is.finite(x_vals)
        ok_use <- ok_expr & ok_x
        n_start <- length(samp_ids)
        n_used <- sum(ok_use)
        n_excl <- n_start - n_used
        cat("[PLOT QC] ", analysis_id, " | ", feat,
            " | n_start=", n_start,
            " n_used=", n_used,
            " n_excluded=", n_excl,
            " (expr_NA=", sum(!ok_expr),
            " x_NA=", sum(ok_expr & !ok_x), ")\n", sep="")

        out_png <- file.path(out_dir_thr, paste0(safe_name(feat), "_scatter_", ts, ".png"))
        # scatter usa x = valor continuo (sin escalar; la escala es interpretativa)
        res <- plot_scatter_cont(
          x = x_cont[samp_ids],
          y = y,
          feat = feat,
          title = der$var,
          out_png = out_png,
          label_outliers = label_outliers,
          n_label = n_label
        )
        if (make_pdf) {
          # repetir en PDF
          plot_scatter_cont(
            x = x_cont[samp_ids],
            y = y,
            feat = feat,
            title = der$var,
            out_png = tempfile(fileext = ".png"),
            label_outliers = FALSE,
            n_label = n_label
          )
          # en PDF no hay que guardar PNG; recreamos el plot directamente
          # (replot: simple)
          ok <- is.finite(x_cont[samp_ids]) & is.finite(y)
          xx <- x_cont[samp_ids][ok]; yy <- y[ok]
          sp <- suppressWarnings(stats::cor.test(xx, yy, method = "spearman", exact = FALSE))
          plot(xx, yy, pch = 16, xlab = der$var, ylab = "logCPM (TMM)", main = feat)
          abline(stats::lm(yy ~ xx), lwd = 2)
          mtext(sprintf("Spearman rho=%.3f, p=%.3g, n=%d", unname(sp$estimate), sp$p.value, length(xx)),
                side = 3, line = 0.2, cex = 0.95)
        }

      } else {
        group <- ann$group
        group <- droplevels(group)

        title <- paste0(analysis_id, " | ", feat)
        subtitle <- paste0("(", levels(group)[2], " vs ", levels(group)[1], ")  [visual QC]")
        out_png <- file.path(out_dir_thr, paste0(safe_name(feat), "_boxjitter_", ts, ".png"))

        ok_expr <- is.finite(y)
        n_start <- length(samp_ids)
        n_used <- sum(ok_expr)
        n_excl <- n_start - n_used
        cat("[PLOT QC] ", analysis_id, " | ", feat,
            " | n_start=", n_start,
            " n_used=", n_used,
            " n_excluded=", n_excl,
            " (expr_NA=", sum(!ok_expr), ")\n", sep="")

        plot_box_jitter(y = y, group = group, title = title, subtitle = subtitle, out_png = out_png)

        if (make_pdf) {
          y_lim <- c(5, 15)
          boxplot(y ~ group, main = title, ylab = "logCPM (TMM)", xlab = "",
                  las = 2, outline = FALSE, ylim = y_lim)
          stripchart(y ~ group, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.8)
          mtext(subtitle, side = 3, line = 0.2, cex = 0.9)
        }
      }
    }

    if (make_pdf) {
      dev.off()
      cat("  -> PDF:", pdf_fp, "\n")
    }

    # guardar lista
    out_list <- file.path(out_dir_thr, paste0("features_", analysis_fs, "_", thr_tag, "_", ts, ".txt"))
    writeLines(feat_list, con = out_list)
    cat("  -> saved:", length(feat_list), "features plots in", out_dir_thr, "\n")
  }

  cat("\n")
}

cat("FIN. Log:", log_fp, "\n")
sink(type="message"); sink(type="output"); close(zz)
