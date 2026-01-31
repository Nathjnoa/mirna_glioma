#!/usr/bin/env Rscript
# ==========================================================================
# GSEA Enrichment Plot for SurvivalRank miRNA GSEA (miEAA results)
# - Uses Cox z scores as ranking metric
# - Builds classic running-sum enrichment plot for a depleted term
# ==========================================================================
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}
as_num_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ---- Inputs ----
in_root <- get_arg("--in_root", file.path("results", "tables", "SurvivalRank_GSEA"))
out_root <- get_arg("--out_root", file.path("results", "figures", "SurvivalRank_GSEA"))
run_tag <- get_arg("--run_tag", "SurvivalRank_CoxZ_miRPathDB")
term_query <- get_arg("--term", "")  # optional; match Subcategory
category_regex <- get_arg("--category_regex", "GO Biological process")

enrichment_dir <- tolower(get_arg("--enrichment", "depleted"))
p_adj_cutoff <- as_num_safe(get_arg("--p_adj_cutoff", "0.05"))
weight_p <- as_num_safe(get_arg("--weight_p", "1"))

preset <- get_arg("--preset", "single_col")
seed <- as.integer(get_arg("--seed", "42"))
stamp <- tolower(get_arg("--stamp", "true")) %in% c("true","t","1","yes","y")

if (!nzchar(run_tag)) stop("run_tag cannot be empty")
if (!is.finite(p_adj_cutoff) || p_adj_cutoff <= 0) stop("p_adj_cutoff invalido")
if (!is.finite(weight_p) || weight_p < 0) stop("weight_p invalido")

# ---- Resolve paths ----
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
proj_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
resolve_path <- function(p) {
  if (file.exists(p) || dir.exists(p)) return(p)
  alt <- file.path(proj_root, p)
  if (file.exists(alt) || dir.exists(alt)) return(alt)
  p
}

in_root <- resolve_path(in_root)
out_root <- resolve_path(out_root)

pick_latest <- function(dir_path, pattern) {
  if (!dir.exists(dir_path)) return(NA_character_)
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
  files[1]
}

run_dir <- file.path(in_root, run_tag)
cox_fp <- get_arg("--cox_table", pick_latest(run_dir, "^Cox_univariate_miRNA_all_.*\\.tsv$"))
mieaa_fp <- get_arg("--mieaa_table", pick_latest(run_dir, "^miEAA_GSEA_all_.*\\.tsv$"))

cox_fp <- resolve_path(cox_fp)
mieaa_fp <- resolve_path(mieaa_fp)

if (!file.exists(cox_fp)) stop("No encuentro Cox table: ", cox_fp)
if (!file.exists(mieaa_fp)) stop("No encuentro miEAA table: ", mieaa_fp)

# ---- Logging ----
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, run_tag), showWarnings = FALSE, recursive = TRUE)
log_fp <- file.path("logs", paste0("survGSEA_enrichment_plot_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("Enrichment plot for SurvivalRank GSEA\n")
cat("Cox table:", cox_fp, "\n")
cat("miEAA table:", mieaa_fp, "\n")
cat("run_tag:", run_tag, "\n")
cat("category_regex:", category_regex, "\n")
cat("enrichment:", enrichment_dir, "\n")
cat("p_adj_cutoff:", p_adj_cutoff, "\n")
cat("weight_p:", weight_p, "\n")
cat("term_query:", if (nzchar(term_query)) term_query else "(auto)", "\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(dplyr)
})

set.seed(seed)

# ---- Publication-ready presets ----
presets <- list(
  single_col = list(
    width_mm = 95, height_mm = 95, base = 8, title = 9, axis = 8,
    legend = 7, ticks = 7, line = 0.6, point = 1.6, spacing_mm = 2,
    margin_mm = 7, dpi_vector = 300, dpi_raster = 600
  ),
  double_col = list(
    width_mm = 190, height_mm = 140, base = 8.5, title = 10, axis = 8.5,
    legend = 7.5, ticks = 7.5, line = 0.7, point = 1.8, spacing_mm = 2.5,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  presentation = list(
    width_mm = 254, height_mm = 160, base = 18, title = 22, axis = 18,
    legend = 16, ticks = 16, line = 1.5, point = 4, spacing_mm = 6,
    margin_mm = 10, dpi_vector = 300, dpi_raster = 300
  ),
  poster = list(
    width_mm = 508, height_mm = 356, base = 28, title = 34, axis = 28,
    legend = 24, ticks = 24, line = 2.0, point = 6, spacing_mm = 10,
    margin_mm = 15, dpi_vector = 300, dpi_raster = 300
  )
)
if (!preset %in% names(presets)) stop("preset invalido: ", preset)
p <- presets[[preset]]

# ---- Theme ----
theme_pub <- function() {
  theme_minimal(base_size = p$base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = p$title, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = p$axis),
      axis.title = element_text(size = p$axis),
      axis.text = element_text(size = p$ticks, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.25),
      plot.margin = margin(p$margin_mm, p$margin_mm + 4, p$margin_mm, p$margin_mm, unit = "mm")
    )
}

# ---- Load data ----
cox_df <- read.table(cox_fp, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
if (!"feature_id" %in% names(cox_df)) stop("Cox table missing feature_id column")
if (!"rank_stat" %in% names(cox_df)) {
  if (!"z" %in% names(cox_df)) stop("Cox table missing rank_stat or z column")
  cox_df$rank_stat <- cox_df$z
}
cox_df$rank_stat <- as_num_safe(cox_df$rank_stat)
cox_df <- cox_df[is.finite(cox_df$rank_stat) & nzchar(cox_df$feature_id), ]
cox_df <- cox_df[order(cox_df$rank_stat, decreasing = TRUE), ]
cox_df <- cox_df[!duplicated(cox_df$feature_id), ]

rank_stats <- cox_df$rank_stat
names(rank_stats) <- cox_df$feature_id

mieaa <- read.table(mieaa_fp, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
# normalize colnames if needed
if (!"Subcategory" %in% names(mieaa) || !"Enrichment" %in% names(mieaa)) {
  stop("miEAA table missing required columns (Subcategory, Enrichment)")
}
if (!"P-adjusted" %in% names(mieaa)) stop("miEAA table missing P-adjusted")

mieaa$P.adjusted <- as_num_safe(mieaa$`P-adjusted`)

sel <- mieaa %>%
  filter(grepl(category_regex, Category, ignore.case = TRUE)) %>%
  filter(tolower(Enrichment) == enrichment_dir) %>%
  filter(is.finite(P.adjusted) & P.adjusted <= p_adj_cutoff)

if (nrow(sel) == 0) {
  cat("WARN: No rows for category_regex; falling back to any category.\n")
  sel <- mieaa %>%
    filter(tolower(Enrichment) == enrichment_dir) %>%
    filter(is.finite(P.adjusted) & P.adjusted <= p_adj_cutoff)
}

if (nrow(sel) == 0) stop("No significant terms found for given filters.")

if (nzchar(term_query)) {
  sel2 <- sel[grepl(term_query, sel$Subcategory, ignore.case = TRUE), , drop = FALSE]
  if (nrow(sel2) == 0) stop("No term matched term_query in Subcategory: ", term_query)
  sel <- sel2
}

# take most significant
sel <- sel[order(sel$P.adjusted, sel$`P-value`), ]
term <- sel[1, , drop = FALSE]
term_name <- as.character(term$Subcategory)
term_cat <- as.character(term$Category)
term_padj <- term$P.adjusted
term_enr <- as.character(term$Enrichment)

# parse miRNA list
if (!"miRNAs/precursors" %in% names(term)) stop("miEAA table missing miRNAs/precursors column")
mir_list <- unlist(strsplit(as.character(term$`miRNAs/precursors`), ";"))
mir_list <- trimws(mir_list)
mir_list <- mir_list[nzchar(mir_list)]
mir_list <- unique(mir_list)

in_set <- names(rank_stats) %in% mir_list
Nh <- sum(in_set)
N <- length(rank_stats)

if (Nh < 5) stop("Term has <5 miRNAs present in ranked list (Nh=", Nh, ")")
if (Nh >= N) stop("Term has all miRNAs in ranked list (Nh=", Nh, ")")

weights <- abs(rank_stats) ^ weight_p
sum_hit <- sum(weights[in_set])
if (sum_hit == 0) stop("Sum of hit weights is 0; check rank_stat values.")

inc <- ifelse(in_set, weights / sum_hit, -1 / (N - Nh))
running <- cumsum(inc)

es_max <- max(running)
es_min <- min(running)
if (abs(es_max) >= abs(es_min)) {
  es <- es_max
  es_pos <- which.max(running)
} else {
  es <- es_min
  es_pos <- which.min(running)
}

plot_df <- data.frame(
  pos = seq_len(N),
  running = running,
  rank_stat = rank_stats
)

hits_df <- data.frame(pos = which(in_set))

# ---- Plot ----
col_es <- "#0072B2"  # Okabe-Ito blue
col_metric <- "#4D4D4D"
col_accent <- "#D55E00"

subtitle_raw <- paste0(term_cat, " | ", term_enr, " | FDR=", signif(term_padj, 3), " | Nh=", Nh)
wrap_width <- if (preset == "single_col") 45 else if (preset == "double_col") 80 else 120
subtitle <- stringr::str_wrap(subtitle_raw, width = wrap_width)

p1 <- ggplot(plot_df, aes(x = pos, y = running)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35) +
  geom_line(color = col_es, linewidth = p$line) +
  geom_vline(xintercept = es_pos, color = col_accent, linetype = "dashed", linewidth = 0.4) +
  labs(title = term_name, subtitle = subtitle, y = "Running ES", x = NULL) +
  theme_pub()

p2 <- ggplot(hits_df, aes(x = pos)) +
  geom_segment(aes(xend = pos, y = 0, yend = 1), color = "black", linewidth = 0.35) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(y = NULL, x = NULL) +
  theme_pub() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(0, p$margin_mm, 0, p$margin_mm, unit = "mm")
  )

p3 <- ggplot(plot_df, aes(x = pos, y = rank_stat)) +
  geom_line(color = col_metric, linewidth = p$line) +
  labs(y = "Cox z", x = "Ranked miRNAs (desc z)") +
  theme_pub()

fig <- p1 / p2 / p3 + plot_layout(heights = c(2, 0.6, 1))

if (stamp) {
  fig <- fig + plot_annotation(
    caption = paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"))
  )
}

# ---- Export ----
out_dir <- file.path(out_root, run_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
base <- file.path(out_dir, paste0("GSEA_enrichmentplot_", safe_name(term_name), "_", ts))

pdf_fp <- paste0(base, ".pdf")
svg_fp <- paste0(base, ".svg")
png_fp <- paste0(base, ".png")

# use cairo_pdf if available
if (capabilities("cairo")) {
  ggsave(pdf_fp, fig, width = p$width_mm, height = p$height_mm, units = "mm", device = cairo_pdf, limitsize = FALSE)
} else {
  ggsave(pdf_fp, fig, width = p$width_mm, height = p$height_mm, units = "mm", limitsize = FALSE)
}

ggsave(svg_fp, fig, width = p$width_mm, height = p$height_mm, units = "mm", limitsize = FALSE)
ggsave(png_fp, fig, width = p$width_mm, height = p$height_mm, units = "mm",
       dpi = p$dpi_raster, limitsize = FALSE)

cat("Saved:\n", pdf_fp, "\n", svg_fp, "\n", png_fp, "\n")
cat("Selected term:", term_name, "\n")
cat("Category:", term_cat, "\n")
cat("Enrichment:", term_enr, "\n")
cat("FDR:", term_padj, "\n")
cat("Nh:", Nh, " N:", N, "\n")
cat("ES:", es, " at pos ", es_pos, "\n", sep = "")
