#!/usr/bin/env Rscript
# ============================================================================
# Pathway redundancy reduction using Jaccard similarity + hierarchical clustering
#
# METHOD:
# Based on Enrichment Map approach (Merico et al. 2010, PLoS ONE)
# - Jaccard coefficient for gene set overlap
# - Hierarchical clustering (UPGMA, average linkage)
# - Dynamic tree cutting with similarity threshold
# - Representative selection by minimum Q-value
#
# SCOPE:
# - GO BP/MF: rrvgo (semantic similarity) - keeps existing approach
# - KEGG/Reactome: Jaccard clustering - NEW robust approach
#
# REFERENCES:
# - Merico D et al. (2010) Enrichment map: a network-based method for
#   gene-set enrichment visualization and interpretation. PLoS ONE 5(11):e13984
# - Yu G et al. (2012) clusterProfiler: an R package for comparing biological
#   themes among gene clusters. OMICS 16(5):284-287
#
# ============================================================================
options(stringsAsFactors = FALSE)

# ---- CLI argument parsing ----
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

# ---- CLI defaults ----
in_root   <- get_arg("--in_root", file.path("results_29s", "tables", "SurvivalRank_GSEA"))
out_root  <- get_arg("--out_root", file.path("results_29s", "figures", "SurvivalRank_GSEA"))
tbl_root  <- get_arg("--tbl_root", file.path("results_29s", "tables", "SurvivalRank_GSEA"))
run_tag   <- get_arg("--run_tag", "SurvivalRank_CoxZ_miRPathDB")
plot_label <- get_arg("--label", run_tag)

# Clustering parameters - database-specific (adaptive thresholds)
# KEGG: Lower threshold due to smaller pathway count and sparser overlap
jaccard_threshold_kegg     <- as_num_safe(get_arg("--jaccard_threshold_kegg", "0.25"))
min_cluster_size_kegg      <- as.integer(get_arg("--min_cluster_size_kegg", "1"))

# Reactome: Standard threshold (Merico et al. default ~0.375-0.5)
jaccard_threshold_reactome <- as_num_safe(get_arg("--jaccard_threshold_reactome", "0.5"))
min_cluster_size_reactome  <- as.integer(get_arg("--min_cluster_size_reactome", "2"))

linkage_method    <- get_arg("--linkage_method", "average")  # UPGMA

# Visualization
n_enriched <- as.integer(get_arg("--n_enriched", "12"))
n_depleted <- as.integer(get_arg("--n_depleted", "12"))
wrap_width <- as.integer(get_arg("--wrap_width", "58"))
preset     <- get_arg("--preset", "double_col")
stamp      <- tolower(get_arg("--stamp", "true")) %in% c("true","t","1","yes","y")
seed       <- as.integer(get_arg("--seed", "42"))

# Use rrvgo for GO terms?
use_rrvgo_for_go <- tolower(get_arg("--use_rrvgo", "true")) %in% c("true","t","1","yes","y")
rrvgo_threshold  <- as_num_safe(get_arg("--rrvgo_threshold", "0.9"))
rrvgo_method     <- get_arg("--rrvgo_method", "Rel")

# ---- Resolve paths ----
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
proj_root  <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
resolve_path <- function(p) {
  if (file.exists(p) || dir.exists(p)) return(p)
  alt <- file.path(proj_root, p)
  if (file.exists(alt) || dir.exists(alt)) return(alt)
  p
}
in_root  <- resolve_path(in_root)
out_root <- resolve_path(out_root)
tbl_root <- resolve_path(tbl_root)

# ---- Logging ----
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
out_tbl_dir <- file.path(tbl_root, run_tag, "reduced_jaccard")
out_fig_dir <- file.path(out_root, run_tag, "reduced_jaccard")
dir.create(out_tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

log_fp <- file.path("logs", paste0("reduce_redundancy_jaccard_", ts, ".txt"))
zz <- file(log_fp, open = "wt")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(zz) }, add = TRUE)

cat("============================================================================\n")
cat("Pathway Redundancy Reduction - Jaccard Clustering (Enrichment Map Method)\n")
cat("============================================================================\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Input root:", in_root, "\n")
cat("Table output:", out_tbl_dir, "\n")
cat("Figure output:", out_fig_dir, "\n")
cat("Run tag:", run_tag, "\n")
cat("KEGG settings:\n")
cat("  Jaccard threshold:", jaccard_threshold_kegg, "\n")
cat("  Min cluster size:", min_cluster_size_kegg, "\n")
cat("Reactome settings:\n")
cat("  Jaccard threshold:", jaccard_threshold_reactome, "\n")
cat("  Min cluster size:", min_cluster_size_reactome, "\n")
cat("Linkage method:", linkage_method, "(UPGMA)\n")
cat("Use rrvgo for GO:", use_rrvgo_for_go, "\n")
if (use_rrvgo_for_go) {
  cat("  rrvgo threshold:", rrvgo_threshold, "\n")
  cat("  rrvgo method:", rrvgo_method, "\n")
}
cat("============================================================================\n\n")

cat("REFERENCE:\n")
cat("  Merico D et al. (2010) Enrichment map: a network-based method for\n")
cat("  gene-set enrichment visualization and interpretation.\n")
cat("  PLoS ONE 5(11):e13984. doi: 10.1371/journal.pone.0013984\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

if (use_rrvgo_for_go) {
  suppressPackageStartupMessages({
    library(GO.db)
    library(rrvgo)
  })
}

set.seed(seed)

# ============================================================================
# PART 1: Load miEAA results
# ============================================================================
cat("=== PART 1: Loading miEAA GSEA results ===\n")

find_latest <- function(root, run_tag, pattern) {
  dir_path <- file.path(root, run_tag)
  if (!dir.exists(dir_path)) return(NA_character_)
  ff <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(ff) == 0) return(NA_character_)
  ff[order(file.info(ff)$mtime, decreasing = TRUE)][1]
}

use_fp <- find_latest(in_root, run_tag, "^miEAA_GSEA_all_miRPathDB_expert_.*\\.tsv$")
if (is.na(use_fp) || !file.exists(use_fp)) stop("No miEAA results found in: ", file.path(in_root, run_tag))
cat("Input file:", use_fp, "\n")

df_raw <- read.delim(use_fp, check.names = FALSE)
names(df_raw) <- trimws(names(df_raw))
if (nrow(df_raw) == 0) stop("Empty miEAA table.")

# Standardize columns
df <- df_raw %>%
  mutate(
    pval      = as_num_safe(`P-value`),
    padj      = as_num_safe(`P-adjusted`),
    qval      = as_num_safe(`Q-value`),
    observed  = as_num_safe(Observed),
    term      = str_squish(as.character(Subcategory)),
    category  = as.character(Category),
    enrichment = tolower(as.character(Enrichment)),
    mirnas_raw = as.character(`miRNAs/precursors`)
  )

cat("Total terms loaded:", nrow(df), "\n")
cat("Enriched:", sum(df$enrichment == "enriched"),
    " | Depleted:", sum(df$enrichment == "depleted"), "\n")

# Assign database labels
df$database <- case_when(
  grepl("GO.*Biological.process", df$category, ignore.case = TRUE) ~ "GO_BP",
  grepl("GO.*Molecular.function", df$category, ignore.case = TRUE) ~ "GO_MF",
  grepl("KEGG", df$category, ignore.case = TRUE) ~ "KEGG",
  grepl("Reactome", df$category, ignore.case = TRUE) ~ "Reactome",
  TRUE ~ "Other"
)

cat("\nTerms per database:\n")
print(table(df$database))
cat("\n")

# ============================================================================
# PART 2: Parse miRNA lists
# ============================================================================
cat("=== PART 2: Parsing miRNA-to-pathway mappings ===\n")

parse_mirnas <- function(mirna_str) {
  if (is.na(mirna_str) || nzchar(mirna_str) == 0) return(character(0))
  mirnas <- unlist(strsplit(as.character(mirna_str), ";", fixed = TRUE))
  mirnas <- trimws(mirnas)
  mirnas <- mirnas[nzchar(mirnas)]
  unique(mirnas)
}

df$mirnas_list <- lapply(df$mirnas_raw, parse_mirnas)
df$n_mirnas <- sapply(df$mirnas_list, length)

cat("miRNA parsing summary:\n")
cat("  Total unique terms:", nrow(df), "\n")
cat("  Terms with ≥1 miRNA:", sum(df$n_mirnas > 0), "\n")
cat("  Median miRNAs per term:", median(df$n_mirnas), "\n")
cat("  Range:", min(df$n_mirnas), "-", max(df$n_mirnas), "\n\n")

# ============================================================================
# PART 3: Jaccard similarity function
# ============================================================================
cat("=== PART 3: Defining Jaccard similarity function ===\n")

#' Calculate Jaccard similarity matrix for gene sets
#' Following Merico et al. (2010) Enrichment Map implementation
calc_jaccard_matrix <- function(gene_sets) {
  n <- length(gene_sets)
  jaccard_mat <- matrix(0, nrow = n, ncol = n,
                        dimnames = list(names(gene_sets), names(gene_sets)))

  for (i in seq_len(n)) {
    jaccard_mat[i, i] <- 1.0
    if (i < n) {
      for (j in (i+1):n) {
        set1 <- gene_sets[[i]]
        set2 <- gene_sets[[j]]

        intersection <- length(intersect(set1, set2))
        union <- length(union(set1, set2))

        jacc <- if (union == 0) 0 else intersection / union

        jaccard_mat[i, j] <- jacc
        jaccard_mat[j, i] <- jacc
      }
    }
  }

  jaccard_mat
}

cat("  Jaccard coefficient: J(A,B) = |A ∩ B| / |A ∪ B|\n")
cat("  Range: [0, 1], where 1 = identical gene sets\n\n")

# ============================================================================
# PART 4: Hierarchical clustering function
# ============================================================================
cat("=== PART 4: Defining hierarchical clustering function ===\n")

#' Cluster pathways by Jaccard similarity
#' @param jaccard_mat Jaccard similarity matrix
#' @param threshold Jaccard threshold for clustering (default: 0.5)
#' @param method Linkage method (default: "average" = UPGMA)
cluster_by_jaccard <- function(jaccard_mat, threshold = 0.5, method = "average") {
  if (nrow(jaccard_mat) < 2) {
    cat("    Too few pathways (<2), no clustering performed.\n")
    return(list(clusters = rep(1L, nrow(jaccard_mat)), n_clusters = 1))
  }

  # Convert similarity to distance
  dist_mat <- as.dist(1 - jaccard_mat)

  # Hierarchical clustering
  hc <- hclust(dist_mat, method = method)

  # Cut tree at height = 1 - threshold
  # (height in distance space corresponds to 1 - similarity)
  clusters <- cutree(hc, h = 1 - threshold)

  list(
    clusters = clusters,
    n_clusters = length(unique(clusters)),
    hclust = hc,
    threshold = threshold
  )
}

cat("  Linkage method:", linkage_method, "(UPGMA)\n")
cat("  Distance metric: D = 1 - Jaccard\n")
cat("  Threshold varies by database (see settings above)\n\n")

# ============================================================================
# PART 5: Apply clustering to KEGG/Reactome
# ============================================================================
cat("=== PART 5: Applying Jaccard clustering ===\n")

# Storage for cluster assignments
df$cluster <- NA_integer_
df$cluster_rep <- FALSE
df$reduction_method <- NA_character_

results_summary <- list()

# Function to process one database × direction combination
process_pathways <- function(df_subset, db_name, dir_label, jaccard_thresh, min_clust_size) {
  cat(sprintf("\n--- %s %s (n=%d) ---\n", db_name, dir_label, nrow(df_subset)))
  cat(sprintf("    Using Jaccard threshold: %.2f, Min cluster size: %d\n",
              jaccard_thresh, min_clust_size))

  if (nrow(df_subset) < 2) {
    cat("  Skipping: too few pathways (<2)\n")
    return(NULL)
  }

  # Build gene sets (miRNA sets in this case)
  gene_sets <- df_subset$mirnas_list
  names(gene_sets) <- df_subset$term

  # Calculate Jaccard matrix
  cat("  Computing Jaccard similarity matrix...\n")
  jaccard_mat <- calc_jaccard_matrix(gene_sets)

  median_sim <- median(jaccard_mat[lower.tri(jaccard_mat)])
  cat(sprintf("    Median Jaccard similarity: %.3f\n", median_sim))

  # Cluster
  cat("  Performing hierarchical clustering...\n")
  clust_result <- cluster_by_jaccard(jaccard_mat,
                                     threshold = jaccard_thresh,
                                     method = linkage_method)

  cat(sprintf("    Number of clusters: %d\n", clust_result$n_clusters))

  # Filter small clusters
  cluster_sizes <- table(clust_result$clusters)
  valid_clusters <- as.integer(names(cluster_sizes[cluster_sizes >= min_clust_size]))

  clusters_filtered <- ifelse(clust_result$clusters %in% valid_clusters,
                               clust_result$clusters,
                               NA_integer_)

  n_removed <- sum(is.na(clusters_filtered))
  cat(sprintf("    Clusters after size filter (≥%d): %d (removed %d singletons)\n",
              min_clust_size,
              length(unique(clusters_filtered[!is.na(clusters_filtered)])),
              n_removed))

  # Select representatives (lowest Q-value per cluster)
  df_subset$cluster_tmp <- clusters_filtered

  representatives <- df_subset %>%
    filter(!is.na(cluster_tmp)) %>%
    group_by(cluster_tmp) %>%
    arrange(qval, padj) %>%
    slice_head(n = 1) %>%
    pull(term)

  cat(sprintf("    Representative pathways: %d (%.1f%% reduction)\n",
              length(representatives),
              100 * (1 - length(representatives) / nrow(df_subset))))

  # Update global df
  idx <- match(df_subset$term, df$term)
  df$cluster[idx] <<- clusters_filtered
  df$cluster_rep[idx] <<- df_subset$term %in% representatives
  df$reduction_method[idx] <<- sprintf("Jaccard_hclust_%s_t%.2f", linkage_method, jaccard_thresh)

  # Return summary
  list(
    database = db_name,
    direction = dir_label,
    method = "Jaccard_hclust",
    n_input = nrow(df_subset),
    n_output = length(representatives),
    reduction_pct = round(100 * (1 - length(representatives) / nrow(df_subset)), 1),
    n_clusters = length(unique(clusters_filtered[!is.na(clusters_filtered)])),
    median_jaccard = median_sim,
    threshold = jaccard_thresh,
    min_cluster_size = min_clust_size
  )
}

# Process KEGG (with KEGG-specific parameters)
for (dir_label in c("enriched", "depleted")) {
  df_sub <- df %>% filter(database == "KEGG", enrichment == dir_label, n_mirnas > 0)

  if (nrow(df_sub) >= 2) {
    res <- process_pathways(df_sub, "KEGG", dir_label,
                           jaccard_thresh = jaccard_threshold_kegg,
                           min_clust_size = min_cluster_size_kegg)
    if (!is.null(res)) {
      results_summary[[paste0("KEGG_", dir_label)]] <- res
    }
  } else {
    cat(sprintf("\nSkipping KEGG %s: too few pathways (n=%d)\n", dir_label, nrow(df_sub)))
  }
}

# Process Reactome (with Reactome-specific parameters)
for (dir_label in c("enriched", "depleted")) {
  df_sub <- df %>% filter(database == "Reactome", enrichment == dir_label, n_mirnas > 0)

  if (nrow(df_sub) >= 2) {
    res <- process_pathways(df_sub, "Reactome", dir_label,
                           jaccard_thresh = jaccard_threshold_reactome,
                           min_clust_size = min_cluster_size_reactome)
    if (!is.null(res)) {
      results_summary[[paste0("Reactome_", dir_label)]] <- res
    }
  } else {
    cat(sprintf("\nSkipping Reactome %s: too few pathways (n=%d)\n", dir_label, nrow(df_sub)))
  }
}

# ============================================================================
# PART 6: Apply rrvgo to GO terms (optional, keeps existing method)
# ============================================================================
if (use_rrvgo_for_go) {
  cat("\n=== PART 6: Applying rrvgo to GO terms ===\n")
  cat("(Using existing rrvgo method for GO BP only - semantic similarity)\n")
  cat("(GO MF is SKIPPED - will not appear in final results)\n\n")

  # Map GO term names to IDs
  go_terms_db <- AnnotationDbi::select(GO.db,
                                       keys = keys(GO.db, "GOID"),
                                       columns = c("GOID", "TERM", "ONTOLOGY"),
                                       keytype = "GOID")
  go_terms_db$term_norm <- tolower(str_squish(go_terms_db$TERM))

  df$go_id <- NA_character_
  is_go <- df$database %in% c("GO_BP")  # Only GO_BP, skip GO_MF

  if (any(is_go)) {
    terms_to_map <- tolower(str_squish(df$term[is_go]))
    match_idx <- match(terms_to_map, go_terms_db$term_norm)

    mapped <- !is.na(match_idx)
    df$go_id[is_go][mapped] <- go_terms_db$GOID[match_idx[mapped]]

    cat("GO terms mapped:", sum(mapped), "/", sum(is_go), "\n")

    # Apply rrvgo per ontology × direction (only BP, skip MF)
    for (ont in c("BP")) {
      ont_db <- "GO_BP"

      for (dir_label in c("enriched", "depleted")) {
        df_go <- df %>% filter(database == ont_db, enrichment == dir_label, !is.na(go_id))

        if (nrow(df_go) < 3) {
          cat(sprintf("  Skip GO %s %s (n=%d, need ≥3)\n", ont, dir_label, nrow(df_go)))
          next
        }

        cat(sprintf("  Processing GO %s %s (n=%d)...\n", ont, dir_label, nrow(df_go)))

        # Calculate semantic similarity
        go_ids <- unique(df_go$go_id)
        scores <- setNames(-log10(pmax(df_go$padj, .Machine$double.xmin)), df_go$go_id)

        sim_matrix <- tryCatch(
          calculateSimMatrix(go_ids,
                             orgdb = "org.Hs.eg.db",
                             ont = ont,
                             method = rrvgo_method),
          error = function(e) NULL
        )

        if (is.null(sim_matrix)) {
          cat("    ERROR: similarity matrix calculation failed\n")
          next
        }

        reduced <- tryCatch(
          reduceSimMatrix(sim_matrix,
                          scores = scores[rownames(sim_matrix)],
                          threshold = rrvgo_threshold,
                          orgdb = "org.Hs.eg.db"),
          error = function(e) NULL
        )

        if (is.null(reduced)) {
          cat("    ERROR: reduction failed\n")
          next
        }

        # Mark representatives
        parent_terms <- unique(reduced$parent)
        go_to_term <- setNames(df_go$term, df_go$go_id)
        rep_terms <- go_to_term[parent_terms]

        idx_go <- match(df_go$term, df$term)
        df$cluster_rep[idx_go] <- df_go$term %in% rep_terms
        df$reduction_method[idx_go] <- sprintf("rrvgo_%s_t%.2f", rrvgo_method, rrvgo_threshold)

        cat(sprintf("    Reduced: %d → %d (%.1f%% reduction)\n",
                    nrow(df_go), length(rep_terms),
                    100 * (1 - length(rep_terms) / nrow(df_go))))

        # Add to summary
        results_summary[[paste0("GO_", ont, "_", dir_label)]] <- list(
          database = paste0("GO_", ont),
          direction = dir_label,
          method = "rrvgo_semantic",
          n_input = nrow(df_go),
          n_output = length(rep_terms),
          reduction_pct = round(100 * (1 - length(rep_terms) / nrow(df_go)), 1),
          n_clusters = NA,
          median_jaccard = NA,
          threshold = rrvgo_threshold
        )
      }
    }
  }
} else {
  cat("\n=== PART 6: Skipping GO terms (use_rrvgo = FALSE) ===\n")
  cat("GO terms will pass through without reduction.\n")
  cat("Set --use_rrvgo true to apply rrvgo to GO BP/MF.\n\n")
}

# ============================================================================
# PART 7: Build reduced table
# ============================================================================
cat("\n=== PART 7: Building reduced table ===\n")

# For terms not processed, mark as representatives
df$cluster_rep[is.na(df$cluster_rep)] <- TRUE

# Exclude GO_MF entirely from final results
# Keep both cluster representatives AND unique pathways (singletons with cluster = NA)
df_reduced <- df %>%
  filter(database != "GO_MF") %>%
  filter(is.na(cluster) | cluster_rep)

cat("Before reduction:", nrow(df), "\n")
cat("After reduction:", nrow(df_reduced), "\n")
cat("Removed:", nrow(df) - nrow(df_reduced), "redundant pathways\n\n")

# Summary by database × direction
summary_tbl <- df_reduced %>%
  group_by(database, enrichment) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(database, enrichment)

cat("Pathways after reduction:\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\n")

# ============================================================================
# PART 8: Export tables
# ============================================================================
cat("=== PART 8: Exporting tables ===\n")

# Reduced table
export_cols <- c(
  intersect(names(df_raw), names(df_reduced)),
  "database", "cluster", "cluster_rep", "reduction_method"
)
export_cols <- intersect(export_cols, names(df_reduced))

# Convert list columns to character strings for export
df_reduced_export <- df_reduced[, export_cols]
list_cols <- sapply(df_reduced_export, function(col) inherits(col, "list") || is.list(col))
if (any(list_cols)) {
  for (col_name in names(df_reduced_export)[list_cols]) {
    df_reduced_export[[col_name]] <- sapply(df_reduced_export[[col_name]], function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })
  }
}

reduced_fp <- file.path(out_tbl_dir, paste0("miEAA_GSEA_reduced_jaccard_", ts, ".tsv"))
write.table(df_reduced_export, reduced_fp,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Reduced table:", reduced_fp, "\n")

# Full annotated table
df_full <- df
# Convert all list-type columns to character strings for export
list_cols <- sapply(df_full, function(col) inherits(col, "list") || is.list(col))
if (any(list_cols)) {
  for (col_name in names(df_full)[list_cols]) {
    df_full[[col_name]] <- sapply(df_full[[col_name]], function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })
  }
}
full_fp <- file.path(out_tbl_dir, paste0("miEAA_GSEA_full_annotated_jaccard_", ts, ".tsv"))
write.table(df_full, full_fp, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Full annotated table:", full_fp, "\n")

# Summary
summary_df <- do.call(rbind, lapply(results_summary, function(x) {
  data.frame(
    Database = x$database,
    Direction = x$direction,
    Method = x$method,
    Input = x$n_input,
    Output = x$n_output,
    Reduction_pct = x$reduction_pct,
    N_clusters = ifelse(is.na(x$n_clusters), "-", as.character(x$n_clusters)),
    Median_Jaccard = ifelse(is.na(x$median_jaccard), "-", sprintf("%.3f", x$median_jaccard)),
    Threshold = x$threshold,
    Min_cluster_size = ifelse(is.null(x$min_cluster_size), "-", as.character(x$min_cluster_size)),
    stringsAsFactors = FALSE
  )
}))

if (!is.null(summary_df) && nrow(summary_df) > 0) {
  summary_fp <- file.path(out_tbl_dir, paste0("reduction_summary_", ts, ".tsv"))
  write.table(summary_df, summary_fp, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Summary table:", summary_fp, "\n")

  cat("\nReduction summary:\n")
  print(summary_df, row.names = FALSE)
}

# Markdown summary
md_fp <- file.path(out_tbl_dir, paste0("reduction_summary_", ts, ".md"))
md_lines <- c(
  "# Pathway Redundancy Reduction Summary",
  "",
  paste0("**Method**: Jaccard similarity + hierarchical clustering (Enrichment Map)"),
  paste0("**Reference**: Merico D et al. (2010) PLoS ONE 5(11):e13984"),
  "",
  paste0("- **Run tag**: ", plot_label),
  paste0("- **Linkage method**: ", linkage_method, " (UPGMA)"),
  paste0("- **Date**: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Database-specific parameters:",
  "",
  "**KEGG pathways**:",
  paste0("  - Jaccard threshold: ", jaccard_threshold_kegg, " (permissive, due to sparse overlap)"),
  paste0("  - Min cluster size: ", min_cluster_size_kegg, " (retain singletons)"),
  "",
  "**Reactome pathways**:",
  paste0("  - Jaccard threshold: ", jaccard_threshold_reactome, " (standard)"),
  paste0("  - Min cluster size: ", min_cluster_size_reactome),
  "",
  "**Rationale**: KEGG pathways exhibit lower inherent gene overlap and smaller",
  "pathway count compared to Reactome, warranting a more permissive threshold to",
  "avoid over-reduction while maintaining biological interpretability.",
  "",
  "## Overall reduction:",
  paste0("- **Before reduction**: ", nrow(df)),
  paste0("- **After reduction**: ", nrow(df_reduced)),
  paste0("- **Removed**: ", nrow(df) - nrow(df_reduced), " redundant pathways (",
         round(100 * (nrow(df) - nrow(df_reduced)) / nrow(df), 1), "%)"),
  "",
  "## Reduction per database/direction",
  ""
)

if (!is.null(summary_df) && nrow(summary_df) > 0) {
  md_lines <- c(md_lines,
                "| Database | Direction | Method | Input | Output | Reduction % | Clusters | Median Jaccard | Threshold | Min Size |",
                "|----------|-----------|--------|-------|--------|-------------|----------|----------------|-----------|----------|")

  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    md_lines <- c(md_lines, sprintf(
      "| %s | %s | %s | %d | %d | %.1f%% | %s | %s | %s | %s |",
      row$Database, row$Direction, row$Method,
      as.integer(row$Input), as.integer(row$Output), as.numeric(row$Reduction_pct),
      as.character(row$N_clusters), as.character(row$Median_Jaccard),
      as.character(row$Threshold), as.character(row$Min_cluster_size)
    ))
  }
}

md_lines <- c(md_lines, "",
              "## Pathways per category (after reduction)",
              "",
              "| Database | Enriched | Depleted |",
              "|----|----:|----:|")

summary_wide <- summary_tbl %>%
  pivot_wider(names_from = enrichment, values_from = n, values_fill = 0)

for (i in seq_len(nrow(summary_wide))) {
  row <- summary_wide[i, ]
  md_lines <- c(md_lines, sprintf("| %s | %d | %d |",
                                   row$database,
                                   if ("enriched" %in% names(row)) row$enriched else 0,
                                   if ("depleted" %in% names(row)) row$depleted else 0))
}

writeLines(md_lines, md_fp)
cat("Summary markdown:", md_fp, "\n\n")

# ============================================================================
# PART 9: Bubble plot visualizations
# ============================================================================
cat("\n=== PART 9: Generating bubble plots ===\n")

# ---- Presets ----
presets <- list(
  single_col = list(
    width_mm = 85, height_mm = 65, base = 8, title = 9, axis = 8,
    legend = 7, ticks = 7, line = 0.6, point = 1.6, spacing_mm = 2,
    margin_mm = 6, dpi_vector = 300, dpi_raster = 600
  ),
  double_col = list(
    width_mm = 220, height_mm = 160, base = 30, title = 35, axis = 30,
    legend = 28, ticks = 28, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  presentation = list(
    width_mm = 254, height_mm = 143, base = 18, title = 22, axis = 18,
    legend = 16, ticks = 16, line = 1.5, point = 4, spacing_mm = 6,
    margin_mm = 10, dpi_vector = 300, dpi_raster = 300
  ),
  poster = list(
    width_mm = 508, height_mm = 356, base = 28, title = 34, axis = 28,
    legend = 24, ticks = 24, line = 2.0, point = 6, spacing_mm = 10,
    margin_mm = 15, dpi_vector = 300, dpi_raster = 300
  ),
  per_database = list(
    width_mm = 220, height_mm = 160, base = 22, title = 26, axis = 22,
    legend = 20, ticks = 20, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  enriched_depleted = list(
    width_mm = 220, height_mm = 160, base = 20, title = 24, axis = 20,
    legend = 18, ticks = 18, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  ),
  multipanel = list(
    width_mm = 220, height_mm = 160, base = 16, title = 20, axis = 16,
    legend = 14, ticks = 14, line = 0.8, point = 2.2, spacing_mm = 3,
    margin_mm = 8, dpi_vector = 300, dpi_raster = 600
  )
)
if (!preset %in% names(presets)) stop("Invalid preset: ", preset)
p <- presets[[preset]]

# ---- Theme ----
theme_pub <- function(preset_params = p) {
  theme_minimal(base_size = preset_params$base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = preset_params$title, face = "bold", hjust = 0),
      axis.title = element_text(size = preset_params$axis, face = "plain"),
      axis.text.x = element_text(size = preset_params$ticks, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = preset_params$axis, color = "black", lineheight = 1.1,
                                  margin = margin(r = 3, unit = "mm")),
      strip.text = element_text(size = preset_params$axis, face = "bold"),
      legend.title = element_text(size = preset_params$legend, face = "bold"),
      legend.text = element_text(size = preset_params$legend),
      legend.position = "right",
      legend.key.size = unit(4, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.border = element_blank(),
      panel.spacing.x = unit(preset_params$spacing_mm, "mm"),
      panel.spacing.y = unit(preset_params$spacing_mm * 2.4, "mm"),
      plot.margin = margin(preset_params$margin_mm, preset_params$margin_mm, preset_params$margin_mm, preset_params$margin_mm + 6, unit = "mm"),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

colors_direction <- c("Enriched" = "#009E73", "Depleted" = "#D55E00")

# ---- Helper functions ----
map_db <- function(x) {
  x <- tolower(x)
  db <- rep(NA_character_, length(x))
  db[grepl("go.*biological.process", x)] <- "GO Biological process (miRPathDB)"
  db[grepl("go.*molecular.function", x)] <- "GO Molecular function (miRPathDB)"
  db[grepl("kegg", x)] <- "KEGG (miRPathDB)"
  db[grepl("reactome", x)] <- "Reactome (miRPathDB)"
  db
}

save_plot <- function(plot, out_base, n_rows = 1, width_mult = 1, height_mult = 1) {
  width_mm  <- p$width_mm * width_mult
  height_mm <- p$height_mm * height_mult * n_rows
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
  pdf_path <- paste0(out_base, ".pdf")
  tryCatch({
    ggsave(pdf_path, plot, width = width_mm, height = height_mm,
           units = "mm", dpi = p$dpi_vector, limitsize = FALSE, device = pdf_dev)
    cat("  Exported:", pdf_path, sprintf("(%d x %d mm)\n", round(width_mm), round(height_mm)))
  }, error = function(e) warning("Failed PDF: ", conditionMessage(e)))
  svg_path <- paste0(out_base, ".svg")
  if (requireNamespace("svglite", quietly = TRUE)) {
    tryCatch({
      ggsave(svg_path, plot, width = width_mm, height = height_mm,
             units = "mm", dpi = p$dpi_vector, limitsize = FALSE, device = svglite::svglite)
      cat("  Exported:", svg_path, "\n")
    }, error = function(e) warning("Failed SVG: ", conditionMessage(e)))
  }
  png_path <- paste0(out_base, ".png")
  tryCatch({
    ggsave(png_path, plot, width = width_mm, height = height_mm,
           units = "mm", dpi = p$dpi_raster, limitsize = FALSE, device = "png")
    cat("  Exported:", png_path, sprintf("(%d dpi)\n", p$dpi_raster))
  }, error = function(e) warning("Failed PNG: ", conditionMessage(e)))
  invisible(NULL)
}

# ---- Prepare plot data ----
df_plot <- df_reduced  # Use the already-created df_reduced from Part 7
df_plot$db <- map_db(df_plot$Category)  # Note: use "Category" not "category"
df_plot$direction <- ifelse(df_plot$Enrichment == "enriched", "Enriched", "Depleted")
df_plot$direction <- factor(df_plot$direction, levels = c("Enriched", "Depleted"))
df_plot$logp <- -log10(pmax(df_plot$`P-adjusted`, .Machine$double.xmin))

# Filter to 3 databases (exclude GO MF per script 16 pattern)
mirpath_levels <- c(
  "GO Biological process (miRPathDB)",
  "KEGG (miRPathDB)",
  "Reactome (miRPathDB)"
)
df_plot <- df_plot[df_plot$db %in% mirpath_levels &
                   is.finite(df_plot$`P-adjusted`) &
                   nzchar(df_plot$Subcategory), ]
df_plot$db <- factor(df_plot$db, levels = mirpath_levels)
df_plot$term <- str_squish(as.character(df_plot$Subcategory))

cat("Plot data prepared:", nrow(df_plot), "pathways (3 databases, excl. GO MF)\n")

# ---- Create output subdirectories ----
dir.create(file.path(out_fig_dir, "enriched_depleted"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_fig_dir, "per_database"), showWarnings = FALSE, recursive = TRUE)

# ---- Generate plots ----
plots_list <- list()
plots_per_db <- list()

for (dir_type in c("Enriched", "Depleted")) {
  n_terms <- if (dir_type == "Depleted") n_depleted else n_enriched

  df_dir <- df_plot %>%
    filter(direction == dir_type) %>%
    group_by(db) %>%
    arrange(`P-adjusted`, .by_group = TRUE) %>%  # Changed from padj
    slice_head(n = n_terms) %>%
    ungroup()

  dir_counts <- table(df_dir$db)
  cat(sprintf("  %s terms (top %d): %s\n", dir_type, n_terms,
              paste(names(dir_counts), dir_counts, sep = "=", collapse = ", ")))

  if (nrow(df_dir) == 0) {
    cat(sprintf("  SKIP %s: no terms available\n", dir_type))
    next
  }

  # Order terms within each DB
  df_dir <- df_dir %>%
    arrange(db, `P-adjusted`) %>%  # Changed from padj
    mutate(term_plot = paste(term, db, sep = "|||"))
  df_dir$term_plot <- factor(df_dir$term_plot, levels = rev(unique(df_dir$term_plot)))
  df_dir$observed <- as.numeric(df_dir$Observed)  # Ensure numeric for size aesthetic

  # Use enriched_depleted preset for multipanel enriched/depleted plots
  p_ed <- presets[["enriched_depleted"]]
  p <- p_ed  # Set p for save_plot function
  size_range_ed <- c(p_ed$point * 0.8, p_ed$point * 2.8)
  dir_color <- colors_direction[[dir_type]]

  # Multipanel plot (all databases)
  p_dir <- ggplot(df_dir, aes(x = logp, y = term_plot)) +
    geom_point(aes(size = observed), alpha = 0.85, color = dir_color) +
    facet_wrap(~db, ncol = 3, drop = FALSE, scales = "free_y",
               labeller = labeller(db = label_wrap_gen(width = 16))) +
    scale_y_discrete(labels = function(x) str_wrap(gsub("\\|\\|\\|.*$", "", x), width = wrap_width),
                     expand = expansion(mult = 0.06, add = 1.2)) +
    scale_size(range = size_range_ed, name = "Observed") +
    labs(
      x = "-log10(P-adjusted)",
      y = NULL,
      title = paste0("miRPathDB ", dir_type, " (Jaccard+rrvgo reduced): ", plot_label)
    ) +
    theme_pub(p_ed)

  # Save multipanel plot in enriched_depleted/
  dir_tag <- tolower(dir_type)
  out_base <- file.path(out_fig_dir, "enriched_depleted",
                        paste0("SurvivalRank_Bubble_miRPathDB_jaccard_", dir_tag, "_",
                               safe_name(plot_label), if (stamp) paste0("_", ts) else ""))
  save_plot(p_dir, out_base, n_rows = 1, width_mult = 4.2, height_mult = 2.2)

  plots_list[[dir_type]] <- p_dir

  # ---- Per-database plots ----
  p_db_preset <- presets[["per_database"]]
  p <- p_db_preset  # Update p for per-database saves
  size_range_db <- c(p_db_preset$point * 0.8, p_db_preset$point * 2.8)

  for (db_name in unique(df_dir$db)) {
    df_db <- df_dir %>% filter(db == db_name)
    if (nrow(df_db) == 0) next

    # Re-order for single DB
    df_db$term_plot <- factor(df_db$term, levels = rev(unique(df_db$term)))

    p_db <- ggplot(df_db, aes(x = logp, y = term_plot)) +
      geom_point(aes(size = observed), alpha = 0.85, color = dir_color) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width),
                       expand = expansion(mult = 0.06, add = 1.0)) +
      scale_size(range = size_range_db, name = "Observed") +
      labs(
        x = "-log10(P-adjusted)",
        y = NULL,
        title = paste0(db_name, " - ", dir_type, " (Jaccard+rrvgo)")
      ) +
      theme_pub(p_db_preset)

    # Save per-database plot
    db_safe <- safe_name(db_name)
    out_base_db <- file.path(out_fig_dir, "per_database",
                             paste0("SurvivalRank_", db_safe, "_jaccard_", dir_tag, "_",
                                    safe_name(plot_label), if (stamp) paste0("_", ts) else ""))
    save_plot(p_db, out_base_db, n_rows = 1, width_mult = 2.0,
              height_mult = max(1.2, nrow(df_db) * 0.24))

    # Store for potential use
    key <- paste(dir_type, db_name, sep = "_")
    plots_per_db[[key]] <- p_db
  }
}

cat("\nBubble plots exported to:\n")
cat("  Multipanel:", file.path(out_fig_dir, "enriched_depleted"), "\n")
cat("  Per-database:", file.path(out_fig_dir, "per_database"), "\n")

# ============================================================================
# DONE
# ============================================================================
cat("============================================================================\n")
cat("DONE - Pathway redundancy reduction + visualization\n")
cat("============================================================================\n")
cat("Reduced table:", reduced_fp, "\n")
cat("Full annotated:", full_fp, "\n")
cat("Summary:", ifelse(exists("summary_fp"), summary_fp, "(none)"), "\n")
cat("Summary MD:", md_fp, "\n")
cat("Figures:", out_fig_dir, "\n")
cat("  - Multipanel plots:", file.path(out_fig_dir, "enriched_depleted"), "\n")
cat("  - Per-database plots:", file.path(out_fig_dir, "per_database"), "\n")
cat("Log:", log_fp, "\n")
cat("============================================================================\n")
