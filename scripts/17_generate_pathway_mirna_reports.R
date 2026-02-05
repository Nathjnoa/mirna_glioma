#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# 17_generate_pathway_mirna_reports.R
# Generate reports of top pathways with miRNAs annotated by rank and statistic
#
# Output format for miRNAs column:
#   hsa-miR-451a (#3, logFC=-0.85); hsa-miR-148a-3p (#6, logFC=-0.72); ...
#
# Two report types:
#   1. DE-based GSEA: uses logFC as ranking statistic
#   2. Survival-based GSEA: uses Cox z-score as ranking statistic
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# --- CLI parsing ---
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[[hit + 1]]
}

# Parameters
spec_fp      <- get_arg("--spec", "config/de_specs.csv")
de_root      <- get_arg("--de_root", "results/DE")
gsea_root    <- get_arg("--gsea_root", "results/tables/MiEAA_GSEA")
surv_root    <- get_arg("--surv_root", "results/tables/SurvivalRank_GSEA")
out_dir      <- get_arg("--out_dir", "results/tables/reports")
run_tag      <- get_arg("--run_tag", "A_conservative")
top_n        <- as.integer(get_arg("--top_n", "12"))
surv_run_tag <- get_arg("--surv_run_tag", "SurvivalRank_CoxZ_miRPathDB")

# Create output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Timestamp
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Log setup
log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE)
log_fp <- file.path(log_dir, paste0("pathway_mirna_reports_", ts, ".txt"))
log_con <- file(log_fp, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat("=== Pathway-miRNA Reports Generator ===\n")
cat("Timestamp:", ts, "\n")
cat("Parameters:\n")
cat("  spec_fp:", spec_fp, "\n")
cat("  de_root:", de_root, "\n")
cat("  gsea_root:", gsea_root, "\n")
cat("  surv_root:", surv_root, "\n")
cat("  out_dir:", out_dir, "\n")
cat("  run_tag:", run_tag, "\n")
cat("  top_n:", top_n, "\n")
cat("  surv_run_tag:", surv_run_tag, "\n\n")

# --- Helper functions ---

read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# Find most recent file matching pattern
find_latest_file <- function(dir_path, pattern) {
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  # Sort by modification time, get most recent
  info <- file.info(files)
  files[which.max(info$mtime)]
}

# Format miRNA with rank and statistic
format_mirna_with_stats <- function(mirnas_str, rank_df, stat_col, stat_name) {
  if (is.na(mirnas_str) || mirnas_str == "") return("")

  mirnas <- trimws(unlist(strsplit(mirnas_str, ";", fixed = TRUE)))
  mirnas <- mirnas[nzchar(mirnas)]

  formatted <- sapply(mirnas, function(m) {
    idx <- which(rank_df$feature_id == m)
    if (length(idx) == 0) {
      # Try matching without hsa- prefix or with variations
      idx <- which(grepl(m, rank_df$feature_id, fixed = TRUE))
    }
    if (length(idx) == 0) {
      return(m)  # Return as-is if not found
    }
    idx <- idx[1]  # Take first match
    rank_val <- rank_df$rank[idx]
    stat_val <- rank_df[[stat_col]][idx]
    sprintf("%s (#%d, %s=%.3f)", m, rank_val, stat_name, stat_val)
  })

  paste(formatted, collapse = "; ")
}

# -----------------------------------------------------------------------------
# PART 1: DE-based GSEA Reports
# -----------------------------------------------------------------------------

cat("\n========== DE-based GSEA Reports ==========\n\n")

# Read spec file
if (!file.exists(spec_fp)) {
  cat("WARNING: Spec file not found:", spec_fp, "\n")
  spec_rows <- data.frame()
} else {
  spec_rows <- read_table(spec_fp)
  cat("Loaded", nrow(spec_rows), "comparisons from spec\n")
}

# Process each comparison
de_results_list <- list()

for (i in seq_len(nrow(spec_rows))) {
  analysis_id <- spec_rows$analysis_id[i]
  cat("\n--- Processing:", analysis_id, "---\n")

  # Find DE results file
  de_dir <- file.path(de_root, analysis_id)
  if (!dir.exists(de_dir)) {
    # Try alternative naming with safe characters
    de_dir_alt <- file.path(de_root, gsub("[^A-Za-z0-9_-]", "_", analysis_id))
    if (dir.exists(de_dir_alt)) {
      de_dir <- de_dir_alt
    } else {
      cat("  DE directory not found, skipping\n")
      next
    }
  }

  de_file <- find_latest_file(de_dir, "^DE_.*_all_.*\\.tsv$")
  if (is.null(de_file)) {
    cat("  DE file not found, skipping\n")
    next
  }
  cat("  DE file:", basename(de_file), "\n")

  # Read DE results and compute ranks
  de_df <- read_table(de_file)
  if (is.null(de_df) || nrow(de_df) == 0) {
    cat("  Empty DE results, skipping\n")
    next
  }

  # Rank by |logFC| descending (most changed first)
  de_df$abs_logFC <- abs(de_df$logFC)
  de_df <- de_df[order(-de_df$abs_logFC), ]
  de_df$rank <- seq_len(nrow(de_df))

  cat("  Loaded", nrow(de_df), "features, top logFC:",
      round(de_df$logFC[1], 3), "\n")

  # Find GSEA results
  gsea_dir <- file.path(gsea_root, analysis_id, run_tag)
  if (!dir.exists(gsea_dir)) {
    # Try with safe name
    analysis_safe <- gsub(">", "", gsub("[^A-Za-z0-9_-]", "_", analysis_id))
    gsea_dir <- file.path(gsea_root, analysis_safe, run_tag)
  }

  gsea_file <- find_latest_file(gsea_dir, "^MiEAA_GSEA_Qlt.*\\.tsv$")
  if (is.null(gsea_file)) {
    # Try all results if filtered not available
    gsea_file <- find_latest_file(gsea_dir, "^MiEAA_GSEA_all.*\\.tsv$")
  }

  if (is.null(gsea_file)) {
    cat("  GSEA file not found in", gsea_dir, ", skipping\n")
    next
  }
  cat("  GSEA file:", basename(gsea_file), "\n")

  # Read GSEA results
  gsea_df <- read_table(gsea_file)
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    cat("  Empty GSEA results, skipping\n")
    next
  }

  # Exclude GO Molecular Function (not biologically relevant for miRNA analysis)
  gsea_df <- gsea_df[!grepl("Molecular function", gsea_df$Category, ignore.case = TRUE), ]

  # Get top N per category and enrichment direction
  categories <- unique(gsea_df$Category)

  for (cat_name in categories) {
    for (enrich_dir in c("enriched", "depleted")) {
      subset_df <- gsea_df[gsea_df$Category == cat_name &
                           gsea_df$Enrichment == enrich_dir, ]

      if (nrow(subset_df) == 0) next

      # Sort by Q-value and take top N
      subset_df <- subset_df[order(subset_df$`Q-value`), ]
      subset_df <- head(subset_df, top_n)

      # Format miRNAs with rank and logFC
      mirna_col <- grep("miRNA|precursor", names(subset_df),
                        ignore.case = TRUE, value = TRUE)[1]
      if (is.na(mirna_col)) mirna_col <- "miRNAs/precursors"

      subset_df$miRNAs_annotated <- sapply(subset_df[[mirna_col]], function(m) {
        format_mirna_with_stats(m, de_df, "logFC", "logFC")
      })

      # Build output row
      for (j in seq_len(nrow(subset_df))) {
        de_results_list[[length(de_results_list) + 1]] <- data.frame(
          comparison = analysis_id,
          Category = cat_name,
          Subcategory = subset_df$Subcategory[j],
          Enrichment = enrich_dir,
          Q_value = subset_df$`Q-value`[j],
          Observed = subset_df$Observed[j],
          miRNAs_with_rank = subset_df$miRNAs_annotated[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  cat("  Added pathways to report\n")
}

# Combine and save DE-based report
if (length(de_results_list) > 0) {
  de_report <- do.call(rbind, de_results_list)

  out_file <- file.path(out_dir,
                        paste0("GSEA_DE_top", top_n, "_annotated_", ts, ".csv"))
  write.csv(de_report, out_file, row.names = FALSE)
  cat("\n\nDE-based report saved:", out_file, "\n")
  cat("Total rows:", nrow(de_report), "\n")
} else {
  cat("\nNo DE-based results to report\n")
}

# -----------------------------------------------------------------------------
# PART 2: Survival-based GSEA Reports
# -----------------------------------------------------------------------------

cat("\n\n========== Survival-based GSEA Reports ==========\n\n")

# Find Cox results
surv_dir <- file.path(surv_root, surv_run_tag)
if (!dir.exists(surv_dir)) {
  cat("Survival directory not found:", surv_dir, "\n")
} else {
  cox_file <- find_latest_file(surv_dir, "^Cox_univariate_miRNA_all.*\\.tsv$")

  if (is.null(cox_file)) {
    cat("Cox results file not found\n")
  } else {
    cat("Cox file:", basename(cox_file), "\n")

    # Read Cox results
    cox_df <- read_table(cox_file)

    if (!is.null(cox_df) && nrow(cox_df) > 0) {
      # Rank by |z| descending
      cox_df$abs_z <- abs(cox_df$z)
      cox_df <- cox_df[order(-cox_df$abs_z), ]
      cox_df$rank <- seq_len(nrow(cox_df))

      cat("Loaded", nrow(cox_df), "features from Cox analysis\n")
      cat("Top z-score:", round(cox_df$z[1], 3), "\n")

      # Find reduced GSEA results (after rrvgo)
      reduced_dir <- file.path(surv_dir, "reduced_rrvgo")

      # Try to find GSEA results
      gsea_file <- NULL
      if (dir.exists(reduced_dir)) {
        gsea_file <- find_latest_file(reduced_dir, "^miEAA_GSEA.*reduced.*\\.tsv$")
      }
      if (is.null(gsea_file)) {
        gsea_file <- find_latest_file(surv_dir, "^miEAA_GSEA.*\\.tsv$")
      }

      if (!is.null(gsea_file)) {
        cat("GSEA file:", basename(gsea_file), "\n")

        gsea_df <- read_table(gsea_file)

        if (!is.null(gsea_df) && nrow(gsea_df) > 0) {
          surv_results_list <- list()

          # Exclude GO Molecular Function
          gsea_df <- gsea_df[!grepl("Molecular function", gsea_df$Category, ignore.case = TRUE), ]

          categories <- unique(gsea_df$Category)

          for (cat_name in categories) {
            for (enrich_dir in c("enriched", "depleted")) {
              subset_df <- gsea_df[gsea_df$Category == cat_name &
                                   gsea_df$Enrichment == enrich_dir, ]

              if (nrow(subset_df) == 0) next

              # Sort by Q-value and take top N
              q_col <- grep("Q-value|q-value|Q_value", names(subset_df),
                           ignore.case = TRUE, value = TRUE)[1]
              if (is.na(q_col)) q_col <- "Q-value"

              subset_df <- subset_df[order(subset_df[[q_col]]), ]
              subset_df <- head(subset_df, top_n)

              # Format miRNAs with rank and z-score
              mirna_col <- grep("miRNA|precursor", names(subset_df),
                               ignore.case = TRUE, value = TRUE)[1]
              if (is.na(mirna_col)) mirna_col <- "miRNAs/precursors"

              subset_df$miRNAs_annotated <- sapply(subset_df[[mirna_col]], function(m) {
                format_mirna_with_stats(m, cox_df, "z", "z")
              })

              # Build output row
              for (j in seq_len(nrow(subset_df))) {
                surv_results_list[[length(surv_results_list) + 1]] <- data.frame(
                  analysis = "SurvivalRank_CoxZ",
                  Category = cat_name,
                  Subcategory = subset_df$Subcategory[j],
                  Enrichment = enrich_dir,
                  Q_value = subset_df[[q_col]][j],
                  Observed = subset_df$Observed[j],
                  miRNAs_with_rank = subset_df$miRNAs_annotated[j],
                  stringsAsFactors = FALSE
                )
              }
            }
          }

          # Combine and save survival-based report
          if (length(surv_results_list) > 0) {
            surv_report <- do.call(rbind, surv_results_list)

            out_file <- file.path(out_dir,
                                  paste0("GSEA_survGSEA_top", top_n, "_annotated_", ts, ".csv"))
            write.csv(surv_report, out_file, row.names = FALSE)
            cat("\nSurvival-based report saved:", out_file, "\n")
            cat("Total rows:", nrow(surv_report), "\n")
          }
        }
      } else {
        cat("GSEA results file not found for survival analysis\n")
      }
    }
  }
}

# Close log
cat("\n\nDone. Log saved:", log_fp, "\n")
sink(type = "message")
sink()
close(log_con)
