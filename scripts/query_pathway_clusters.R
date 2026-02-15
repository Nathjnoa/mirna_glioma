#!/usr/bin/env Rscript
# ============================================================================
# Query Pathway Clustering Results
#
# Explore which pathways were grouped together in the redundancy reduction
# and which representative was selected for each cluster
# ============================================================================

library(dplyr)
library(tidyr)

# ---- Load full annotated table ----
full_file <- "results_29s/tables/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/reduced_jaccard/miEAA_GSEA_full_annotated_jaccard_20260215_142223.tsv"

if (!file.exists(full_file)) {
  # Try to find latest file
  full_file <- list.files("results_29s/tables/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/reduced_jaccard/",
                          pattern = "miEAA_GSEA_full_annotated_jaccard.*\\.tsv$",
                          full.names = TRUE)
  if (length(full_file) == 0) stop("No full annotated file found")
  full_file <- full_file[1]
}

cat("Loading:", full_file, "\n\n")
df <- read.delim(full_file, stringsAsFactors = FALSE)

# ---- Filter to clustered pathways (Jaccard method) ----
df_jaccard <- df %>%
  filter(reduction_method %in% c("Jaccard_hclust_average_t0.25", "Jaccard_hclust_average_t0.50"),
         !is.na(cluster))

cat("Total pathways with Jaccard clustering:", nrow(df_jaccard), "\n")
cat("  Representatives:", sum(df_jaccard$cluster_rep), "\n")
cat("  Absorbed:", sum(!df_jaccard$cluster_rep), "\n\n")

# ---- Summary by database ----
cluster_summary <- df_jaccard %>%
  group_by(database, enrichment) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    n_pathways = n(),
    n_representatives = sum(cluster_rep),
    n_absorbed = sum(!cluster_rep),
    .groups = "drop"
  ) %>%
  arrange(database, enrichment)

cat("=== Clustering Summary by Database ===\n")
print(as.data.frame(cluster_summary), row.names = FALSE)
cat("\n")

# ---- Function to show cluster details ----
show_cluster <- function(db, direction, cluster_id) {
  cat(sprintf("\n=== Cluster %s: %s %s ===\n", cluster_id, db, direction))

  cluster_df <- df_jaccard %>%
    filter(database == db, enrichment == direction, cluster == cluster_id) %>%
    arrange(cluster_rep %>% desc(), qval) %>%
    select(term, qval, cluster_rep)

  if (nrow(cluster_df) == 0) {
    cat("  No pathways found\n")
    return(invisible(NULL))
  }

  # Representative
  rep_pathway <- cluster_df %>% filter(cluster_rep) %>% pull(term)
  rep_qval <- cluster_df %>% filter(cluster_rep) %>% pull(qval)

  cat(sprintf("  Representative: %s (Q = %.2e)\n", rep_pathway, rep_qval))
  cat(sprintf("  Absorbed: %d pathways\n", sum(!cluster_df$cluster_rep)))

  if (nrow(cluster_df) > 1) {
    cat("\n  Absorbed pathways:\n")
    absorbed <- cluster_df %>% filter(!cluster_rep)
    for (i in 1:nrow(absorbed)) {
      cat(sprintf("    - %s (Q = %.2e)\n", absorbed$term[i], absorbed$qval[i]))
    }
  }

  invisible(cluster_df)
}

# ---- Show examples ----
cat("\n=== Example Clusters ===\n")

# KEGG depleted cluster 1 (largest cluster)
show_cluster("KEGG", "depleted", 1)

# Reactome enriched cluster 1
if (any(df_jaccard$database == "Reactome" & df_jaccard$enrichment == "enriched")) {
  show_cluster("Reactome", "enriched", 1)
}

# ---- Export cluster membership table ----
cluster_table <- df_jaccard %>%
  select(database, enrichment, cluster, term, qval, cluster_rep, reduction_method) %>%
  arrange(database, enrichment, cluster, cluster_rep %>% desc(), qval)

out_file <- "results_29s/tables/SurvivalRank_GSEA/SurvivalRank_CoxZ_miRPathDB/reduced_jaccard/cluster_membership.tsv"
write.table(cluster_table, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nCluster membership table exported to:\n")
cat(" ", out_file, "\n")

cat("\n=== Usage Examples ===\n")
cat("# Show specific cluster:\n")
cat('show_cluster("KEGG", "depleted", 1)\n')
cat('show_cluster("Reactome", "enriched", 5)\n')
cat("\n# Show all clusters for a database:\n")
cat('df_jaccard %>% filter(database == "KEGG", enrichment == "depleted") %>%\n')
cat('  group_by(cluster) %>%\n')
cat('  summarise(n_pathways = n(), representative = term[cluster_rep][1])\n')
