#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(igraph)
  library(dplyr)
  library(readr)
})

options(stringsAsFactors = FALSE)

# ===================== PATHS ===================== #
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"
out_base <- file.path(base_dir, "1_1_Centrality_comparison_rho")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# ===================== NETWORKS (manual, readable) ===================== #
networks <- list(
  PCC   = c("ROSMAP_PCC_counts_AD_topN200000.tsv",
            "ROSMAP_PCC_counts_control_topN200000.tsv"),
  DLPFC = c("ROSMAP_DLPFC_counts_AD_topN200000.tsv",
            "ROSMAP_DLPFC_counts_control_topN200000.tsv"),
  HCN   = c("ROSMAP_HCN_counts_AD_topN200000.tsv",
            "ROSMAP_HCN_counts_control_topN200000.tsv"),
  MAYO_CRB = c("Mayo_CRB_counts_AD_topN200000.tsv",
               "Mayo_CRB_counts_control_topN200000.tsv"),
  MAYO_TC  = c("Mayo_TC_counts_AD_topN200000.tsv",
               "Mayo_TC_counts_control_topN200000.tsv")
)

centralities <- c("degree", "pagerank", "betweenness")

# ============================================================
# 1) LOAD (edgelist -> igraph)
# ============================================================
load_graph <- function(file) {
  df <- vroom(
    file.path(base_dir, "1_Edgelist", file),
    col_names = c("GenA", "GenB", "MI")
  )
  
  # ---- FIX: drop any non-gene placeholder like "gene1" ----
  # This usually comes from a header/extra row sneaking in.
  df <- df %>%
    filter(
      !is.na(GenA), !is.na(GenB),
      GenA != "gene1", GenB != "gene1",
      GenA != "Gene1", GenB != "Gene1"
    )
  
  graph_from_data_frame(df, directed = FALSE)
}

# ============================================================
# 2) COMPUTE CENTRALITY (igraph -> named numeric vector)
# ============================================================
compute_centrality <- function(g, type) {
  if (type == "degree") {
    degree(g)
  } else if (type == "pagerank") {
    page_rank(g)$vector
  } else if (type == "betweenness") {
    betweenness(g, directed = FALSE, normalized = TRUE)
  } else {
    stop("Centrality not supported: ", type)
  }
}

# ============================================================
# 2b) Percentile ranks (order-only, robust to long tails)
# ============================================================
to_percentile <- function(x) {
  rank(x, ties.method = "average") / length(x)  # in (0,1]
}

# ============================================================
# 3) BUILD TABLES (AD vec + CTL vec -> FULL + COMMON)
#    Adds: percentile ranks + delta_percentile
# ============================================================
build_tables <- function(vec_AD, vec_CTL) {
  
  df_full <- full_join(
    tibble(gene = names(vec_AD),  AD  = as.numeric(vec_AD)),
    tibble(gene = names(vec_CTL), CTL = as.numeric(vec_CTL)),
    by = "gene"
  ) %>%
    # ---- FIX: ensure placeholder doesn't propagate even if present ----
  filter(!is.na(gene), gene != "gene1", gene != "Gene1")
  
  df_common <- df_full %>%
    filter(!is.na(AD) & !is.na(CTL)) %>%
    mutate(
      # Raw change
      delta = AD - CTL,
      
      # Descriptive Z-score of delta (NOT a p-value)
      z_change = (delta - mean(delta)) / sd(delta),
      
      # Percentile ranks within each condition (order-based)
      perc_AD  = to_percentile(AD),
      perc_CTL = to_percentile(CTL),
      
      # Change in percentile rank (interpretable, distribution-free)
      delta_perc = perc_AD - perc_CTL,
      
      HigherNetwork = ifelse(AD > CTL, "AD", "Control")
    )
  
  list(full = df_full, common = df_common)
}

# ============================================================
# 4) SAVE TABLES (same filenames as before)
# ============================================================
save_tables <- function(tables, net_name, centrality) {
  
  out_dir <- file.path(out_base, centrality)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_tsv(
    tables$full,
    file.path(out_dir, paste0(net_name, "_", centrality, "_FULL.tsv"))
  )
  
  write_tsv(
    tables$common,
    file.path(out_dir, paste0(net_name, "_", centrality, "_COMMON.tsv"))
  )
}

# ============================================================
# 5) Spearman per network (common genes only)
#    Saves one TSV per network with rho and p-value per centrality.
# ============================================================
save_spearman <- function(net_name, spearman_df) {
  out_dir <- file.path(out_base, "spearman")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_tsv(
    spearman_df,
    file.path(out_dir, paste0(net_name, "_spearman.tsv"))
  )
}

# ============================================================
# RUN (load -> compute -> tables -> save) + Spearman
# ============================================================
for (net in names(networks)) {
  
  cat("\n=============================\n")
  cat("Network:", net, "\n")
  cat("=============================\n")
  
  # 1) Load graphs
  g_AD  <- load_graph(networks[[net]][1])
  g_CTL <- load_graph(networks[[net]][2])
  
  cat("Loaded graphs | AD nodes:", vcount(g_AD), "edges:", ecount(g_AD),
      "| CTL nodes:", vcount(g_CTL), "edges:", ecount(g_CTL), "\n")
  
  # Collect Spearman results for this network
  spearman_res <- tibble(
    network = character(),
    centrality = character(),
    n_common_genes = integer(),
    spearman_rho = double(),
    spearman_pvalue = double()
  )
  
  # 2) For each centrality: compute -> tables -> save
  for (cent in centralities) {
    
    cat("  - Centrality:", cent, "\n")
    
    vec_AD  <- compute_centrality(g_AD,  cent)
    vec_CTL <- compute_centrality(g_CTL, cent)
    
    cat("    computed | AD n=", length(vec_AD), "CTL n=", length(vec_CTL), "\n")
    
    tables <- build_tables(vec_AD, vec_CTL)
    save_tables(tables, net, cent)
    
    # Spearman on common genes (rank stability between conditions)
    dfc <- tables$common
    sp <- suppressWarnings(cor.test(dfc$AD, dfc$CTL, method = "spearman", exact = FALSE))
    
    spearman_res <- bind_rows(
      spearman_res,
      tibble(
        network = net,
        centrality = cent,
        n_common_genes = nrow(dfc),
        spearman_rho = unname(sp$estimate),
        spearman_pvalue = sp$p.value
      )
    )
    
    cat("    saved | FULL:", nrow(tables$full), "COMMON:", nrow(tables$common), "\n")
    cat("    spearman | rho=", round(unname(sp$estimate), 4), " p=", signif(sp$p.value, 3), "\n")
  }
  
  # Save Spearman summary for this network
  save_spearman(net, spearman_res)
}

cat("\nDONE\n")
