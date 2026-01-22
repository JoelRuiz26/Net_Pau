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
out_base <- file.path(base_dir, "1_1_Centrality_comparison")
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
    file.path(base_dir,"1_Edgelist", file),
    col_names = c("GenA", "GenB", "MI")
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
# 3) BUILD TABLES (AD vec + CTL vec -> FULL + COMMON)
# ============================================================
build_tables <- function(vec_AD, vec_CTL) {
  
  df_full <- full_join(
    tibble(gene = names(vec_AD),  AD  = as.numeric(vec_AD)),
    tibble(gene = names(vec_CTL), CTL = as.numeric(vec_CTL)),
    by = "gene"
  )
  
  df_common <- df_full %>%
    filter(!is.na(AD) & !is.na(CTL)) %>%
    mutate(
      delta = AD - CTL,
      z_change = (delta - mean(delta)) / sd(delta),
      HigherNetwork = ifelse(AD > CTL, "AD", "Control")
    )
  
  list(full = df_full, common = df_common)
}

# ============================================================
# 4) SAVE TABLES
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
# RUN (load -> compute -> tables -> save)
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
  
  # 2) For each centrality: compute -> tables -> save
  for (cent in centralities) {
    
    cat("  - Centrality:", cent, "\n")
    
    vec_AD  <- compute_centrality(g_AD,  cent)
    vec_CTL <- compute_centrality(g_CTL, cent)
    
    # “que se vea el resultado” (mini check)
    cat("    computed | AD n=", length(vec_AD), "CTL n=", length(vec_CTL), "\n")
    
    tables <- build_tables(vec_AD, vec_CTL)
    save_tables(tables, net, cent)
    
    cat("    saved | FULL:", nrow(tables$full), "COMMON:", nrow(tables$common), "\n")
  }
}

cat("\nDONE\n")
