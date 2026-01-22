#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(igraph)
  library(dplyr)
  library(readr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 0) Paths and inputs
# ============================================================
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"
out_base <- file.path(base_dir, "1_1_Centrality_comparison")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Networks (manual, readable)
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

# Permutation settings (keep it simple; increase on HPC if needed)
n_perm_labels <- 1000
n_perm_rewire <- 200          # for betweenness you may want 50-100 if too slow
rewire_mult   <- 10           # niter = rewire_mult * ecount(g)
top_n         <- 50           # report top genes by |Δ| (change to 10/100 as needed)
set.seed(1)

# ============================================================
# 1) Load edgelist -> igraph
# ============================================================
load_graph <- function(file) {
  df <- vroom(
    file.path(base_dir, file),
    col_names = c("GenA", "GenB", "MI")
  )
  graph_from_data_frame(df, directed = FALSE)
}

# ============================================================
# 2) Compute centrality (igraph -> named numeric vector)
# ============================================================
compute_centrality <- function(g, type) {
  if (type == "degree") {
    degree(g)
  } else if (type == "pagerank") {
    page_rank(g)$vector
  } else if (type == "betweenness") {
    # normalized=TRUE helps comparability across graphs, but still requires a null for significance
    betweenness(g, directed = FALSE, normalized = TRUE)
  } else {
    stop("Centrality not supported: ", type)
  }
}

# ============================================================
# 3) Build FULL + COMMON tables and basic deltas
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
      abs_delta = abs(delta),
      z_change = (delta - mean(delta)) / sd(delta),
      HigherNetwork = ifelse(AD > CTL, "AD", "Control")
    )
  
  list(full = df_full, common = df_common)
}

# ============================================================
# 4) Empirical p-values: label permutation (topology fixed)
#    - This tests whether a gene-specific |Δ| is larger than expected
#      by random re-assignment of centrality values to gene labels.
# ============================================================
p_empirical_label <- function(obs_AD, obs_CTL, n_perm = 1000) {
  # obs_AD/obs_CTL are numeric vectors aligned to the SAME gene order
  abs_obs <- abs(obs_AD - obs_CTL)
  m <- length(abs_obs)
  
  counts <- integer(m)
  
  for (b in seq_len(n_perm)) {
    perm_AD  <- sample(obs_AD,  size = m, replace = FALSE)
    perm_CTL <- sample(obs_CTL, size = m, replace = FALSE)
    abs_perm <- abs(perm_AD - perm_CTL)
    counts <- counts + (abs_perm >= abs_obs)
  }
  
  # add-one smoothing to avoid zeros
  (counts + 1) / (n_perm + 1)
}

# ============================================================
# 5) Empirical p-values: degree-preserving rewiring (topology null)
#    - This tests whether observed |Δ| is larger than expected under
#      random graphs with the SAME degree sequence (per condition).
# ============================================================
p_empirical_rewire <- function(g_AD, g_CTL, genes_common, type,
                               n_perm = 200, rewire_mult = 10) {
  
  # observed centrality (aligned to genes_common)
  obs_AD_all  <- compute_centrality(g_AD,  type)
  obs_CTL_all <- compute_centrality(g_CTL, type)
  
  obs_AD  <- as.numeric(obs_AD_all[genes_common])
  obs_CTL <- as.numeric(obs_CTL_all[genes_common])
  
  abs_obs <- abs(obs_AD - obs_CTL)
  m <- length(abs_obs)
  
  counts <- integer(m)
  
  niter_AD  <- max(1, rewire_mult * ecount(g_AD))
  niter_CTL <- max(1, rewire_mult * ecount(g_CTL))
  
  for (b in seq_len(n_perm)) {
    gA <- rewire(g_AD,  with = keeping_degseq(niter = niter_AD))
    gC <- rewire(g_CTL, with = keeping_degseq(niter = niter_CTL))
    
    cA_all <- compute_centrality(gA, type)
    cC_all <- compute_centrality(gC, type)
    
    cA <- as.numeric(cA_all[genes_common])
    cC <- as.numeric(cC_all[genes_common])
    
    abs_perm <- abs(cA - cC)
    counts <- counts + (abs_perm >= abs_obs)
  }
  
  (counts + 1) / (n_perm + 1)
}

# ============================================================
# 6) Save outputs per centrality
# ============================================================
save_outputs <- function(net_name, centrality, tables, top_tbl) {
  out_dir <- file.path(out_base, centrality)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_tsv(tables$full,   file.path(out_dir, paste0(net_name, "_", centrality, "_FULL.tsv")))
  write_tsv(tables$common, file.path(out_dir, paste0(net_name, "_", centrality, "_COMMON.tsv")))
  write_tsv(top_tbl,       file.path(out_dir, paste0(net_name, "_", centrality, "_TOP.tsv")))
}

# ============================================================
# 7) Run: load -> compute -> tables -> p-values -> top -> save
# ============================================================
for (net in names(networks)) {
  
  cat("\n=============================\n")
  cat("Network:", net, "\n")
  cat("=============================\n")
  
  # Load graphs
  g_AD  <- load_graph(networks[[net]][1])
  g_CTL <- load_graph(networks[[net]][2])
  
  cat("Loaded graphs | AD nodes:", vcount(g_AD), "edges:", ecount(g_AD),
      "| CTL nodes:", vcount(g_CTL), "edges:", ecount(g_CTL), "\n")
  
  for (cent in centralities) {
    
    cat("  - Centrality:", cent, "\n")
    
    # Compute observed centralities
    vec_AD  <- compute_centrality(g_AD,  cent)
    vec_CTL <- compute_centrality(g_CTL, cent)
    
    # Build tables (FULL + COMMON + deltas)
    tables <- build_tables(vec_AD, vec_CTL)
    df_common <- tables$common
    
    # Align vectors for permutation tests (common genes only)
    genes_common <- df_common$gene
    obs_AD  <- df_common$AD
    obs_CTL <- df_common$CTL
    
    cat("    common genes:", length(genes_common), "\n")
    
    # Empirical p-values (label permutation)
    cat("    p_emp label permutations:", n_perm_labels, "\n")
    p_label <- p_empirical_label(obs_AD, obs_CTL, n_perm = n_perm_labels)
    
    # Empirical p-values (degree-preserving rewiring)
    cat("    p_emp degree-preserving rewiring:", n_perm_rewire, "\n")
    p_rewire <- p_empirical_rewire(
      g_AD = g_AD, g_CTL = g_CTL,
      genes_common = genes_common,
      type = cent,
      n_perm = n_perm_rewire,
      rewire_mult = rewire_mult
    )
    
    # Add p-values (+ optional FDR per network/centrality)
    df_common <- df_common %>%
      mutate(
        p_emp_label  = p_label,
        p_emp_rewire = p_rewire,
        FDR_label    = p.adjust(p_emp_label,  method = "BH"),
        FDR_rewire   = p.adjust(p_emp_rewire, method = "BH")
      )
    
    # Update tables$common with p-values
    tables$common <- df_common
    
    # Top genes by |Δ| (report both empirical p-values)
    top_tbl <- df_common %>%
      arrange(desc(abs_delta)) %>%
      slice_head(n = top_n) %>%
      dplyr::select(gene, AD, CTL, delta, abs_delta, z_change, HigherNetwork,
             p_emp_label, FDR_label, p_emp_rewire, FDR_rewire)
    
    # Save everything
    save_outputs(net, cent, tables, top_tbl)
    
    cat("    saved | FULL:", nrow(tables$full),
        "COMMON:", nrow(tables$common),
        "TOP:", nrow(top_tbl), "\n")
  }
}

cat("\nDONE\n")
