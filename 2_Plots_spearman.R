#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(biomaRt)
  library(grid)    # unit()
  library(scales)  # pseudo_log_trans, label_number
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau/1_1_Centrality_comparison_rho"

centralities <- c("degree", "pagerank", "betweenness")
networks     <- c("PCC", "DLPFC", "HCN", "MAYO_CRB", "MAYO_TC")

TOP_PROP <- 0.001  # top 1% absolute by |Δ percentile rank|

OUT_DIR <- file.path(base_dir, "QC_plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CENT_LABELS <- c(degree = "Degree", pagerank = "PageRank", betweenness = "Betweenness")

# ===================== TRANSFORM (log-like, safe with zeros) ===================== #
SIGMA <- 1e-8
plog  <- scales::pseudo_log_trans(base = 10, sigma = SIGMA)

# Breaks that won't crash with free scales, and helps prevent label overlap
safe_breaks <- function(n = 2) {
  function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(numeric(0))
    rng <- range(x)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) return(rng[1])
    scales::pretty_breaks(n = n)(x)
  }
}

# Compact axis labels (k, M, B, T) using label_number(scale_cut=...)
fmt_si <- scales::label_number(
  accuracy = 1,
  scale_cut = scales::cut_si("")
)

cent_pretty <- function(cent_name) {
  out <- unname(CENT_LABELS[cent_name])
  if (is.na(out) || length(out) != 1) return(cent_name)
  out
}

# ===================== HELPERS ===================== #
read_common_tbl <- function(net, cent) {
  f <- file.path(base_dir, cent, paste0(net, "_", cent, "_COMMON.tsv"))
  if (!file.exists(f)) return(NULL)
  
  df <- vroom::vroom(f, show_col_types = FALSE)
  
  need <- c("gene", "AD", "CTL", "perc_AD", "perc_CTL", "delta_perc", "HigherNetwork")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", f, ": ", paste(miss, collapse = ", "))
  
  df %>%
    mutate(network = net, centrality = cent) %>%
    dplyr::select(network, centrality, everything())
}

read_spearman_tbl <- function(net) {
  f <- file.path(base_dir, "spearman", paste0(net, "_spearman.tsv"))
  if (!file.exists(f)) return(NULL)
  
  df <- vroom::vroom(f, show_col_types = FALSE)
  
  need <- c("network", "centrality", "n_common_genes", "spearman_rho", "spearman_pvalue")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", f, ": ", paste(miss, collapse = ", "))
  
  df %>%
    mutate(
      network    = as.character(network),
      centrality = as.character(centrality)
    ) %>%
    dplyr::select(network, centrality, n_common_genes, spearman_rho, spearman_pvalue)
}

# --- Ensembl -> HGNC symbol (biomaRt) ---
map_ensembl_to_hgnc <- function(ens_ids) {
  ens_ids <- unique(ens_ids)
  ens_ids <- sub("\\..*$", "", ens_ids)
  
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = ens_ids,
    mart       = mart
  ) %>%
    mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", NA_character_, hgnc_symbol)) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
}

# ===================== LOAD ALL COMMON TABLES ===================== #
all_common <- purrr::map_dfr(
  networks,
  \(net) purrr::map_dfr(centralities, \(cent) read_common_tbl(net, cent))
)
if (nrow(all_common) == 0) stop("No *_COMMON.tsv files were loaded. Check paths and filenames.")

# ===================== BIOMART ANNOTATION ===================== #
cat("Annotating ENSG -> HGNC with biomaRt...\n")
ann <- map_ensembl_to_hgnc(all_common$gene)

all_common <- all_common %>%
  mutate(gene_key = sub("\\..*$", "", gene)) %>%
  left_join(ann, by = c("gene_key" = "ensembl_gene_id")) %>%
  mutate(gene_label = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, NA_character_))

# ===================== TOP GENES PER PANEL (top 1% by |Δ percentile|) ===================== #
top_tbl <- all_common %>%
  group_by(network, centrality) %>%
  mutate(abs_dperc = abs(delta_perc)) %>%
  slice_max(order_by = abs_dperc, prop = TOP_PROP, with_ties = FALSE) %>%
  ungroup()

# ===================== LOAD SPEARMAN RESULTS ===================== #
spearman_all <- purrr::map_dfr(networks, read_spearman_tbl)
if (nrow(spearman_all) == 0) stop("No spearman/<NET>_spearman.tsv files were loaded. Check paths and filenames.")

# ===================== RHO LABEL POSITION (TOP-RIGHT, robust under transforms + free scales) ===================== #
rho_pos <- all_common %>%
  group_by(network, centrality) %>%
  summarise(
    x_raw_min = suppressWarnings(min(AD,  na.rm = TRUE)),
    x_raw_max = suppressWarnings(max(AD,  na.rm = TRUE)),
    y_raw_min = suppressWarnings(min(CTL, na.rm = TRUE)),
    y_raw_max = suppressWarnings(max(CTL, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(spearman_all, by = c("network", "centrality")) %>%
  mutate(
    x_raw_min = ifelse(is.finite(x_raw_min), x_raw_min, 0),
    x_raw_max = ifelse(is.finite(x_raw_max), x_raw_max, x_raw_min),
    y_raw_min = ifelse(is.finite(y_raw_min), y_raw_min, 0),
    y_raw_max = ifelse(is.finite(y_raw_max), y_raw_max, y_raw_min),
    
    x_raw_max = ifelse(x_raw_max == x_raw_min, x_raw_max + 1e-12, x_raw_max),
    y_raw_max = ifelse(y_raw_max == y_raw_min, y_raw_max + 1e-12, y_raw_max),
    
    rho_abs   = abs(spearman_rho),
    rho_label = sprintf("\u03c1 = %.3f", spearman_rho),
    
    x_tr_min = plog$transform(x_raw_min),
    x_tr_max = plog$transform(x_raw_max),
    y_tr_min = plog$transform(y_raw_min),
    y_tr_max = plog$transform(y_raw_max),
    
    x_tr = x_tr_min + 0.78 * (x_tr_max - x_tr_min),
    y_tr = y_tr_min + 0.90 * (y_tr_max - y_tr_min),
    
    x = plog$inverse(x_tr),
    y = plog$inverse(y_tr)
  )

# ===================== PLOT A: HEATMAP of Spearman rho ===================== #
p_rho_tile <- ggplot(spearman_all, aes(x = network, y = centrality, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", spearman_rho)), size = 4) +
  scale_y_discrete(
    limits = rev(centralities),
    labels = unname(CENT_LABELS[rev(centralities)])
  ) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-1, 1), name = expression(rho)
  ) +
  labs(
    title = "Rank correlation of gene centrality between AD and control networks (Spearman \u03c1)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_rho_tile)
ggsave(file.path(OUT_DIR, "SpearmanRho_tile.pdf"), p_rho_tile, width = 8.5, height = 3.2)
ggsave(file.path(OUT_DIR, "SpearmanRho_tile.png"), p_rho_tile, width = 8.5, height = 3.2, dpi = 300)

# ===================== PLOT B: ONE GRID PER CENTRALITY ===================== #
make_scatter_grid_for_cent <- function(cent_name) {
  
  df_cent  <- all_common %>% filter(centrality == cent_name)
  top_cent <- top_tbl     %>% filter(centrality == cent_name)
  rho_cent <- rho_pos     %>% filter(centrality == cent_name)
  
  top_cent_lab <- top_cent %>% filter(!is.na(gene_label) & gene_label != "")
  cent_title <- cent_pretty(cent_name)
  
  ggplot(df_cent, aes(x = AD, y = CTL)) +
    geom_point(aes(color = HigherNetwork), alpha = 0.22, size = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey45") +
    geom_point(data = top_cent, aes(color = HigherNetwork), alpha = 0.95, size = 2.8) +
    
    geom_text_repel(
      data = top_cent_lab,
      aes(label = gene_label),
      color = "black",
      fontface = "bold",
      size = 3.3,
      box.padding = 0.55,
      point.padding = 0.35,
      segment.alpha = 0.55,
      segment.size = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 1.8,
      force_pull = 0.25,
      seed = 1
    ) +
    
    geom_label(
      data = rho_cent,
      aes(x = x, y = y, label = rho_label, fill = rho_abs),
      inherit.aes = FALSE,
      label.size = 0.25,
      label.r = unit(0.20, "lines"),
      label.padding = unit(0.28, "lines"),
      alpha = 0.98,
      size = 3.9,
      color = "black"
    ) +
    
    facet_wrap(~ network, scales = "free", nrow = 2) +
    scale_color_manual(values = c("AD" = "#E41A1C", "Control" = "#377EB8")) +
    
    scale_fill_gradient(
      limits = c(0, 1),
      low  = "#C9A227",
      high = "#FFF7CC",
      guide = "none"
    ) +
    
    labs(
      title = paste0(cent_title, " centrality: AD vs control"),
      subtitle = paste0("Top ", TOP_PROP * 100, "% genes by |Δ percentile rank|"),
      x = paste0(cent_title, " (AD network)"),
      y = paste0(cent_title, " (Control network)"),
      color = "Higher\nnetwork"
    ) +
    
    scale_x_continuous(
      trans  = plog,
      breaks = safe_breaks(n = 2),
      labels = fmt_si
    ) +
    scale_y_continuous(
      trans  = plog,
      breaks = safe_breaks(n = 2),
      labels = fmt_si
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10.8, hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(margin = margin(t = 4)),
      axis.text.y = element_text(margin = margin(r = 4))
    )
}

for (cent in centralities) {
  p <- make_scatter_grid_for_cent(cent)
  print(p)
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("Scatter_", cent, "_Top", TOP_PROP * 100, "pct_by_absDeltaPercentile.pdf")),
    plot = p, width = 16, height = 8
  )
  ggsave(
    filename = file.path(OUT_DIR, paste0("Scatter_", cent, "_Top", TOP_PROP * 100, "pct_by_absDeltaPercentile.png")),
    plot = p, width = 16, height = 8, dpi = 300
  )
}

cat("\nDONE\n")
