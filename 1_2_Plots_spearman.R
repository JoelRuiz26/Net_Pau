#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(biomaRt)
  library(grid)   # unit()
})

options(stringsAsFactors = FALSE)

# ============================================================
# CONFIG
# ============================================================
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"

# DGE limma
dge_dir <- file.path(base_dir, "2_1_DGE_limma")

# Degree delta percentile (COMMON tables)
deg_dir <- file.path(base_dir, "1_1_Centrality_comparison_rho", "degree")

# Output
out_dir <- file.path(base_dir, "Figures_Grids", "Quadrants_INTERSECTION_only")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Thresholds (match your reporting)
FDR_CUT <- 0.05
LFC_CUT <- 0.5

# Top 10% by |delta_perc| per region
TOP_PCT <- 0.10

# Labels
LABEL_TOP_N <- 12  # per region, within hit-only
USE_HGNC <- TRUE

# Saving
out_png <- file.path(out_dir, "Quadrants_LFC_vs_DeltaDegreePercentile_INTERSECTION.png")
out_pdf <- file.path(out_dir, "Quadrants_LFC_vs_DeltaDegreePercentile_INTERSECTION.pdf")
out_tsv <- file.path(out_dir, "INTERSECTION_table.tsv")

# HGNC cache
map_cache <- file.path(out_dir, "ensembl_to_hgnc_cache.rds")

# Region order
REG_ORDER <- c("PCC","HCN","DLPFC","TC","CRB")

# ============================================================
# HELPERS
# ============================================================
read_tsv_quiet <- function(path) readr::read_tsv(path, col_types = cols(), progress = FALSE)

norm_key <- function(x) {
  x %>%
    str_replace(regex("^ROSMAP_", ignore_case = TRUE), "") %>%
    str_replace(regex("^MAYO_",   ignore_case = TRUE), "") %>%
    toupper()
}

top_abs_cut <- function(x, top_pct = 0.10) {
  stats::quantile(abs(x), probs = 1 - top_pct, na.rm = TRUE, names = FALSE, type = 7)
}

map_ensembl_to_hgnc <- function(ens_ids) {
  ens_ids <- unique(sub("\\..*$", "", ens_ids))
  
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

# ============================================================
# 1) DISCOVER FILES + PAIR BY NORMALIZED REGION KEY
# ============================================================
dge_files <- list.files(dge_dir, pattern = "_limma_AD_vs_Control_allGenes\\.tsv$", full.names = TRUE)
deg_files <- list.files(deg_dir, pattern = "_degree_COMMON\\.tsv$", full.names = TRUE)

if (length(dge_files) == 0) stop("No limma TSVs found in: ", dge_dir)
if (length(deg_files) == 0) stop("No degree COMMON TSVs found in: ", deg_dir)

dge_tbl <- tibble(
  dge_path = dge_files,
  dge_key_raw = basename(dge_files) %>% str_replace("_limma_AD_vs_Control_allGenes\\.tsv$", ""),
  join_key = norm_key(dge_key_raw)
)

deg_tbl <- tibble(
  deg_path = deg_files,
  deg_key_raw = basename(deg_files) %>% str_replace("_degree_COMMON\\.tsv$", ""),
  join_key = norm_key(deg_key_raw)
)

pairs <- inner_join(dge_tbl, deg_tbl, by = "join_key") %>% arrange(join_key)

if (nrow(pairs) == 0) {
  stop(
    "No matching region keys between DGE and degree COMMON after normalization.\n",
    "DGE join keys: ", paste(head(unique(dge_tbl$join_key), 10), collapse = ", "), "\n",
    "DEG join keys: ", paste(head(unique(deg_tbl$join_key), 10), collapse = ", ")
  )
}

message("Matched regions: ", paste(pairs$join_key, collapse = ", "))
print(pairs %>% dplyr::select(join_key, dge_key_raw, deg_key_raw))

# ============================================================
# 2) LOAD + MERGE (per region)
# ============================================================
merged_list <- pmap(
  list(pairs$dge_path, pairs$deg_path, pairs$join_key),
  function(dge_path, deg_path, join_key) {
    
    dge <- read_tsv_quiet(dge_path) %>%
      transmute(
        gene,
        gene_clean = sub("\\..*$", "", gene),
        logFC      = as.numeric(logFC),
        adj.P.Val  = as.numeric(adj.P.Val)
      ) %>%
      mutate(
        pass_dge = (adj.P.Val <= FDR_CUT & abs(logFC) >= LFC_CUT)
      )
    
    deg <- read_tsv_quiet(deg_path) %>%
      transmute(
        gene,
        gene_clean = sub("\\..*$", "", gene),
        delta_perc = as.numeric(delta_perc)
      )
    
    # merge by gene (exact). If your ids differ by versions, join by gene_clean.
    df <- inner_join(dge, deg, by = "gene") %>%
      mutate(region = join_key)
    
    # if inner_join by gene yields too few due to version mismatch, fallback to gene_clean join:
    if (nrow(df) < 100) {
      df <- inner_join(dge, deg, by = "gene_clean", suffix = c("_dge", "_deg")) %>%
        transmute(
          gene = gene_dge,
          gene_clean,
          logFC,
          adj.P.Val,
          pass_dge,
          delta_perc,
          region = join_key
        )
    }
    
    # per-region cutoff for top10% |delta_perc|
    cut <- top_abs_cut(df$delta_perc, top_pct = TOP_PCT)
    
    df %>%
      mutate(
        delta_abs_cut = cut,
        top_delta10   = abs(delta_perc) >= cut,
        hit           = pass_dge & top_delta10
      )
  }
)

all_df <- bind_rows(merged_list) %>%
  mutate(region = factor(region, levels = REG_ORDER))

# sanity check
message("\n[CHECK] Intersection counts per region:")
print(all_df %>%
        group_by(region) %>%
        summarise(
          n_merged = n(),
          n_pass_dge = sum(pass_dge, na.rm = TRUE),
          n_top10 = sum(top_delta10, na.rm = TRUE),
          n_hit = sum(hit, na.rm = TRUE),
          .groups = "drop"
        ) %>% arrange(desc(n_hit)))

hit_df <- all_df %>% filter(hit)

if (nrow(hit_df) == 0) {
  stop("Intersection is empty (n_hit=0). Check thresholds or whether degree COMMON and limma genes match.")
}

# ============================================================
# 3) HGNC ANNOTATION (only for hit_df)
# ============================================================
if (USE_HGNC) {
  if (file.exists(map_cache)) {
    message("\nLoading cached HGNC mapping: ", map_cache)
    map_tbl <- readRDS(map_cache)
  } else {
    message("\nQuerying biomaRt for HGNC symbols (hit genes only)...")
    map_tbl <- map_ensembl_to_hgnc(hit_df$gene_clean)
    saveRDS(map_tbl, map_cache)
    message("Saved mapping cache: ", map_cache)
  }
  
  hit_df <- hit_df %>%
    left_join(map_tbl, by = c("gene_clean" = "ensembl_gene_id")) %>%
    mutate(label = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, gene_clean))
} else {
  hit_df <- hit_df %>% mutate(label = gene_clean)
}

# ============================================================
# 4) PICK LABELS (top by |delta_perc| then adj.P.Val)
# ============================================================
label_df <- hit_df %>%
  group_by(region) %>%
  arrange(desc(abs(delta_perc)), adj.P.Val) %>%
  slice_head(n = LABEL_TOP_N) %>%
  ungroup()

# ============================================================
# 5) PLOT (INTERSECTION ONLY) — nice colors + quadrants
# ============================================================
# palette
COL_POS <- "#D1495B"   # red-ish
COL_NEG <- "#00798C"   # blue/teal
COL_LAB <- "grey10"

hit_df <- hit_df %>%
  mutate(
    quad = case_when(
      logFC >= 0 & delta_perc >= 0 ~ "Q1 (+LFC, +Δperc)",
      logFC <  0 & delta_perc >= 0 ~ "Q2 (-LFC, +Δperc)",
      logFC <  0 & delta_perc <  0 ~ "Q3 (-LFC, -Δperc)",
      TRUE                         ~ "Q4 (+LFC, -Δperc)"
    ),
    sign_delta = ifelse(delta_perc >= 0, "Δperc ≥ 0", "Δperc < 0")
  )

# per region cutoff lines (top10% cutoff differs by region)
cuts_df <- hit_df %>% distinct(region, delta_abs_cut)

p <- ggplot(hit_df, aes(x = logFC, y = delta_perc)) +
  
  # quadrant lines
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey35") +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "grey35") +
  
  # threshold lines (optional but informative)
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.40, color = "grey45") +
  geom_hline(
    data = cuts_df, aes(yintercept =  delta_abs_cut),
    inherit.aes = FALSE, linetype = "dashed", linewidth = 0.40, color = "grey45"
  ) +
  geom_hline(
    data = cuts_df, aes(yintercept = -delta_abs_cut),
    inherit.aes = FALSE, linetype = "dashed", linewidth = 0.40, color = "grey45"
  ) +
  
  # points (ONLY intersection)
  geom_point(aes(color = sign_delta), size = 2.2, alpha = 0.85) +
  
  # emphasize labeled points
  geom_point(data = label_df, aes(color = sign_delta), size = 3.0, alpha = 0.98) +
  
  # labels
  ggrepel::geom_label_repel(
    data = label_df,
    aes(label = label),
    fill = "white",
    label.size = 0.25,
    label.r = unit(0.15, "lines"),
    label.padding = unit(0.18, "lines"),
    color = COL_LAB,
    size = 3.6,
    box.padding = 0.55,
    point.padding = 0.35,
    min.segment.length = 0,
    segment.alpha = 0.55,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 1
  ) +
  
  facet_wrap(~ region, ncol = 3, scales = "free_y") +
  
  scale_color_manual(
    values = c("Δperc ≥ 0" = COL_POS, "Δperc < 0" = COL_NEG),
    name = "Δ degree percentile"
  ) +
  
  labs(
    title = "Intersection-only: significant DGE ∩ top10% |Δ degree percentile|",
    subtitle = sprintf(
      "DGE: FDR ≤ %.2g and |log2FC| ≥ %.2g;  Δperc cutoff: top %.0f%% by |Δperc| within each region. Points shown = intersection only.",
      FDR_CUT, LFC_CUT, 100*TOP_PCT
    ),
    x = "log2 fold-change (AD vs Control)",
    y = "Δ degree percentile (perc_AD − perc_Control)"
  ) +
  
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 20)
  )

# ============================================================
# 6) SAVE
# ============================================================
write_tsv(hit_df, out_tsv)

ggsave(out_png, p, width = 13.5, height = 8.5, dpi = 600)
ggsave(out_pdf, p, width = 13.5, height = 8.5)

message("\nSaved:")
message(" - ", out_png)
message(" - ", out_pdf)
message(" - ", out_tsv)
message("DONE ✅")
