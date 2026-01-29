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
in_dir   <- file.path(base_dir, "2_1_DGE_limma")
out_dir  <- file.path(base_dir, "2_1_DGE_limma/2_2_DGE_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR_CUT <- 0.05
LFC_CUT <- 0.5
TOP_N_PER_SIDE <- 8

USE_FDR_ON_X <- TRUE  # TRUE: -log10(adj.P.Val); FALSE: -log10(P.Value)
P_FLOOR <- 1e-300

out_png <- file.path(out_dir, "VolcanoRow_byRegion_AD_vs_Control.png")
out_pdf <- file.path(out_dir, "VolcanoRow_byRegion_AD_vs_Control.pdf")

# ============================================================
# HELPERS
# ============================================================
read_limma_tsv <- function(path) {
  df <- read_tsv(path, col_types = cols(), progress = FALSE)
  needed <- c("gene","logFC","P.Value","adj.P.Val")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", basename(path), ": ", paste(miss, collapse = ", "))
  
  region_raw <- basename(path) %>%
    str_replace("_limma_AD_vs_Control_allGenes\\.tsv$", "")
  
  df %>% mutate(region_raw = region_raw)
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
# 1) LOAD
# ============================================================
files <- list.files(in_dir, pattern = "_limma_AD_vs_Control_allGenes\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No limma TSVs found in: ", in_dir)

message("Found ", length(files), " region files:")
print(basename(files))

all_df <- map_dfr(files, read_limma_tsv)

region_map <- tibble(
  region_raw = c("Mayo_CRB","Mayo_TC","ROSMAP_DLPFC","ROSMAP_HCN","ROSMAP_PCC"),
  region     = c("CRB","TC","DLPFC","HCN","PCC")
)

all_df <- all_df %>%
  left_join(region_map, by = "region_raw") %>%
  mutate(
    region = ifelse(is.na(region), region_raw, region),
    region = factor(region, levels = c("CRB","TC","DLPFC","HCN","PCC"))
  )

message("\n[CHECK] Genes per region:")
print(all_df %>% count(region, name = "n_genes") %>% arrange(desc(n_genes)))

# ============================================================
# 2) AXES + SIGNIFICANCE
# ============================================================
all_df <- all_df %>%
  mutate(
    gene_clean = sub("\\..*$", "", gene),
    
    p_for_x    = if (USE_FDR_ON_X) adj.P.Val else P.Value,
    p_for_x    = pmax(p_for_x, P_FLOOR),
    x_sig      = -log10(p_for_x),
    y_fc       = as.numeric(logFC),
    
    pass_fdr   = adj.P.Val <= FDR_CUT,
    pass_lfc   = abs(logFC) >= LFC_CUT,
    pass_both  = pass_fdr & pass_lfc,
    
    # ✅ professional labels
    direction = case_when(
      pass_both & logFC > 0 ~ "Overexpressed",
      pass_both & logFC < 0 ~ "Underexpressed",
      TRUE ~ "Not significant"
    )
  )

message("\n[CHECK] Significant counts per region (FDR & |logFC|):")
print(all_df %>%
        group_by(region) %>%
        summarise(
          n_genes = n(),
          n_BOTH  = sum(pass_both),
          .groups = "drop"
        ) %>% arrange(desc(n_BOTH)))

# ============================================================
# 3) ANNOTATION (HGNC only)
# ============================================================
map_cache <- file.path(out_dir, "ensembl_to_hgnc_cache.rds")

if (file.exists(map_cache)) {
  message("\nLoading cached HGNC mapping: ", map_cache)
  map_tbl <- readRDS(map_cache)
} else {
  message("\nQuerying biomaRt for HGNC symbols...")
  map_tbl <- map_ensembl_to_hgnc(all_df$gene_clean)
  saveRDS(map_tbl, map_cache)
  message("Saved mapping cache: ", map_cache)
}

all_df <- all_df %>%
  left_join(map_tbl, by = c("gene_clean" = "ensembl_gene_id")) %>%
  mutate(
    has_hgnc = !is.na(hgnc_symbol),
    label    = hgnc_symbol
  )

# ============================================================
# 4) LABELS: TOP N UP + TOP N DOWN PER REGION (HGNC-only)
# ============================================================
label_df <- all_df %>%
  filter(pass_both, has_hgnc) %>%
  group_by(region) %>%
  group_modify(~{
    up <- .x %>%
      filter(logFC > 0) %>%
      arrange(adj.P.Val, desc(logFC)) %>%
      slice_head(n = TOP_N_PER_SIDE)
    
    down <- .x %>%
      filter(logFC < 0) %>%
      arrange(adj.P.Val, logFC) %>%
      slice_head(n = TOP_N_PER_SIDE)
    
    bind_rows(up, down)
  }) %>%
  ungroup()

message("\n[CHECK] Labels per region (HGNC-only):")
print(label_df %>% count(region, name = "n_labels") %>% arrange(desc(n_labels)))

# ✅ needed because you use label_over / label_under in the plot
label_over  <- label_df %>% filter(direction == "Overexpressed")
label_under <- label_df %>% filter(direction == "Underexpressed")

# ============================================================
# 5) PLOT: ONE ROW, SHARED Y AXIS (clean) + colored segments
# ============================================================
COL_OVER  <- "#D7301F"
COL_UNDER <- "#2C7FB8"
COL_NS    <- "grey75"

x_cut_line <- -log10(pmax(FDR_CUT, P_FLOOR))

p <- ggplot(all_df, aes(x = x_sig, y = y_fc)) +
  # thresholds
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey35") +
  geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.35, color = "grey35") +
  geom_vline(xintercept = x_cut_line, linetype = "dashed", linewidth = 0.35, color = "grey35") +
  
  # points: NS lighter; sig visible
  geom_point(
    data = all_df %>% filter(direction == "Not significant"),
    aes(color = direction),
    size = 0.7, alpha = 0.22
  ) +
  geom_point(
    data = all_df %>% filter(direction != "Not significant"),
    aes(color = direction),
    size = 0.85, alpha = 0.50
  ) +
  
  # anchor points for labeled genes
  geom_point(
    data = label_df,
    aes(color = direction),
    size = 1.6, alpha = 0.90
  ) +
  
  # labels Overexpressed (red segments)
  ggrepel::geom_label_repel(
    data = label_over,
    aes(label = label),
    color = "grey10",
    segment.color = COL_OVER,
    show.legend = FALSE,
    size = 3.05,
    fill = NA,
    label.size = 0,
    label.padding = grid::unit(0.12, "lines"),
    direction = "both",
    force = 3.5,
    force_pull = 0.35,
    box.padding = 0.35,
    point.padding = 0.45,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.75,
    segment.size = 0.35,
    seed = 2
  ) +
  
  # labels Underexpressed (blue segments)
  ggrepel::geom_label_repel(
    data = label_under,
    aes(label = label),
    color = "grey10",
    segment.color = COL_UNDER,
    show.legend = FALSE,
    size = 3.05,
    fill = NA,
    label.size = 0,
    label.padding = grid::unit(0.12, "lines"),
    direction = "both",
    force = 3.5,
    force_pull = 0.35,
    box.padding = 0.35,
    point.padding = 0.45,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.75,
    segment.size = 0.35,
    seed = 2
  ) +
  
  facet_grid(. ~ region) +
  
  # legend order: Overexpressed, Not significant, Underexpressed
  scale_color_manual(
    values = c(
      "Overexpressed"   = COL_OVER,
      "Not significant" = COL_NS,
      "Underexpressed"  = COL_UNDER
    ),
    breaks = c("Overexpressed", "Not significant", "Underexpressed")
  ) +
  
  guides(color = guide_legend(override.aes = list(shape = 16, size = 3, alpha = 1))) +
  
  labs(
    title = "Differential expression across brain regions",
    x = expression(-log[10](FDR)),
    y = expression(log[2]~FC),
    color = "Expression"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text.x = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ============================================================
# 6) SAVE
# ============================================================
ggsave(out_png, p, width = 12.8, height = 6.0, dpi = 600)
ggsave(out_pdf, p, width = 12.8, height = 6.0)
