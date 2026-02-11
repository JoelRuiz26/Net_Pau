#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(tibble)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"
out_dir  <- file.path(base_dir, "4_DGE_centrality")  # donde guardaste las tablas nuevas

# Percentiles a evaluar: 1% ... 50% (top X% por |Δperc|)
PCTS <- 1:50

# ===================== LOAD ===================== #
deg_paths <- list.files(out_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(deg_paths) == 0) stop("No TSV files found in: ", out_dir)

deg_list2 <- set_names(
  map(deg_paths, ~ read_tsv(.x, col_types = cols(), progress = FALSE)),
  tools::file_path_sans_ext(basename(deg_paths))
)

# ===================== REGION PARSING ===================== #
get_region <- function(name_or_path) {
  x <- toupper(basename(name_or_path))
  case_when(
    str_detect(x, "DLPFC") ~ "DLPFC",
    str_detect(x, "HCN")   ~ "HCN",
    str_detect(x, "PCC")   ~ "PCC",
    str_detect(x, "TC")    ~ "TC",
    str_detect(x, "CRB")   ~ "CRB",
    TRUE ~ NA_character_
  )
}

# ===================== BUILD CURVES ===================== #
# Para cada región y cada percentil p (1..50):
# threshold = 1 - p/100  (top p% => delta_abs_percentil >= threshold)
# n_dge_at_cut = n_dge al final de ese subconjunto (equivale a cuántos DGE entran en ese corte)
curve_df <- imap_dfr(
  deg_list2,
  function(df, nm) {
    
    if (!all(c("delta_abs_percentil", "n_dge") %in% colnames(df))) {
      stop("Missing required columns in: ", nm, " (need delta_abs_percentil and n_dge)")
    }
    
    region <- get_region(nm)
    
    # Asegura orden (por si acaso)
    df <- df %>% arrange(desc(delta_abs_percentil), desc(abs(delta_perc)))
    
    map_dfr(PCTS, function(p) {
      thr <- 1 - p/100  # top p%
      sub <- df %>% filter(delta_abs_percentil >= thr)
      
      tibble(
        region   = region,
        pct_top  = p,
        threshold = thr,
        n_dge    = if (nrow(sub) == 0) 0L else as.integer(tail(sub$n_dge, 1))
      )
    })
  }
) %>%
  filter(!is.na(region)) %>%
  mutate(region = factor(region, levels = c("PCC","HCN","DLPFC","TC","CRB")))

# ===================== PLOT ===================== #
p <- ggplot(curve_df, aes(x = pct_top, y = n_dge)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  facet_wrap(~ region, scales = "free_y") +
  labs(
    title = "Cumulated DGE vs |Δ degree percentile|",
    subtitle = "top percentil (1% a 50%)",
    x = "|Δperc|",
    y = "n_dge"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

print(p)

ggsave(file.path(out_dir, "n_dge_vs_percentil_top1_50.png"), p, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "n_dge_vs_percentil_top1_50.pdf"), p, width = 12, height = 7)
