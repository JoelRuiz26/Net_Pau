#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(scales)   
  library(ggsci)  
})

#install.packages("ggsci")

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"
in_dir   <- file.path(base_dir, "4_DGE_centrality")
out_dir  <- file.path(in_dir, "Quadrants_top10pct_LFC_vs_DeltaPerc_DGEonly_ONEPLOT")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TOP_PCT <- 0.10
LABEL_TOP_N <- 10
REG_ORDER <- c("PCC","HCN","DLPFC","TC","CRB")

# shapes per region (distinct)
REG_SHAPES <- c(
  "PCC"   = 15,  # square
  "HCN"   = 17,  # triangle
  "DLPFC" = 16,  # circle
  "TC"    = 18,  # diamond
  "CRB"   = 8    # star
)

# --------
lan <- ggsci::pal_lancet("lanonc")(9)
REG_COLORS <- c(
  "PCC"   = lan[1],
  "HCN"   = lan[2],
  "DLPFC" = lan[3],
  "TC"    = lan[4],
  "CRB"   = lan[5]
)

# ---- X-axis "cut" (compress gap) around LFC thresholds ----
gap_trans <- function(left = -0.5, right = 0.5, compress = 0.25) {
  stopifnot(left < right, compress > 0, compress <= 1)
  
  forward <- function(x) {
    x <- as.numeric(x)
    out <- x
    in_gap <- !is.na(x) & x >= left & x <= right
    right_side <- !is.na(x) & x > right
    out[in_gap] <- left + (x[in_gap] - left) * compress
    out[right_side] <- left + (right - left) * compress + (x[right_side] - right)
    out
  }
  
  inverse <- function(y) {
    y <- as.numeric(y)
    out <- y
    gap_right_t <- left + (right - left) * compress
    in_gap_t <- !is.na(y) & y >= left & y <= gap_right_t
    right_side_t <- !is.na(y) & y > gap_right_t
    out[in_gap_t] <- left + (y[in_gap_t] - left) / compress
    out[right_side_t] <- right + (y[right_side_t] - gap_right_t)
    out
  }
  
  scales::trans_new(
    name      = sprintf("gap_%g_%g_%g", left, right, compress),
    transform = forward,
    inverse   = inverse
  )
}

# ===================== LOAD ===================== #
paths <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(paths) == 0) stop("No TSV files found in: ", in_dir)

deg_list2 <- set_names(
  map(paths, ~ read_tsv(.x, col_types = cols(), progress = FALSE)),
  tools::file_path_sans_ext(basename(paths))
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

# ===================== BUILD DATA (DGE only; keep ALL percentiles) ===================== #
thr <- 1 - TOP_PCT

df_all <- imap_dfr(deg_list2, function(df, nm) {
  region <- get_region(nm)
  if (is.na(region)) return(NULL)
  
  needed <- c("delta_abs_percentil", "delta_perc", "logFC", "label", "gene_clean", "is_dge")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", nm, ": ", paste(miss, collapse = ", "))
  
  df %>%
    mutate(region = region) %>%
    filter(is_dge) %>%
    mutate(
      is_top10 = delta_abs_percentil >= thr
    )
}) %>%
  mutate(
    region = factor(region, levels = REG_ORDER),
    row_id = row_number()  # <- NEW: id estable para separar etiquetados vs no-etiquetados
  )

if (nrow(df_all) == 0) stop("No rows after filtering to DGE-only.")

# ===================== TOP10 DATA + LABELS ===================== #
df_top <- df_all %>% filter(is_top10)
if (nrow(df_top) == 0) stop("No rows in top10% among DGE genes. (Revisa delta_abs_percentil).")

# NEW: score combinado para priorizar BOTH |Δperc| y |LFC|
df_top <- df_top %>%
  mutate(score_combo = abs(delta_perc) * abs(logFC))

# Labels = TOP 10 por región según score_combo (no solo degree)
label_df <- df_top %>%
  filter(!is.na(label) & label != "") %>%
  group_by(region) %>%
  arrange(desc(score_combo), desc(abs(delta_perc)), desc(abs(logFC))) %>%
  slice_head(n = LABEL_TOP_N) %>%
  ungroup()

# NEW: top10% pero SIN etiqueta (para hacerlos un poco más visibles)
df_top_nolabel <- df_top %>%
  filter(!(row_id %in% label_df$row_id))

# ===================== helpers for break marks + breaks ===================== #
y_rng <- range(df_all$delta_perc, na.rm = TRUE)
y_span <- diff(y_rng)
y_break <- y_rng[1] - 0.06 * y_span

x_rng <- range(df_all$logFC, na.rm = TRUE)
x_breaks <- seq(floor(x_rng[1]), ceiling(x_rng[2]), by = 1)

# ===================== PLOT ===================== #
p <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey35") +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "grey35") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.35, color = "grey55") +
  
  # non-top10 DGE genes (gris; con formas por región)
  geom_point(
    data = df_all %>% filter(!is_top10),
    aes(x = logFC, y = delta_perc, shape = region),
    color = "grey80",
    alpha = 0.55,
    size = 2.0,
    show.legend = FALSE
  ) +
  
  # top10% DGE genes SIN etiqueta (más visibles)
  geom_point(
    data = df_top_nolabel,
    aes(x = logFC, y = delta_perc, color = region, shape = region),
    alpha = 0.95,     
    size  = 2.8       
  ) +
  
  # top10% DGE genes CON etiqueta (igual o apenas más)
  geom_point(
    data = label_df,
    aes(x = logFC, y = delta_perc, color = region, shape = region),
    alpha = 1.0,
    size  = 3.0
  ) +
  
  # labels (solo TOP 10 por score_combo)
  ggrepel::geom_label_repel(
    data = label_df,
    aes(x = logFC, y = delta_perc, label = label, color = region),
    fill = "white",
    label.size = 0.25,
    size = 3.4,
    box.padding = 0.55,
    point.padding = 0.35,
    min.segment.length = 0,
    segment.alpha = 0.6,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 1,
    show.legend = FALSE
  ) +
  
  # marcas del “corte” visual
  annotate("segment", x = -0.5, xend = -0.47, y = y_break, yend = y_break + 0.03*y_span, linewidth = 0.5, color = "grey35") +
  annotate("segment", x = -0.47, xend = -0.44, y = y_break, yend = y_break + 0.03*y_span, linewidth = 0.5, color = "grey35") +
  annotate("segment", x =  0.44, xend =  0.47, y = y_break, yend = y_break + 0.03*y_span, linewidth = 0.5, color = "grey35") +
  annotate("segment", x =  0.47, xend =  0.50, y = y_break, yend = y_break + 0.03*y_span, linewidth = 0.5, color = "grey35") +
  
  scale_shape_manual(values = REG_SHAPES, guide = "none") +
  scale_color_manual(values = REG_COLORS, breaks = REG_ORDER, name = "region") +
  
  scale_x_continuous(
    trans  = gap_trans(left = -0.5, right = 0.5, compress = 0.1),
    breaks = x_breaks,
    labels = scales::label_number(accuracy = 1)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 9)) +
  
  labs(
    title = "Integrated expression and network centrality changes in Alzheimer’s disease",
    #    subtitle = "Highlighted genes show the strongest combined shifts in expression and connectivity.",
    x = "Expression change (log2 fold-change)",
    y = "Connectivity change (Delta degree)"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = REG_SHAPES[levels(df_all$region)],
        size  = 3,
        alpha = 1
      )
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(10, 20, 25, 15)
  )

print(p)

# ===================== SAVE  ===================== #

ggsave(file.path(out_dir, "4_2_1_Top10pct_LFC_vs_DeltaPerc.png"),
       p, width = 12.5, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "4_2_1__Top10pct_LFC_vs_DeltaPerc.pdf"),
       p, width = 12.5, height = 7.5)

saveRDS(deg_list2, file = file.path(out_dir, "4_2_1_deg_list2_all_regions.rds"))
saveRDS(label_df, file = file.path(out_dir, "4_2_1_top10_genes.rds"))

