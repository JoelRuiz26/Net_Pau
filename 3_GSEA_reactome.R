#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(clusterProfiler)
  library(msigdbr)
  library(ggplot2)
})

base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"
in_dir   <- file.path(base_dir, "2_1_DGE_limma")
out_dir  <- file.path(base_dir, "3_GSEA_REACTOME")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pattern_in <- "_limma_AD_vs_Control_allGenes\\.tsv$"

TERM2GENE_REACT <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME") |>
  transmute(term = gs_name, gene = db_ensembl_gene) |>
  filter(!is.na(gene), gene != "") |>
  distinct()

run_one <- function(path){
  region <- str_replace(basename(path), pattern_in, "")
  df <- read_tsv(path, col_types = cols(), progress = FALSE) |>
    transmute(gene = sub("\\..*$","", gene), stat = as.numeric(t)) |>
    filter(!is.na(gene), gene != "", is.finite(stat)) |>
    group_by(gene) |>
    summarise(stat = stat[which.max(abs(stat))], .groups="drop")
  
  geneList <- df$stat; names(geneList) <- df$gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  as_tibble(as.data.frame(
    suppressWarnings(GSEA(
      geneList     = geneList,
      TERM2GENE    = TERM2GENE_REACT,
      pvalueCutoff = 1,
      minGSSize    = 10,
      maxGSSize    = 500,
      eps          = 0,
      verbose      = FALSE
    ))
  )) |>
    mutate(region = region, .before = 1)
}

files <- list.files(in_dir, pattern = pattern_in, full.names = TRUE)
gsea_full <- map_dfr(files, run_one)
#[1] 6553

write_tsv(gsea_full, file.path(out_dir, "3_1_Full_terms_Region"))


#########################################################################
#Filter results

x <- gsea_full %>%
  mutate(
    region = str_remove(region, "^(Mayo_|ROSMAP_)"),
    term_simple = str_replace_all(str_remove(Description, "^REACTOME_"), "_", " ")
  ) %>%
  filter(p.adjust < 0.05)
#[1] 1218
length(unique(x$term_simple))
#[1] 867

top5_up_down_by_region <- x %>%
  mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
  group_by(region, direction) %>%
  arrange(
    region,
    direction,
    desc(abs(NES))   # prioriza términos más enriquecidos
  ) %>%
  slice_head(n = 4) %>%
  ungroup()

# guardar tabla
write_tsv(
  top5_up_down_by_region,
  file.path(out_dir, "3_2_Top20_terms.tsv")
)
length(unique(top5_up_down_by_region$term_simple))
#[1] 35   / 50


########################################################
#Plot

REG_ORDER <- c("CRB","TC","DLPFC","HCN","PCC")

plot_df2 <- top5_up_down_by_region %>%
  group_by(term_simple) %>%
  mutate(absNES_max = max(abs(NES), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    region      = factor(region, levels = REG_ORDER),
    sig         = -log10(p.adjust),
    term_simple = reorder(term_simple, absNES_max),
    
    # categorías de FDR (ajústalas si quieres)
    fdr_cat = dplyr::case_when(
      p.adjust <= 1e-3 ~ "≤1e-3",
      p.adjust <= 1e-2 ~ "≤1e-2",
      p.adjust <= 5e-2 ~ "≤5e-2",
      TRUE             ~ ">5e-2"
    ),
    
    # estrellas (clásico y MUY legible)
    star = dplyr::case_when(
      p.adjust <= 1e-3 ~ "***",
      p.adjust <= 1e-2 ~ "**",
      p.adjust <= 5e-2 ~ "*",
      TRUE             ~ ""
    ),
    
    # alpha continuo (suaviza no-significativos sin “tapar” el NES)
    alpha_sig = scales::rescale(sig, to = c(0.25, 1.0), from = range(sig, finite = TRUE))
  )

p_altA <- ggplot(plot_df2, aes(region, term_simple)) +
  geom_tile(
    aes(fill = NES, alpha = alpha_sig),
    color = "grey92", linewidth = 0.4
  ) +
  geom_text(
    aes(label = star),
    fontface = "bold", size = 4, color = "grey10"
  ) +
  scale_fill_gradient2(
    name = "NES",
    low  = "#2C7FB8",
    mid  = "#F7F7F7",
    high = "#D7301F",
    midpoint = 0
  ) +
  scale_alpha(guide = "none") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    
    # 👇 sutil y legible
    axis.text.y  = element_text(size = 10, face = "bold"),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold"),
    
    plot.caption = element_text(size = 10, face = "bold", hjust = 0),
    
    # ✅ "celdas grises" tipo volcano (strip gris #grey92)
    axis.text.y.right = element_text(face = "bold"),
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(face = "bold")
  ) +
  # ✅ esto crea el "strip" gris tipo volcano para la región
  facet_grid(. ~ region, switch = "x")

# Márgenes + guardar (menos agresivo que r = 60)
p2 <- p_altA + theme(plot.margin = margin(t = 8, r = 25, b = 8, l = 8))

save.image(file.path(out_dir, "3_4_Image_GSEA.RData"))
#load("/STORAGE/csbig/jruiz/Redes_Pau/3_GSEA_REACTOME/3_4_Image_GSEA.RData")

saveRDS(p2, file.path(out_dir, "3_3_Heatmap_Stars_alpha_plot.rds"))

ggsave(
  file.path(out_dir, "3_3_Heatmap_Stars_alpha.png"),
  p2, width = 16, height = 9, units = "in", dpi = 600, limitsize = FALSE)

ggsave(
  file.path(out_dir, "3_3_Heatmap_Stars_alpha.pdf"),
  p2, width = 16, height = 9, units = "in",
  limitsize = FALSE, device = cairo_pdf)
