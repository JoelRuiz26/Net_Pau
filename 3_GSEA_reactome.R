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
  slice_head(n = 5) %>%
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

plot_df <- top5_up_down_by_region %>%
  group_by(term_simple) %>%
  mutate(absNES_max = max(abs(NES), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    region      = factor(region, levels = REG_ORDER),
    sig         = -log10(p.adjust),
    term_simple = reorder(term_simple, absNES_max)
  )

p <- ggplot(plot_df, aes(region, term_simple)) +
  geom_point(aes(size = sig, fill = NES), shape = 21, color = "grey25") +
  scale_size(
    range = c(0, 8),
    name  = "-log10(FDR)",
    guide = guide_legend(
      override.aes = list(fill = "black", color = "black")
    )
  ) +
  scale_fill_gradient2(name = "NES") +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    axis.text.y  = element_text(size = 8, face = "bold"),
    axis.text.x  = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "3_3_Bubble_terms.png"),
  plot = p,
  width = 8.5,
  height = 10,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "3_3_Bubble_terms.pdf"),
  plot = p,
  width = 8.5,
  height = 10
)



