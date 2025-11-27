###############################################
## DEGREE COMPARISON BETWEEN AD AND CTL NETWORKS
## - Compute degree in original graphs (all nodes)
## - Build global degree table
## - Restrict to common genes
## - Δdegree + Z-score for common genes
## - Plot all genes + top 10 by |Z-score| labeled
###############################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

## =========================
## 0) Paths and output dir
## =========================
OUT_DIR <- "~/Net_Pau/Degree_comparison/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## =========================
## 1) Load graphs (original, unfiltered)
## =========================
g_AD  <- readRDS("/STORAGE/csbig/jruiz/Redes_Pau/g_AD_200k.rds")
g_CTL <- readRDS("/STORAGE/csbig/jruiz/Redes_Pau/g_ctl_200k.rds")

## =========================
## 2) Degree in original graphs (all nodes)
## =========================
deg_AD_full  <- igraph::degree(g_AD)
deg_CTL_full <- igraph::degree(g_CTL)

AD_names  <- names(deg_AD_full)
CTL_names <- names(deg_CTL_full)

## =========================
## 3) Global degree table (all genes)
## =========================
deg_full_df <- full_join(
  tibble(gene = AD_names,  degree_AD  = as.numeric(deg_AD_full)),
  tibble(gene = CTL_names, degree_CTL = as.numeric(deg_CTL_full)),
  by = "gene"
)

saveRDS(deg_full_df, file.path(OUT_DIR, "Degree_AD_CTL_FULL_table.rds"))
write_tsv(deg_full_df, file.path(OUT_DIR, "Degree_AD_CTL_FULL_table.tsv"))

## =========================
## 4) Restrict to common genes (present in both networks)
## =========================
deg_common_df <- deg_full_df %>%
  filter(!is.na(degree_AD) & !is.na(degree_CTL))

## Optional: remove genes with degree 0 in both
deg_common_nofloat <- deg_common_df %>%
  filter(!(degree_AD == 0 & degree_CTL == 0))

## =========================
## 5) Δdegree + Z-score (only for common genes)
## =========================
deg_common_df <- deg_common_df %>%
  mutate(
    delta_degree     = degree_AD - degree_CTL,
    abs_delta_degree = abs(delta_degree),
    mean_degree      = (degree_AD + degree_CTL) / 2,
    z_change         = (delta_degree - mean(delta_degree)) / sd(delta_degree)
  )

deg_common_nofloat <- deg_common_df %>%
  filter(!(degree_AD == 0 & degree_CTL == 0)) %>%
  mutate(
    HigherNetwork = ifelse(degree_AD > degree_CTL, "AD", "Control")
  )

saveRDS(deg_common_df,
        file.path(OUT_DIR, "Degree_AD_CTL_commonGenes_table.rds"))
write_tsv(deg_common_df,
          file.path(OUT_DIR, "Degree_AD_CTL_commonGenes_table.tsv"))

## =========================
## 6) Top 10 genes by |Z-score|
## =========================
top10_z <- deg_common_nofloat %>%
  arrange(desc(abs(z_change))) %>%
  slice_head(n = 10)

## =========================
## 7) Plot: all genes + top10 by |Z| labeled
## =========================
p_deg_all_top10 <- ggplot(deg_common_nofloat,
                          aes(x = degree_AD,
                              y = degree_CTL)) +
  # todos los genes (fondo)
  geom_point(aes(color = HigherNetwork),
             alpha = 0.25, size = 1.5) +
  # top 10 resaltados
  geom_point(data = top10_z,
             aes(color = HigherNetwork),
             alpha = 0.95, size = 4) +
  # diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  # labels solo para top 10, en negro
  geom_text_repel(
    data = top10_z,
    aes(label = gene),
    color = "black",
    size = 3.5,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.alpha = 0.7
  ) +
  scale_color_manual(values = c("AD" = "#E41A1C", "Control" = "#377EB8")) +
  labs(
    title = "Degree comparison (AD vs Control)\nTop 10 genes by |Z-score of Δdegree|",
    x = "Degree in AD network",
    y = "Degree in Control network",
    color = "Higher-degree\nnetwork"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "Degree_AD_vs_Control_all_top10_Zscore.pdf"),
       p_deg_all_top10, width = 7, height = 6)
ggsave(file.path(OUT_DIR, "Degree_AD_vs_Control_all_top10_Zscore.png"),
       p_deg_all_top10, width = 7, height = 6, dpi = 300)

###############################################
## END
###############################################
