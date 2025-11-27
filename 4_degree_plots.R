###############################################
## DEGREE COMPARISON BETWEEN AD AND CTL NETWORKS
## - Compute degree in original graphs
## - Restrict to common genes
## - Δdegree + Z-score
## - Elegant plots
###############################################

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

## =========================
## 0) Paths and output dir
## =========================
OUT_DIR <- "~/Pau/Degree_comparison/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## =========================
## 1) Load graphs (original, sin recortar)
## =========================
g_AD  <- readRDS("/STORAGE/csbig/jruiz/Redes_Pau/g_AD_200k.rds")
g_CTL <- readRDS("/STORAGE/csbig/jruiz/Redes_Pau/g_ctl_200k.rds")

## =========================
## 2) Degree in original graphs
## =========================
deg_AD_full  <- igraph::degree(g_AD)
deg_CTL_full <- igraph::degree(g_CTL)

AD_names  <- names(deg_AD_full)
CTL_names <- names(deg_CTL_full)

## Genes comunes (existentes en ambas redes)
common_genes <- intersect(AD_names, CTL_names)
length(common_genes)  # Quick check

## =========================
## 3) Build degree table (usando grados originales)
## =========================
deg_df <- tibble(
  gene       = common_genes,
  degree_AD  = as.numeric(deg_AD_full[common_genes]),
  degree_CTL = as.numeric(deg_CTL_full[common_genes])
) %>%
  mutate(
    delta_degree     = degree_AD - degree_CTL,
    abs_delta_degree = abs(delta_degree),
    mean_degree      = (degree_AD + degree_CTL) / 2
  )

## Opcional: quitar genes con grado 0 en ambas redes
deg_df_nofloat <- deg_df %>%
  filter(!(degree_AD == 0 & degree_CTL == 0))

## =========================
## 4) Z-score del cambio de grado
## =========================
deg_df <- deg_df %>%
  mutate(
    z_change = (delta_degree - mean(delta_degree)) / sd(delta_degree)
  )

deg_df_nofloat <- deg_df %>%
  filter(!(degree_AD == 0 & degree_CTL == 0))

## Guardar tabla
saveRDS(deg_df, file.path(OUT_DIR, "Degree_AD_CTL_table.rds"))
readr::write_tsv(deg_df, file.path(OUT_DIR, "Degree_AD_CTL_table.tsv"))

## =========================
## 5) Top genes que más cambian
## =========================
top_changed_abs <- deg_df %>%
  arrange(desc(abs_delta_degree)) %>%
  slice_head(n = 50)

top_changed_z <- deg_df %>%
  arrange(desc(abs(z_change))) %>%
  slice_head(n = 50)

readr::write_tsv(top_changed_abs,
                 file.path(OUT_DIR, "Top50_genes_abs_delta_degree.tsv"))
readr::write_tsv(top_changed_z,
                 file.path(OUT_DIR, "Top50_genes_z_change.tsv"))

## =========================
## 6) Plots elegantes
## =========================

## 6.1) Distribución de grados (log10) AD vs CTL
deg_long <- deg_df_nofloat %>%
  select(gene, degree_AD, degree_CTL) %>%
  pivot_longer(cols = starts_with("degree_"),
               names_to = "network",
               values_to = "degree") %>%
  mutate(
    network = recode(network,
                     "degree_AD"  = "AD",
                     "degree_CTL" = "CTL")
  )

p_deg_density <- ggplot(deg_long, aes(x = degree + 1, color = network, fill = network)) +
  geom_density(alpha = 0.15, adjust = 1.2) +
  scale_x_log10() +
  labs(
    title = "Degree distribution per network",
    x = "Degree (log10(degree + 1))",
    y = "Density",
    color = "Network",
    fill  = "Network"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(file.path(OUT_DIR, "Degree_distribution_AD_vs_CTL.pdf"),
       p_deg_density, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "Degree_distribution_AD_vs_CTL.png"),
       p_deg_density, width = 7, height = 5, dpi = 300)

## 6.2) Scatter degree_CTL vs degree_AD (log–log), coloreado por Z-score
p_scatter <- ggplot(deg_df_nofloat,
                    aes(x = degree_CTL + 1, y = degree_AD + 1, color = z_change)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradient2(
    low  = "#2166AC",
    mid  = "grey90",
    high = "#B2182B",
    midpoint = 0
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Node degree comparison (AD vs CTL)",
    x = "Degree in CTL (log10 + 1)",
    y = "Degree in AD (log10 + 1)",
    color = "Z-change\n(Δdegree)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "Degree_scatter_AD_vs_CTL.pdf"),
       p_scatter, width = 6.5, height = 6)
ggsave(file.path(OUT_DIR, "Degree_scatter_AD_vs_CTL.png"),
       p_scatter, width = 6.5, height = 6, dpi = 300)

## 6.3) Volcano-style: Δdegree vs mean_degree, coloreado por z_change
p_volcano <- ggplot(deg_df_nofloat,
                    aes(x = mean_degree + 1, y = delta_degree, color = z_change)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 1) +
  scale_x_log10() +
  scale_color_gradient2(
    low  = "#2166AC",
    mid  = "grey90",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    title = "Δdegree vs mean degree (AD – CTL)",
    x = "Mean degree (AD & CTL, log10 + 1)",
    y = "Δdegree (AD – CTL)",
    color = "Z-change"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "Degree_volcano_AD_minus_CTL.pdf"),
       p_volcano, width = 7, height = 5.5)
ggsave(file.path(OUT_DIR, "Degree_volcano_AD_minus_CTL.png"),
       p_volcano, width = 7, height = 5.5, dpi = 300)

## 6.4) Scatter con labels para top genes (por |z_change|)
top_label <- deg_df_nofloat %>%
  arrange(desc(abs(z_change))) %>%
  slice_head(n = 30)

p_scatter_labels <- p_scatter +
  geom_text_repel(
    data = top_label,
    aes(label = gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.alpha = 0.6
  ) +
  labs(title = "Node degree comparison (AD vs CTL) – top |Z-change| genes")

ggsave(file.path(OUT_DIR, "Degree_scatter_AD_vs_CTL_topZ_labels.pdf"),
       p_scatter_labels, width = 7.5, height = 7)
ggsave(file.path(OUT_DIR, "Degree_scatter_AD_vs_CTL_topZ_labels.png"),
       p_scatter_labels, width = 7.5, height = 7, dpi = 300)

###############################################
## END
###############################################
