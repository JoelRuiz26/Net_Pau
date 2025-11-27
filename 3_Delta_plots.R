# ===================== LIBRERÍAS ===================== #
library(tidyverse)
library(ggrepel)

# ===================== CARGA DE DATOS ===================== #
load("/STORAGE/csbig/jruiz/Redes_Pau/Delta_image.RData")  

# ===================== CONTRIBUCIÓN POR GEN ===================== #
Delta_diff <- sqrt(S_AD) - sqrt(S_CTL)
gene_contrib <- rowSums(Delta_diff^2)
total_contrib <- sum(gene_contrib)

# Tabla con contribuciones relativas y Z robusto
contrib_df <- tibble(
  Gene = rownames(S_AD),
  Contribution = gene_contrib,
  Percent = 100 * gene_contrib / total_contrib
) %>%
  mutate(Z_robust = (Contribution - median(Contribution)) / mad(Contribution)) %>%
  arrange(desc(Contribution))

# Selecciona genes que más cambian (Z_robust > 3)
top_genes <- contrib_df %>% filter(Z_robust > 4) %>% pull(Gene)

# ===================== PLOT DE CAMBIO DE AFINIDAD ===================== #
Affinity_df <- tibble(
  Gene = top_genes,
  AD = rowSums(sqrt(S_AD)[top_genes, ]),      # afinidad en AD
  CTL = rowSums(sqrt(S_CTL)[top_genes, ])     # afinidad en Control
) %>%
  mutate(Group = ifelse(AD > CTL, "AD", "CTL"))

dim(Affinity_df)
#[1] 377   4
# ===================== PLOT ===================== #

p_affinity <- ggplot(Affinity_df, aes(x = AD, y = CTL)) + 
  geom_point(aes(color = Group), size = 5, alpha = 0.9) +
  geom_text_repel(
    aes(label = Gene),
    color = "black",
    fontface = "bold",
    max.overlaps = 80,
    size = 3.5,
    segment.color = "gray55",
    segment.alpha = 0.6,
    box.padding = 0.3
  ) +
  scale_color_manual(values = c("AD" = "#E41A1C", "CTL" = "#377EB8")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Top Genes With Strongest Affinity Changes (Robust Z > 4)",
    x = "Affinity in AD Network",
    y = "Affinity in Control Network",
    color = "Higher Affinity\nNetwork"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# ===================== SAVE PLOT ===================== #
OUT_DIR <- "~/Net_Pau/Affinity_plots/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(OUT_DIR, "Affinity_AD_vs_Control.pdf"),
       p_affinity, width = 10, height = 8)

ggsave(file.path(OUT_DIR, "Affinity_AD_vs_Control.png"),
       p_affinity, width = 10, height = 8, dpi = 300)

