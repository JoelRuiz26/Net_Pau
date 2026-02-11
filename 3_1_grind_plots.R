#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(grid)
  library(cowplot)
  library(stringr)
})

# ===================== CONFIG ===================== #
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"

volcano_rds <- file.path(base_dir, "2_1_DGE_limma/2_2_DGE_plots/VolcanoRow_byRegion_AD_vs_Control_plot.rds")
heat_rds    <- file.path(base_dir, "3_GSEA_REACTOME/3_3_Heatmap_Stars_alpha_plot.rds")

out_dir <- file.path(base_dir, "3_1_DGE_GSEA_plot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(out_dir, "DGE_GSEA.png")
out_pdf <- file.path(out_dir, "DGE_GSEA.pdf")
out_rds <- file.path(out_dir, "DGE_GSEA.rds")

REG_ORDER  <- c("CRB","TC","DLPFC","HCN","PCC")
STRIP_FILL <- "grey92"

# ===================== SIZES (TUNED) ===================== #
# Volcano
VOL_BASE       <- 14.5
VOL_STRIP_SZ   <- 15.5   # títulos de panel del volcano (PCC/HCN/...)
VOL_AXIS_TXT   <- 12.5
VOL_AXIS_TTL   <- 13.5

# Heatmap
HEAT_BASE      <- 14.5
HEAT_STRIP_SZ  <- 15.0   # títulos de panel del heatmap (CRB/TC/...) en “strip”
HEAT_Y_SZ      <- 10.2   # procesos
STAR_SZ        <- 4.4

# Legend typography
LEG_TITLE_SZ <- 13
LEG_TEXT_SZ  <- 12

# Wrap width (más ancho = menos saltos de línea)
WRAP_WIDTH <- 62

# Output device size
OUT_W_IN <- 24
OUT_H_IN <- 14.5
OUT_DPI  <- 600

# ===================== LOAD ===================== #
p_volcano0 <- readRDS(volcano_rds)
p_heat_old <- readRDS(heat_rds)

# ===================== HEATMAP DATA ===================== #
# Importante: mantener region en factor con el mismo orden
df2 <- p_heat_old$data %>%
  mutate(
    region       = factor(region, levels = REG_ORDER),
    term_wrapped = stringr::str_wrap(term_simple, width = WRAP_WIDTH),
    x1           = 1
  )

# ===================== HEATMAP (WITH GREY STRIP + NO “EMPTY NECK”) ===================== #
# Clave para:
# 1) tener strip gris como el volcano => facet_grid(. ~ region) + strip.background
# 2) NO tener “cuello” (espacio vacío enorme) => y como factor con drop=TRUE y recorte de panel
# 3) labels de términos a la derecha => scale_y_discrete(position="right")
p_heat <- ggplot(df2, aes(x = x1, y = term_wrapped)) +
  geom_tile(
    aes(fill = NES, alpha = alpha_sig),
    color = "grey92", linewidth = 0.4
  ) +
  geom_text(
    aes(label = star),
    fontface = "bold", size = STAR_SZ, color = "grey10"
  ) +
  facet_grid(. ~ region, drop = TRUE, scales = "free_x", space = "free_x") +
  scale_y_discrete(
    position = "right",
    expand = expansion(mult = c(0.01, 0.01))   # <- menos aire vertical
  ) +
  # 🔥 ESTO QUITA EL “CUELLO” HORIZONTAL
  scale_x_continuous(
    NULL,
    breaks = NULL,
    limits = c(0.5, 1.5),                      # <- fuerza el panel a la columna del tile
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_gradient2(
    name = "NES",
    low  = "#2C7FB8",
    mid  = "#F7F7F7",
    high = "#D7301F",
    midpoint = 0
  ) +
  scale_alpha(guide = "none") +
  theme_bw(base_size = HEAT_BASE) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    
    axis.text.y  = element_text(
      size = HEAT_Y_SZ,
      face = "bold",
      lineheight = 1.10,
      margin = margin(l = 8)
    ),
    axis.ticks.y = element_line(linewidth = 0.3),
    
    strip.background = element_rect(fill = STRIP_FILL, color = NA),
    strip.text.x     = element_text(face = "bold", size = HEAT_STRIP_SZ),
    
    panel.spacing.x  = unit(0.15, "lines"),
    plot.margin      = margin(t = 8, r = 10, b = 8, l = 8),
    
    legend.position  = "none"
  )

# ===================== VOLCANO (STYLE + NO LEGEND IN PANEL) ===================== #
p_volcano <- p_volcano0 +
  theme_bw(base_size = VOL_BASE) +
  theme(
    strip.background = element_rect(fill = STRIP_FILL, color = NA),
    strip.text       = element_text(face = "bold", size = VOL_STRIP_SZ),
    
    axis.text  = element_text(size = VOL_AXIS_TXT),
    axis.title = element_text(size = VOL_AXIS_TTL, face = "bold"),
    
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8),
    legend.position = "none"
  )

# ===================== LEGENDS COLUMN (CENTERED) ===================== #
leg_volcano <- cowplot::get_legend(
  p_volcano0 +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = LEG_TITLE_SZ),
      legend.text  = element_text(face = "bold", size = LEG_TEXT_SZ),
      legend.key.size = unit(0.9, "lines")
    )
)

p_heat_leg <- ggplot(df2, aes(x = x1, y = term_simple)) +
  geom_tile(aes(fill = NES)) +
  scale_fill_gradient2(
    name = "NES",
    low  = "#2C7FB8",
    mid  = "#F7F7F7",
    high = "#D7301F",
    midpoint = 0
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = LEG_TITLE_SZ),
    legend.text  = element_text(face = "bold", size = LEG_TEXT_SZ),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width  = unit(0.9, "lines")
  )

leg_heat <- cowplot::get_legend(p_heat_leg)

p_legends <- (
  plot_spacer() /
    wrap_elements(leg_volcano) /
    plot_spacer() /
    plot_spacer() /
    wrap_elements(leg_heat) /
    plot_spacer()
) +
  plot_layout(heights = c(0.35, 1, 0.35, 0.35, 1, 0.35)) &
  theme(plot.margin = margin(0, 0, 0, 0))

# ===================== COMBINE ===================== #
p_combined <- (p_volcano | p_legends | p_heat) +
  plot_layout(widths = c(2.55, 0.55, 2.10)) +
  plot_annotation(NULL) &
  theme(plot.margin = margin(t = 6, r = 8, b = 6, l = 8))

# ===================== SAVE ===================== #
saveRDS(p_combined, out_rds)

png(out_png, width = OUT_W_IN, height = OUT_H_IN, units = "in",
    res = OUT_DPI, type = "cairo-png")
print(p_combined)
dev.off()

pdf(out_pdf, width = OUT_W_IN, height = OUT_H_IN, onefile = TRUE)
print(p_combined)
dev.off()

message("Saved:\n- ", out_png, "\n- ", out_pdf, "\n- ", out_rds)
