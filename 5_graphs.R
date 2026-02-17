#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(stringr)
})

# ===================== CONFIG ===================== #
edge_dir    <- "/STORAGE/csbig/jruiz/Redes_Pau/1_Edgelist"
df_top_path <- "/STORAGE/csbig/jruiz/Redes_Pau/4_DGE_centrality/4_2_plot_fout_scatter/4_2_1_top10_genes.rds"
common_dir  <- "/STORAGE/csbig/jruiz/Redes_Pau/4_DGE_centrality"

out_dir <- "/STORAGE/csbig/jruiz/Redes_Pau/5_graph_interest_genes"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

REG_ORDER <- c("PCC","HCN","DLPFC","TC","CRB")

region2common <- c(
  PCC   = "PCC_degree_COMMON.tsv",
  HCN   = "HCN_degree_COMMON.tsv",
  DLPFC = "DLPFC_degree_COMMON.tsv",
  TC    = "MAYO_TC_degree_COMMON.tsv",
  CRB   = "MAYO_CRB_degree_COMMON.tsv"
)

# ===================== HELPERS ===================== #
read_edgelist <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  setnames(dt, old = names(dt), new = c("gene1","gene2","MI")[seq_along(names(dt))])
  dt
}

parse_edge_filename <- function(f) {
  bn <- basename(f)
  reg <- str_match(bn, "_(PCC|HCN|DLPFC|TC|CRB)_")[,2]
  cond <- case_when(
    grepl("_AD_", bn) ~ "AD",
    grepl("_control_", bn) ~ "control",
    TRUE ~ NA_character_
  )
  list(region = reg, condition = cond)
}

filter_edges_1hop <- function(dt_edges, genes_interest) {
  genes_interest <- unique(na.omit(genes_interest))
  if (length(genes_interest) == 0) return(dt_edges[0])
  dt_edges[gene1 %chin% genes_interest | gene2 %chin% genes_interest]
}

# Carga COMMON por región: AHORA incluye is_dge
load_common_map <- function(region) {
  f <- file.path(common_dir, region2common[[region]])
  if (!file.exists(f)) stop(sprintf("No existe COMMON para %s: %s", region, f))
  
  dt <- fread(f, sep = "\t", header = TRUE, showProgress = FALSE)
  
  req <- c("gene_clean","label","logFC","is_dge")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop(sprintf("COMMON %s no tiene columnas: %s", f, paste(miss, collapse=", ")))
  
  dt <- dt[, .(
    ensg  = gene_clean,
    label = as.character(label),
    logFC = as.numeric(logFC),
    is_dge = as.logical(is_dge)
  )]
  
  unique(dt)
}

# ===================== LOAD INPUTS ===================== #
df_top <- readRDS(df_top_path)

genes_by_region <- df_top %>%
  mutate(region = as.character(region)) %>%
  group_by(region) %>%
  summarise(genes_interest = list(unique(na.omit(gene_clean))), .groups = "drop") %>%
  filter(region %in% REG_ORDER)

# ===================== 1) BUILD FILTERED EDGE LISTS ===================== #
edge_files <- list.files(edge_dir, pattern = "\\.tsv$", full.names = TRUE)

filtered_edges <- setNames(vector("list", length(REG_ORDER)), REG_ORDER)
summary_rows <- list()

for (f in edge_files) {
  info <- parse_edge_filename(f)
  reg  <- info$region
  cond <- info$condition
  if (is.na(reg) || is.na(cond) || !(reg %in% REG_ORDER)) next
  
  genes_interest <- genes_by_region$genes_interest[[which(genes_by_region$region == reg)]]
  
  message(sprintf(">> FILTER %s | %s | %s", reg, cond, basename(f)))
  dt_edges <- read_edgelist(f)
  dt_sub   <- filter_edges_1hop(dt_edges, genes_interest)
  
  if (is.null(filtered_edges[[reg]])) filtered_edges[[reg]] <- list()
  filtered_edges[[reg]][[cond]] <- dt_sub
  
  fwrite(dt_sub, file.path(out_dir, sprintf("%s_%s_edgelist_filtered.tsv", reg, cond)), sep = "\t")
  
  nodes_sub <- sort(unique(c(dt_sub$gene1, dt_sub$gene2)))
  fwrite(data.table(node = nodes_sub), file.path(out_dir, sprintf("%s_%s_nodes_filtered.tsv", reg, cond)), sep = "\t")
  
  summary_rows[[length(summary_rows)+1]] <- data.table(
    region = reg, condition = cond,
    input_edges = nrow(dt_edges),
    filtered_edges = nrow(dt_sub),
    filtered_nodes = length(nodes_sub),
    n_interest = length(genes_interest),
    n_interest_present = sum(genes_interest %chin% nodes_sub),
    n_interest_missing = sum(!(genes_interest %chin% nodes_sub))
  )
}

summary_df <- rbindlist(summary_rows, fill = TRUE)
summary_df[, region := factor(region, levels = REG_ORDER)]
setorder(summary_df, region, condition)
print(as.data.frame(summary_df))

# ===================== 2) BUILD ONE GRAPH PER REGION ===================== #
graphs_by_region <- list()
graphs_summary   <- list()

for (reg in REG_ORDER) {
  
  e_ad  <- filtered_edges[[reg]][["AD"]]
  e_ctl <- filtered_edges[[reg]][["control"]]
  
  if (is.null(e_ad) && is.null(e_ctl)) {
    message(sprintf("[GraphML] %s: no filtered edges, skipping.", reg))
    next
  }
  
  if (!is.null(e_ad))  e_ad[, cond := "AD"]
  if (!is.null(e_ctl)) e_ctl[, cond := "control"]
  
  e_all <- rbindlist(Filter(Negate(is.null), list(e_ad, e_ctl)), use.names = TRUE, fill = TRUE)
  
  e_all[, a := pmin(gene1, gene2)]
  e_all[, b := pmax(gene1, gene2)]
  
  e_wide <- dcast(
    e_all,
    a + b ~ cond,
    value.var = "MI",
    fun.aggregate = function(x) if (length(x) == 0) NA_real_ else max(x, na.rm = TRUE)
  )
  
  e_wide[, present := fifelse(!is.na(AD) & !is.na(control), "both",
                              fifelse(!is.na(AD), "AD",
                                      fifelse(!is.na(control), "control", NA_character_)))]
  e_wide[, MI := fifelse(!is.na(AD), AD, control)]
  
  g <- graph_from_data_frame(
    d = as.data.frame(e_wide[, .(from = a, to = b, MI = MI, MI_AD = AD, MI_control = control, present = present)]),
    directed = FALSE
  )
  
  # --- base ---
  V(g)$ensg <- as.character(V(g)$name)
  
  nodes_ad  <- if (!is.null(e_ad))  unique(c(e_ad$gene1,  e_ad$gene2))  else character(0)
  nodes_ctl <- if (!is.null(e_ctl)) unique(c(e_ctl$gene1, e_ctl$gene2)) else character(0)
  
  V(g)$member <- ifelse(V(g)$ensg %in% nodes_ad & V(g)$ensg %in% nodes_ctl, "both",
                        ifelse(V(g)$ensg %in% nodes_ad, "AD",
                               ifelse(V(g)$ensg %in% nodes_ctl, "control", NA_character_)))
  
  # --- COMMON: label + logFC + is_dge ---
  common_map <- load_common_map(reg)
  
  lab_map <- setNames(common_map$label, common_map$ensg)
  lfc_map <- setNames(common_map$logFC, common_map$ensg)
  dge_map <- setNames(common_map$is_dge, common_map$ensg)
  
  V(g)$display_label <- ifelse(V(g)$ensg %in% names(lab_map), lab_map[V(g)$ensg], V(g)$ensg)
  
  # >>> CAMBIO: si is_dge FALSE -> logFC = 0; si TRUE -> logFC real; si NA -> NA
  V(g)$is_dge <- ifelse(V(g)$ensg %in% names(dge_map), as.logical(dge_map[V(g)$ensg]), NA)
  raw_lfc <- ifelse(V(g)$ensg %in% names(lfc_map), as.numeric(lfc_map[V(g)$ensg]), NA_real_)
  V(g)$logFC <- ifelse(is.na(V(g)$is_dge), raw_lfc,
                       ifelse(V(g)$is_dge, raw_lfc, 0))
  
  # interés (desde df_top)
  genes_interest <- genes_by_region$genes_interest[[which(genes_by_region$region == reg)]]
  V(g)$is_interest <- V(g)$ensg %chin% genes_interest
  
  # renombra a label (único)
  V(g)$name <- make.unique(as.character(V(g)$display_label), sep = "_")
  V(g)$label <- V(g)$name
  
  graphs_by_region[[reg]] <- g
  
  graphs_summary[[reg]] <- data.table(
    region = reg,
    nodes = vcount(g),
    edges = ecount(g),
    n_interest_in_graph = sum(V(g)$is_interest, na.rm = TRUE),
    member_AD = sum(V(g)$member == "AD", na.rm = TRUE),
    member_control = sum(V(g)$member == "control", na.rm = TRUE),
    member_both = sum(V(g)$member == "both", na.rm = TRUE),
    edges_AD = sum(E(g)$present == "AD", na.rm = TRUE),
    edges_control = sum(E(g)$present == "control", na.rm = TRUE),
    edges_both = sum(E(g)$present == "both", na.rm = TRUE)
  )
  
  out_gml <- file.path(out_dir, sprintf("%s_AD_Control.graphml", reg))
  write_graph(g, out_gml, format = "graphml")
  message(sprintf("[GraphML] saved: %s", out_gml))
}

print(head(graphs_by_region))

graphs_summary_df <- rbindlist(graphs_summary, fill = TRUE)
graphs_summary_df[, region := factor(region, levels = REG_ORDER)]
setorder(graphs_summary_df, region)
print(as.data.frame(graphs_summary_df))

saveRDS(filtered_edges, file.path(out_dir, "5_1_Edgelist_byRegion_topGenes.rds"))
saveRDS(graphs_by_region, file.path(out_dir, "5_2_graphs_by_region_ADplusControl.rds"))
