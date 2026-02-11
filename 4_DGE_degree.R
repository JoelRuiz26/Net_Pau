suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(biomaRt)
})

options(stringsAsFactors = FALSE)

# ============================================================
# CONFIG
# ============================================================
base_dir <- "/STORAGE/csbig/jruiz/Redes_Pau"

FDR_CUT <- 0.05
LFC_CUT <- 0.5

out_dir <- file.path(base_dir, "4_DGE_centrality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map_cache <- file.path(out_dir, "ensembl_to_hgnc_cache.rds")

# ============================================================
# HELPERS
# ============================================================
norm_key <- function(x) {
  x %>%
    toupper() %>%
    str_replace("^ROSMAP_", "") %>%
    str_replace("^MAYO_", "") %>%
    str_replace("_DEGREE_COMMON$", "") %>%
    str_replace("_LIMMA_AD_VS_CONTROL_ALLGENES$", "") %>%
    str_replace("_AD_VS_CONTROL_ALLGENES$", "") %>%
    str_replace("_ALLGENES$", "") %>%
    str_trim()
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

make_deg_with_n_dge <- function(deg_df, dge_df, fdr_cut = 0.05, lfc_cut = 0.5) {
  # DGE completo (para traer logFC por gen; no solo significativos)
  dge_small <- dge_df %>%
    transmute(
      gene_clean = sub("\\..*$", "", gene),
      logFC = as.numeric(logFC),
      adj.P.Val = as.numeric(adj.P.Val)
    ) %>%
    distinct(gene_clean, .keep_all = TRUE)
  
  # set de significativos para is_dge
  sig_set <- dge_small %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    filter(adj.P.Val <= fdr_cut, abs(logFC) >= lfc_cut) %>%
    pull(gene_clean)
  
  deg_df %>%
    mutate(
      gene_clean = sub("\\..*$", "", gene),
      delta_abs_percentil = rank(abs(delta_perc), ties.method = "average") / n()
    ) %>%
    # agrega logFC (LFC) por gen
    left_join(dge_small %>% dplyr::select(gene_clean, logFC), by = "gene_clean") %>%
    mutate(
      is_dge = gene_clean %in% sig_set
    ) %>%
    arrange(desc(delta_abs_percentil), desc(abs(delta_perc))) %>%
    mutate(n_dge = cumsum(is_dge))
}

# ============================================================
# 1) LOAD DGE LIST
# ============================================================
dge_files <- list.files(
  file.path(base_dir, "2_1_DGE_limma"),
  pattern = "_limma_AD_vs_Control_allGenes\\.tsv$",
  full.names = TRUE
)
dge_list <- set_names(
  map(dge_files, ~ read_tsv(.x, col_types = cols(), progress = FALSE)),
  tools::file_path_sans_ext(basename(dge_files))
)

# ============================================================
# 2) LOAD DEG LIST
# ============================================================
deg_files <- list.files(
  file.path(base_dir, "1_1_Centrality_comparison_rho", "degree"),
  pattern = "_degree_COMMON\\.tsv$",
  full.names = TRUE
)
deg_list <- set_names(
  map(deg_files, ~ read_tsv(.x, col_types = cols(), progress = FALSE)),
  tools::file_path_sans_ext(basename(deg_files))
)

# ============================================================
# 3) MATCH REGIONS
# ============================================================
deg_key_tbl <- tibble(deg_name = names(deg_list)) %>% mutate(join_key = norm_key(deg_name))
dge_key_tbl <- tibble(dge_name = names(dge_list)) %>% mutate(join_key = norm_key(dge_name))
pairs <- inner_join(deg_key_tbl, dge_key_tbl, by = "join_key")
if (nrow(pairs) == 0) stop("No pude emparejar DEG vs DGE por región. Revisa nombres en deg_list y dge_list.")

# ============================================================
# 4) BUILD deg_list2 (ahora incluye logFC)
# ============================================================
deg_list2 <- list()
for (i in seq_len(nrow(pairs))) {
  deg_name <- pairs$deg_name[i]
  dge_name <- pairs$dge_name[i]
  deg_list2[[deg_name]] <- make_deg_with_n_dge(
    deg_df = deg_list[[deg_name]],
    dge_df = dge_list[[dge_name]],
    fdr_cut = FDR_CUT,
    lfc_cut = LFC_CUT
  )
}

# ============================================================
# 5) HGNC ANNOTATION (MINIMAL)
# ============================================================
all_ens <- unique(unlist(map(deg_list2, ~ sub("\\..*$", "", .x$gene_clean))))
all_ens <- all_ens[!is.na(all_ens) & all_ens != ""]

if (file.exists(map_cache)) {
  message("Loading cached HGNC mapping: ", map_cache)
  map_tbl <- readRDS(map_cache)
} else {
  message("Querying biomaRt for HGNC symbols...")
  map_tbl <- map_ensembl_to_hgnc(all_ens)
  saveRDS(map_tbl, map_cache)
  message("Saved mapping cache: ", map_cache)
}

deg_list2 <- map(
  deg_list2,
  ~ .x %>%
    left_join(map_tbl, by = c("gene_clean" = "ensembl_gene_id")) %>%
    mutate(label = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, gene_clean))
)

# ============================================================
# 6) SAVE
# ============================================================
walk2(
  deg_list2,
  names(deg_list2),
  ~ write_tsv(.x, file.path(out_dir, paste0(.y, ".tsv")))
)

message("DONE ✅ Saved annotated tables (with logFC) in: ", out_dir)
