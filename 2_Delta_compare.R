# SCRIPT: DeltaCon Similarity Analysis between AD and CTL Graphs (Robust + Parallel)
# Author: Jose Joel Ruiz Hernandez | Date: 2025-06-05
#REFERENCE:
  # Koutra D, Vogelstein JT, Faloutsos C. DeltaCon: A Principled Massive-Graph Similarity Function.
  # In: Proceedings of the 2013 SIAM International Conference on Data Mining. 2013.

start_time <- Sys.time()

# -------------------------- LIBRARIES ----------------------------------------
library(igraph)
library(Matrix)
library(furrr)
library(future)

# ---------------------- PARALLEL PLAN (safe) ---------------------------------
plan(multisession, workers = 2)

# --------------------------- LOAD GRAPHS -------------------------------------
g_AD  <- readRDS("/STORAGE/csbig/jruiz/Redes_Pau/g_AD_200k.rds")
g_CTL <- readRDS("/STORAGE/csbig/jruiz/Redes_Pau/g_ctl_200k.rds")

AD_names <- V(g_AD)$name
ctl_names<- V(g_CTL)$name
common_genes <- sort(intersect(AD_names, ctl_names))

#cut graph
g_AD  <- igraph::induced_subgraph(g_AD,  vids = common_genes)
g_CTL <- igraph::induced_subgraph(g_CTL, vids = common_genes)

# Asegurar que ambos grafos tengan mismos nodos y mismo orden
nodes_order <- sort(union(igraph::V(g_AD)$name, igraph::V(g_CTL)$name))

g_AD  <- igraph::permute(g_AD,  base::match(nodes_order, igraph::V(g_AD)$name))
g_CTL <- igraph::permute(g_CTL, base::match(nodes_order, igraph::V(g_CTL)$name))

# ----------------- FUNCTION: Compute Affinity Matrix (DeltaCon) --------------
compute_similarity_matrix_parallel <- function(graph) {
  message("🧠 Calculando matriz de afinidad con ΔCon...")
  
  # epsilon dinámico según conectividad del grafo
  epsilon <- 1 / (1 + max(igraph::degree(graph)))
  
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
  D <- Diagonal(x = Matrix::rowSums(A))
  I <- Diagonal(n = nrow(A))
  M <- I + (epsilon^2) * D - epsilon * A
  
  n <- ncol(M)
  workers <- future::nbrOfWorkers()
  block_size <- ceiling(n / workers)
  blocks <- split(seq_len(n), ceiling(seq_len(n) / block_size))
  
  solve_block <- function(cols) {
    tryCatch({
      solve(M, I[, cols, drop = FALSE])
    }, error = function(e) {
      message("❌ Error resolviendo columnas: ", paste(head(cols, 5), "..."))
      message("→ ", e$message)
      return(matrix(NA, nrow = nrow(M), ncol = length(cols)))
    })
  }
  
  result_blocks <- future_map(blocks, solve_block,
                              .options = furrr_options(seed = TRUE))
  S <- do.call(cbind, result_blocks)
  rownames(S) <- colnames(S) <- igraph::V(graph)$name
  return(S)
}

# -------------------- COMPUTE AFFINITY MATRICES ------------------------------
S_AD  <- compute_similarity_matrix_parallel(g_AD)
S_CTL <- compute_similarity_matrix_parallel(g_CTL)

# -------------------- ALIGN MATRICES TO SAME ORDER ---------------------------
genes <- sort(intersect(rownames(S_AD), rownames(S_CTL)))
S_AD  <- S_AD[genes, genes]
S_CTL <- S_CTL[genes, genes]

# -------------------- DELTACON DISTANCE --------------------------------------
delta_con_distance <- function(S1, S2) {
  sqrt(sum((sqrt(S1) - sqrt(S2))^2, na.rm = TRUE))
}

delta_con <- delta_con_distance(S_AD, S_CTL)
delta_con_similarity <- 1 / (1 + delta_con)

# ------------------------- OUTPUT --------------------------------------------
cat("✅ DeltaCon Distance (AD vs CTL): ", round(delta_con, 4), "\n")
cat("✅ DeltaCon Similarity (AD vs CTL): ", round(delta_con_similarity, 6), "\n")

# ------------------------- SAVE OUTPUT ---------------------------------------
out_dir <- "~/Pau/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(S_AD,              file.path(out_dir, "S_AD_matrix.rds"))
saveRDS(S_CTL,             file.path(out_dir, "S_CTL_matrix.rds"))
saveRDS(delta_con,         file.path(out_dir, "DeltaCon_distance_AD_CTL.rds"))
saveRDS(delta_con_similarity,
        file.path(out_dir, "DeltaCon_similarity_AD_CTL.rds"))
save.image("~/Pau/Delta_image.RData")

end_time <- Sys.time()
cat("🕒 Total execution time:", round(end_time - start_time, 2), "seconds\n")
