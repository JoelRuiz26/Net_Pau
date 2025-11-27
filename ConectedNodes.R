# --------------------------------------------------------
# Determine the minimal number of edges required to 
# cover all nodes in Co-expression Networks (ARACNE sorted)
#
# Author: Joel Ruiz Hernandez
#
# Description:
#   - For each network: identifies the row (edge count) 
#     at which all unique genes are included at least once.
# Input: 
#   -Edge list (.sort) with highest mutual information 
# Output:
#   - Reports summary statistics to the console.
#   - Saves the maximum cutoff value between both networks
#     as a single .rds object.
# --------------------------------------------------------

# --------------------------------------------------------
# Load required libraries
# --------------------------------------------------------
library(dplyr)
library(vroom)
library(igraph)
library(tools)

# --------------------------------------------------------
# Function to analyze a single network file
# Assumes the input file is pre-sorted by MI
# --------------------------------------------------------
#                               Set a huge slice of edgelist that could have all nodes
analyze_network <- function(file_path, max_rows = 100000000) {
  
  # Load the edge list (up to max_rows)
  edge_list <- vroom(
    file_path, 
    col_names = c("GenA", "GenB", "MI"),
    show_col_types = FALSE
  ) %>%
    na.omit() %>%
    slice(1:max_rows)
  
  # Extract all unique nodes (genes)
  all_nodes <- unique(c(edge_list$GenA, edge_list$GenB))
  total_nodes <- length(all_nodes)
  
  # Determine the cutoff row where all nodes are covered
  covered_nodes <- character()
  cutoff_row <- NA
  
  for (i in 1:nrow(edge_list)) {
    covered_nodes <- unique(c(covered_nodes, edge_list$GenA[i], edge_list$GenB[i]))
    
    if (length(covered_nodes) == total_nodes) {
      cutoff_row <- i
      break
    }
  }
  
  # Return structured summary
  return(list(
    network_name = file_path_sans_ext(basename(file_path)),
    total_nodes = total_nodes,
    cutoff_row = cutoff_row,
    covered_nodes = length(covered_nodes),
    coverage_percent = round(length(covered_nodes) / total_nodes * 100, 2),
    sample_size = nrow(edge_list)
  ))
}


# --------------------------------------------------------
# Analyze both networks
# --------------------------------------------------------
results_control <- analyze_network("/STORAGE/csbig/ROSMAP_HCN_counts_control.sort")

results_AD <- analyze_network("/STORAGE/csbig/ROSMAP_HCN_counts_AD.sort")

# --------------------------------------------------------
# Function to generate and print a detailed report
# --------------------------------------------------------
generate_report <- function(results) {
  report <- paste0(
    "\n════════════════════════════════════\n",
    " NETWORK ANALYSIS REPORT - ", results$network_name, "\n",
    "════════════════════════════════════\n",
    "• Total unique nodes: ", results$total_nodes, "\n",
    "• Rows analyzed: ", results$sample_size, "\n",
    "• Cutoff row found: ", 
    ifelse(is.na(results$cutoff_row), "NOT COMPLETE", results$cutoff_row), "\n",
    "• Nodes covered: ", results$covered_nodes, "/", results$total_nodes, "\n",
    "• Coverage percent: ", results$coverage_percent, "%\n",
    "════════════════════════════════════\n"
  )
  
  cat(report)
  invisible(report)
}

# --------------------------------------------------------
# Print summary reports to console
# --------------------------------------------------------
generate_report(results_control)
generate_report(results_AD)

# --------------------------------------------------------
# Save INDIVIDUAL cutoff values (A7 y A9)
# --------------------------------------------------------
cutoff_control <- results_control$cutoff_row
cutoff_AD <- results_AD$cutoff_row

saveRDS(
  cutoff_control, 
  file = "/STORAGE/csbig/control_cutoff_links.rds"
)

saveRDS(
  cutoff_AD, 
  file = "/STORAGE/csbig/AD_cutoff_links.rds"
)

cat("\n════════════════════════════════════\n")
cat(" INDIVIDUAL THRESHOLDS SAVED\n")
cat("════════════════════════════════════\n")
cat("A7 cutoff row: ", cutoff_control, "\n")
cat("A9 cutoff row: ", cutoff_AD, "\n")
cat("════════════════════════════════════\n")

save.image("/STORAGE/csbig/cut.RData")

