#Script for compare node centrality between networks
library(igraph)
library(dplyr)
library(vroom)
###load edge list 

AD <-  vroom("/STORAGE/csbig/jruiz/Redes_Pau/ROSMAP_PCC_counts_AD.sort",col_names = c("GenA", "GenB", "MI")) %>% dplyr::slice(1:200000)
ctl <-  vroom("/STORAGE/csbig/jruiz/Redes_Pau/ROSMAP_PCC_counts_control.sort" ,col_names = c("GenA", "GenB", "MI")) %>% dplyr::slice(1:200000)

#########
#Extract all unique and shared nodes
AD_nodes <- unique(c(AD$GenA, AD$GenB))
ctl_nodes <- unique(c(ctl$GenA, ctl$GenB))
shared_nodes <- intersect(AD_nodes, ctl_nodes)

AD_unique<- setdiff(AD_nodes, shared_nodes)
ctl_unique<- setdiff(ctl_nodes, shared_nodes)

#########
#Get edgelist only with common nodes
AD_shared  <- AD  %>% dplyr::filter(GenA %in% shared_nodes & GenB %in% shared_nodes)
ctl_shared <- ctl %>% dplyr::filter(GenA %in% shared_nodes & GenB %in% shared_nodes)

#######
#Build igraph objects
#g_AD_shared  <- graph_from_data_frame(AD_shared,  directed = FALSE)
#g_ctl_shared <- graph_from_data_frame(ctl_shared, directed = FALSE)

g_AD_shared  <- graph_from_data_frame(AD,  directed = FALSE)
g_ctl_shared <- graph_from_data_frame(ctl, directed = FALSE)


########
#Save to RDS
saveRDS(g_AD_shared,  "/STORAGE/csbig/jruiz/Redes_Pau/g_AD_200k.rds")
saveRDS(g_ctl_shared, "/STORAGE/csbig/jruiz/Redes_Pau/g_ctl_200k.rds")
saveRDS(shared_nodes, "/STORAGE/csbig/jruiz/Redes_Pau/shared_nodes.rds")
