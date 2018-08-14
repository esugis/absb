library(data.table)
#library(tidyverse)
library(dplyr)
library(igraph)
library(tidyr)

dt_nodes <- fread("~/absb/data/graph/nodes_alz_nonalz_gwas_filt.txt", header=T)
#head(dt_nodes)
dt_interactions <- fread("~/absb/data/graph/interactions_alz_nonalz_gwas.txt", header=T)

cleaned_nodes <- dt_nodes %>%
  group_by(ensg, node_type) %>%
  tally %>%
  spread(node_type, n, fill= 0) %>%
  mutate(total = alz + nonalz) %>%
  mutate(node_type = ifelse(alz!=0, 1, 0))
#mutate(node_type = ifelse(nonalz > alz, 0, 1))

# how to deal with features, mostly features duplicated apart fromn snp_id
dt_nodes$SNP_id <- NULL
dt_nodes$GWAS_pvalue <- NULL
dt_nodes$node_type <- NULL

dt_nodes <- distinct(dt_nodes) # 46519 down to 36004
dt_nodes <- left_join(cleaned_nodes, dt_nodes, by='ensg') # cleaned nodes!

dt_interactions$data_source <- NULL
dt_interactions$score <- NULL
dt_interactions$alz_edge <- NULL

dt_interactions <- distinct(dt_interactions) # 122358 down to 116177
g <- graph_from_data_frame(dt_interactions, directed = F, vertices=dt_nodes$ensg)

alz_types = dt_nodes[,c(1,5)]
two_hop_neigh <- function(x) {
  bfs_inter <- bfs(g, root=x, order=TRUE, rank=TRUE, dist=TRUE)$dist
  filtered_dist <- bfs_inter[bfs_inter==2]
  node_names <- names(filtered_dist)
  
  alz_count = 0
  for(item in node_names){
    alz_count = alz_count + alz_types$node_type[alz_types$ensg == item]
  }
  
  ratio = alz_count / length(node_names)
  print(item)
  return(c(alz_count, ratio))
}

dt_features <- alz_types %>%
  mutate(hop_count = two_hop_neigh(ensg)[1], hop_ratio=two_hop_neigh(ensg)[2])
print("Calculations are done!")

write.table(dt_features, "~/absb/results/graph/two_hop_alz.csv", col.names = T, row.names = F)
