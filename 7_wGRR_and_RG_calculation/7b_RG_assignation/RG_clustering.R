####
library(data.table)
library(tidyverse)
library(igraph)

#Load RGs
RGdf = fread("Pipolins_MGEs_RefSeq_RGs_25.tsv",select=c("qseqid", "sseqid"))

#Get graph objects
edge_tbl = RGdf %>% select(qseqid, sseqid) 
vert_tbl = RGdf %>% select(Gene_ID=qseqid) %>% 
                     bind_rows(RGdf%>% ungroup %>% select(Gene_ID=sseqid)) %>% distinct()

#Make graph
gene_network = graph_from_data_frame(d = edge_tbl, vertices = vert_tbl)

#Cluster them if they are connected (single linkage)
components = data.frame( cluster = components(gene_network, mode = "weak")$membership) %>% rownames_to_column("protein_id")

#Make gene families
clustered_gene = components %>% 
  left_join(components %>% group_by(cluster) %>% slice(1) %>% select(cluster, gene_family = protein_id), by = "cluster")%>%
  select(RG_family=gene_family, Gene_ID = protein_id) %>% group_by(RG_family) %>% mutate(n_members = n()) %>% ungroup()

#Add vert MGE type data
MGE_type = fread("MGE_info_table.tsv")
clustered_gene_names_info = clustered_gene %>% mutate(MGE_id = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>%
                                                left_join(MGE_type, by=c("MGE_id"))%>%
                                                group_by(RG_family) %>% 
                                                mutate(in_MGE = paste(unique(MGE), collapse = ", ")) %>%
                                                select(-MGE_id)

write_tsv(clustered_gene_names_info, "Pipolins_MGEs_RefSeq_RG_25_families.tsv")
