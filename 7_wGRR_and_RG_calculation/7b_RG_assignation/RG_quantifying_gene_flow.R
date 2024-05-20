#libraries
library("data.table")
library("readxl")
library("RColorBrewer")
library("tidyverse")
library(circlize)
library(reshape2)

### Load input files
#MGEs information
MGE_tbl = fread("MGE_info_table.tsv")

#Gene annotations
RG_NRGs = fread("RG_NRG_list_25.tsv") %>%filter(Recombination_Type %in% c("RG")) #, !(Locus_Tag %in% "noTag") removed (not needed)

#RG assignments 
BBH_08_pident_01_wgrr = fread("Pipolins_MGEs_RefSeq_RGs_25.tsv") 

#RG families
RG_fams= fread("Pipolins_MGEs_RefSeq_RG_25_families.tsv")

# attach MGE and RG information on BBH tbl
BBH_08_pident_gfams = BBH_08_pident_01_wgrr %>% 
  
  #This is to keep RGs (redundant?)
  semi_join(RG_NRGs, by = c("qseqid"="Gene_ID")) %>%
  semi_join(RG_NRGs, by = c("sseqid"="Gene_ID"))%>%
  
  # assign_MGEs_to 
  mutate(MGE_from = case_when(query  %in% (MGE_tbl %>% filter(MGE == "pipolin"))$MGE_id ~ "pipolin", 
                              query  %in% (MGE_tbl %>% filter(MGE == "phage"))$MGE_id ~ "phage",
                              query  %in% (MGE_tbl %>% filter(MGE == "plasmid"))$MGE_id ~ "plasmid",
                              query  %in% (MGE_tbl %>% filter(MGE == "ciMGE"))$MGE_id ~ "ciMGE")
         ,MGE_to = case_when(subject  %in% (MGE_tbl %>% filter(MGE == "pipolin"))$MGE_id ~ "pipolin", 
                             subject  %in% (MGE_tbl %>% filter(MGE == "phage"))$MGE_id ~ "phage",
                             subject  %in% (MGE_tbl %>% filter(MGE == "plasmid"))$MGE_id ~ "plasmid",
                             subject  %in% (MGE_tbl %>% filter(MGE == "ciMGE"))$MGE_id ~ "ciMGE")) %>% 
  
  # assign RG families
  left_join(RG_fams %>% select(RG_family_from=RG_family, Member_from=Gene_ID), by = c("qseqid"="Member_from")) %>%
  left_join(RG_fams %>% select(RG_family_to=RG_family, Member_to=Gene_ID), by = c("sseqid"="Member_to")) %>%

  # remove redundant events/counts
  group_by(RG_family_from,RG_family_to, MGE_from, MGE_to) %>% 
  summarise(n_hits = n()) %>% #, mean_pident = mean(pident), mean_cov_from = mean(coverage_from), mean_cov_to = mean(coverage_to)
  ungroup()



#Count gene flow events (cluster/family level)
flow_count_MGE_gfams = BBH_08_pident_gfams %>% 
  select(RG_family_from,RG_family_to, MGE_from, MGE_to) %>% ungroup() %>%
  group_by(RG_family_from,RG_family_to) %>% mutate(counts = n(), MGE_from_to = str_c(MGE_from, MGE_to, sep = "_") ) %>% 
  ungroup() %>% dplyr::count(MGE_from, MGE_to) %>% group_by(MGE_from) %>%
   mutate(n_gene_fams = sum(n), fraction = round(n*100/n_gene_fams,1)) %>% ungroup() %>%
  mutate(n_exchanges = sum(n), fraction_exchanges = round(n*100/n_exchanges,1))


#transform to matrix
matrix_flow_MGE <- dcast(flow_count_MGE_gfams, MGE_from ~ MGE_to, value.var = "n")
row.names(matrix_flow_MGE) <- matrix_flow_MGE$MGE_from
matrix_flow_MGE$MGE_from <- NULL

# Make the circular plot
chordDiagram(as.matrix(matrix_flow_MGE), transparency = 0.5) #not very nice
chordDiagram(as.matrix(matrix_flow_MGE[3,]), transparency = 0.5) #only pipolins, this is even worse


#heatmap
ggplot(flow_count_MGE_gfams, aes(x = MGE_to, y = MGE_from, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "#FF6666") +
  theme_classic() +
  labs(x = "MGE to", y = "MGE from", title = "MGE from vs MGE to Heatmap")
