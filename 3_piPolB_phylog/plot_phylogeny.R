library(ggtree)
library(phyloseq)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(treeio)
library(tidytree)
library(ggstar)
library(ggnewscale)
library(svglite)
library(phytools)
packageVersion("phytools")

# Read tree
pipolb_tree <- read.tree("pipolb_seq_filter800aa_mafft-linsi_trimAl_auto.aln.contree")
outG <- match("tr|Q6X3W4|Q6X3W4_BP35C", pipolb_tree$tip.label)
pipolb_tree <- reroot(pipolb_tree, outG, 0.5) #Root on outG (Bam35 DNAP)

#Read data -> PIPOLIN DATA
pipolin_data_extra <- read.table("../pipolin_summary_new.tsv", header = TRUE, sep="\t") %>%
  group_by(genus) %>%
  mutate(genus_n = ifelse(n() < 10 | genus=="", "Other", as.character(genus))) %>%
  ungroup()

pipolin_data <- read.table("../pipolin_phylogeny_correlations/table_info_pipolins_phylo.txt", header = TRUE, sep="\t") %>%
               left_join(pipolin_data_extra, by=c("pipolin"="Pipolin_ID")) %>%
               mutate(integration_site_grouped = ifelse(grepl("SsrA", integration_site), "SsrA tmRNA", ifelse(grepl(",", integration_site)| !grepl("tRNA", integration_site), NA, integration_site)),
                      Att_type_grouped = ifelse(Att_type=="[de novo]", "De novo", ifelse(Att_type=="[pipolin conserved]", "E. coli like (BLASTn)", NA))) %>%
               select(piPolB_repre, pipolin, order=order.x, genus_n, atts, Att_type, Att_type_grouped, integration_site, integration_site_grouped, 
                      piPolB_len_sum, Accession, species)


#Color coding
pipolin_data$genus <- relevel(factor(pipolin_data$genus_n), "Other")
order_other_names_vector <- levels(pipolin_data$genus)
order_color_vector <- c("#3b3b3b", "#8a2b32", "#8e2abf", "#9c2c6b", "#f5828f", 
                        "#fa7223", "#ff94f1", "#ed5611", "#a83258", "#eb3f4d",
                        "#c97190", "#9e0606", "#e32bac", "#4a47e6", "#4781e6",
                        "#83a838","#5bad51","#f0ab8b","#e3c371","#87ffb3",
                        "#47c483","#c782ff","#d41717","#6b1212","#25d9c1",
                        "#ffa947") 
#plot tree
p <- ggtree(pipolb_tree, size=0.2) #, layout="fan", open.angle=5
p_circ <- ggtree(pipolb_tree, size=0.2, layout="fan", open.angle=5)

p_info <- p_circ %<+% pipolin_data +
  geom_tippoint(aes(color=genus), size=0.02, alpha=0.6) + 
  scale_color_manual(values=order_color_vector) + ggtree::vexpand(.1, -1) + 
  theme(legend.position="none")  + geom_tiplab(aes(label=paste(label, pipolin, Accession, species, sep="   ")),
                                                   size = 0.2) +
  geom_nodelab(aes(label = label), nudge_x = -0.03, nudge_y = 0.5, size = 0.2) 

svglite("Fig_S_pipolb_tree_800aa_with_info_circ.svg", width = 50, height = 50)
p_info
dev.off()

#integration site info
int_df <- data.frame(Integration_site=factor(pipolin_data$integration_site_grouped), Type=factor(pipolin_data$Att_type_grouped)) 
rownames(int_df) <- pipolin_data$piPolB_repre


#svglite("Fig_S_pipolb_tree_800aa_with_info_IntSite.svg", width = 50, height = 50)
gheatmap(p_info, int_df, color=NA, offset=0.1, width=.1, colnames = FALSE) +
  scale_fill_manual(values=c("coral","lightblue",
                             "#0E294B","#641C34","#F75E25","#89AC76","#C51D34","#2271B3","#00BB2D",
                             "#BDECB6","#6C4675","#3E5F8A","brown1","#955F20","#7FB5B5","#F39F18",
                             "#9D9101"))
#dev.off()

#alter: linear
p_info_lin <- p %<+% pipolin_data +
  geom_tippoint(aes(color=genus), size=0.02, alpha=0.6) + 
  scale_color_manual(values=order_color_vector) + ggtree::vexpand(.1, -1) + 
  theme(legend.position="none")  + geom_tiplab(aes(label=paste(label, pipolin, Accession, species, sep="   ")),
                                               size = 0.8) +
  geom_nodelab(aes(label = label), nudge_x = -0.03, nudge_y = 0.5, size = 0.5) 