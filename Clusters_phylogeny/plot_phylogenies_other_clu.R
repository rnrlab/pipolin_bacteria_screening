#library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(phytools)
library(treeio)
library(tidytree)
library(ggstar)
library(ggnewscale)
library(svglite)
library(gridExtra)

#Read data -> PIPOLIN DATA
pipolin_data <- read.table("../../pipolin_summary_new.tsv", header = TRUE, sep="\t") %>%
  group_by(order) %>%
  mutate(order_n = ifelse(n() < 20 | order=="", "Other", as.character(order))) %>%
  ungroup() %>%
  group_by(genus) %>%
  mutate(genus_n = ifelse(n() < 15 | genus=="", "Other", as.character(genus))) %>%
  ungroup()

pipolin_data$genus <- relevel(factor(pipolin_data$genus_n), "Other")
genus_other_names_vector <- levels(pipolin_data$genus)
genus_color_vector <- c("#3b3b3b",
                        "chartreuse2", "#25d9c1","yellow", "#9d47ff","#fa1919",
                        "#4781e6", "#ff7ad7", "darkgreen", "darkorange2",
                        "#3b3b3b")#",


### Cluster 2 - IntSXT ###
Clu2_tree <- read.tree("Cluster_2_IntSXT/Cluster_2_out_mafft-linsi_trimal.aln.contree")
outG <- match("sp|P22886|TNR6_ENTFL", Clu2_tree$tip.label)
Clu2_tree_rooted <- reroot(Clu2_tree, outG, 2) 
p_Clu2 <- ggtree(Clu2_tree_rooted, size=0.2, layout = "fan", open.angle=1)

Clu2_info <- data.frame(Clu2_id = Clu2_tree_rooted$tip.label,
                       Pipolin_ID = sub("_[^_]*$", "", Clu2_tree_rooted$tip.label)) %>%
             left_join(pipolin_data, by=c("Pipolin_ID"))

genus_color_vector_Clu2 <- c("#3b3b3b",
                             "chartreuse2", "#25d9c1","yellow", "#9d47ff","#fa1919",
                             "#4781e6", "#ff7ad7", "darkgreen", "darkorange2",
                             "#3b3b3b")#", 

p_Clu2_plot <- p_Clu2 %<+% Clu2_info +
                      geom_tippoint(aes(color=genus), size=0.02, alpha=0.4)+ 
                      scale_color_manual(values=genus_color_vector_Clu2) +
                      geom_treescale(x = -1)
p_Clu2_plot

### Cluster 3 - DUF2787 ### 
Clu3_tree <- read.tree("Cluster_3_DUF2787/Cluster_3_out_mafft-linsi_trimal.aln.contree")
Clu3_outG <- match("tr|Q9KUL0|Q9KUL0_VIBCH", Clu3_tree$tip.label)
Clu3_tree_rooted <- reroot(Clu3_tree, Clu3_outG, 0.5) 
p_Clu3 <- ggtree(Clu3_tree_rooted, size=0.2, layout = "fan", open.angle=1)

Clu3_info <- data.frame(Clu3_id = Clu3_tree_rooted$tip.label,
                        Pipolin_ID = sub("_[^_]*$", "", Clu3_tree_rooted$tip.label)) %>%
  left_join(pipolin_data, by=c("Pipolin_ID"))

genus_color_vector_Clu3 <- c("#3b3b3b",
                              "chartreuse2", "#25d9c1","yellow", "#9d47ff","#fa1919",
                              "#4781e6", "#ff7ad7", "darkgreen", "darkorange2",
                              "#3b3b3b")#",  

p_Clu3_plot <- p_Clu3 %<+% Clu3_info +
  geom_tippoint(aes(color=genus), size=0.02, alpha=0.4)+ 
  scale_color_manual(values=genus_color_vector_Clu3 ) +
  geom_treescale(x = -1)
p_Clu3_plot

### Cluster 5 - UDG ###
Clu5_tree <- read.tree("Cluster_5_UDG/Cluster_5_out_mafft-linsi_trimal.aln.contree")
Clu5_outG <- match("CDK41219.1", Clu5_tree$tip.label)
Clu5_tree_rooted <- reroot(Clu5_tree, Clu5_outG, 0.5) 
p_Clu5 <- ggtree(Clu5_tree_rooted, size=0.2, layout = "fan", open.angle=1)

Clu5_info <- data.frame(Clu5_id = Clu5_tree_rooted$tip.label,
                        Pipolin_ID = sub("_[^_]*$", "", Clu5_tree_rooted$tip.label)) %>%
  left_join(pipolin_data, by=c("Pipolin_ID"))

genus_color_vector_Clu5 <- c("#3b3b3b",
                              "#25d9c1","yellow", "#9d47ff","#fa1919",
                              "#ff7ad7", "darkgreen")

p_Clu5_plot <- p_Clu5 %<+% Clu5_info +
  geom_tippoint(aes(color=genus), size=0.00001, alpha=0.4)+ 
  scale_color_manual(values=genus_color_vector_Clu5 ) +
  geom_treescale(x = -1)
p_Clu5_plot

### Cluster 11 - IntP2 ###  
Clu11_tree <- read.tree("Cluster_11_IntP2/Cluster_11_out_mafft-linsi_trimal.aln.contree")
outG <- match("sp|P22886|TNR6_ENTFL", Clu11_tree$tip.label)
Clu11_tree_rooted <- reroot(Clu11_tree, outG, 0.5) 
p_Clu11 <- ggtree(Clu11_tree_rooted, size=0.2, layout = "fan", open.angle=1)

Clu11_info <- data.frame(Clu11_id = Clu11_tree_rooted$tip.label,
                        Pipolin_ID = sub("_[^_]*$", "", Clu11_tree_rooted$tip.label)) %>%
  left_join(pipolin_data, by=c("Pipolin_ID"))

genus_color_vector_Clu11 <- c("#3b3b3b",
                             "chartreuse2", "#25d9c1","yellow", "#9d47ff","#fa1919",
                             "#4781e6", "#ff7ad7", "darkgreen", "darkorange2",
                             "#3b3b3b")#", 

p_Clu11_plot <- p_Clu11 %<+% Clu11_info +
  geom_tippoint(aes(color=genus), size=0.00001, alpha=0.4)+ 
  scale_color_manual(values=genus_color_vector_Clu11) +
  geom_treescale(x = -1) 
p_Clu11_plot




### Merge plots
full_plot <- grid.arrange(p_Clu2_plot,p_Clu3_plot,p_Clu5_plot, p_Clu11_plot,nrow=2)
ggsave(file="Cluster_phylogenies.svg", plot=full_plot, width=50, height=50, limitsize = FALSE)
