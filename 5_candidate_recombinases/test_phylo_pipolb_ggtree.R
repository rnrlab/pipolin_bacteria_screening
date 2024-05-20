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


# Read trees
pipolb_tree <- read.tree("../piPolB_analysis/pipolb_seq_filter800aa_mafft-linsi_trimAl_auto.aln.contree")
outG <- match("tr|Q6X3W4|Q6X3W4_BP35C", pipolb_tree$tip.label)
pipolb_tree <- reroot(pipolb_tree, outG, 0.5) #Root on outG (Bam35 DNAP)

#Read data -> PIPOLIN DATA
pipolin_data <- read.table("table_info_pipolins_phylo.txt", header = TRUE, sep="\t") %>%
                            group_by(order) %>%
                            mutate(order_n = ifelse(n() < 20 | order=="", "Other", as.character(order))) %>%
                            ungroup()

#Color coding
pipolin_data$order <- relevel(factor(pipolin_data$order_n), "Other")
order_other_names_vector <- levels(pipolin_data$order)
order_color_vector <- c("#3b3b3b",
                        "#8a2b32", "#fa7223", "#25d9c1", "#eb3f4d", "#c782ff", 
                        "#87ffb3", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",
                        "#5bad51", "#ffa947") 

#Read data -> GENE CLUSTERING
pipolb_data <- read.table("piPolB_context_gene_presence.txt", header = TRUE)
pipolb_data[nrow(pipolb_data) + 1,] <- c(rep(0, ncol(pipolb_data)))
pipolb_data[nrow(pipolb_data),]$piPolB_repre <- "tr|Q6X3W4|Q6X3W4_BP35C"
rownames(pipolb_data) <- pipolb_data$piPolB_repre
colnames(pipolb_data) <- gsub("Cluster_", "", colnames(pipolb_data))
pipolb_data_matrix <- as.matrix(pipolb_data[,2:ncol(pipolb_data)])

#Filter only diff freq > X
pipolb_associated_clusters <- read.table("../Cluster_association/list_pipolb_assoc_rec_clusters_list.txt", header = TRUE)
pipolb_data_matrix_pipolb_associated <- pipolb_data_matrix[,unique(pipolb_associated_clusters$Index)]
pipolb_data_matrix_pipolb_associated_clu1 <- pipolb_data_matrix[,unique(pipolb_associated_clusters[pipolb_associated_clusters$Near_piPolB=="Cluster_1",]$Index)]
pipolb_data_matrix_pipolb_associated_noclu1 <- pipolb_data_matrix[,unique(pipolb_associated_clusters[pipolb_associated_clusters$Near_piPolB!="Cluster_1",]$Index)]

#Color vector for heatmap (heatmap is for order)
pipolin_data$order_col <- pipolin_data$order
levels(pipolin_data$order_col) <- order_color_vector
color_row_vector <- c(as.vector(pipolin_data$order_col),"#3b3b3b") #Last is phi29 (should not be here)

pipolb_data_matrix_heatmap <- heatmap(pipolb_data_matrix_pipolb_associated_clu1, scale="none") #Rowv = NA, , RowSideColors=color_row_vector
pipolb_data_matrix_heatmap_ord <- pipolb_data_matrix_pipolb_associated_clu1[,pipolb_data_matrix_heatmap$colInd]

pipolb_data_matrix_heatmap_Gcol <- heatmap(pipolb_data_matrix_pipolb_associated_noclu1, scale="none")
pipolb_data_matrix_heatmap_Gcol_ord <- pipolb_data_matrix_pipolb_associated_noclu1[,pipolb_data_matrix_heatmap_Gcol$colInd]
pipolb_data_matrix_heatmap_Gcol_ord_custom <- pipolb_data_matrix[,c(99,1165,470,1018,611,929,
                                                                    2047,705,1152,1212,218,1643,
                                                                    398,474,67,767,1082,1251,
                                                                    388,988,433,719)]         

#Read data -> PROPHAGES AND ASSEMBLY GAPS
#prophage_data <-  read.table("pipolin_provirus_data.txt", header = TRUE)
#prophage_data$Assembly_gaps <- as.factor(prophage_data$Assembly_gaps)

#Plot basic tree
p <- ggtree(pipolb_tree, size=0.2) #, layout="fan", open.angle=5
#print node id? ->  geom_text(aes(label=node), size=1, hjust=-.3)

offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
p_Gcol <- scaleClade(p, 9751, .001) %>% collapse(9747, 'min', fill="gray")   


# Add order info
p <- p %<+% pipolin_data +
  geom_tippoint(aes(color=order), size=0.05, alpha=0.6) + 
  scale_color_manual(values=order_color_vector) + ggtree::vexpand(.1, -1)

p_Gcol <- p_Gcol %<+% pipolin_data +
  geom_tippoint(aes(color=order), size=0.05, alpha=0.6) + 
  scale_color_manual(values=order_color_vector)

#Plot clusters
p_plot <- gheatmap(p, pipolb_data_matrix_pipolb_associated_clu1, color = NA, width = 0.75, 
                   colnames = FALSE, hjust=6) + 
              scale_fill_continuous(high="#363636",low="white")+  # , colnames_offset_y = 0, hjust=1
              theme(legend.position = "none") +
              geom_treescale(x = -0.5, y = 0.5) + 
              coord_flip() 
              
p_plot 
ggsave(file="phylogeny_cluster1_wRec.svg", plot=p_plot, width=10, height=3.5)


p_Gcol_plot <- gheatmap(p_Gcol, pipolb_data_matrix_heatmap_Gcol_ord_custom, color = NA, width = 5,
                        colnames = FALSE, offset = 2) + 
                   scale_fill_continuous(high="#363636",low="white") +
                   theme(legend.position = "none") +
                   geom_treescale(x = 5, y = 0.5) + 
                   coord_flip()
p_Gcol_plot 
ggsave(file="phylogeny_clusterOther_wRec.svg", plot=p_Gcol_plot, width=10, height=5)









### DEPRECATED ####

# ### Alternative:
# p <- p + new_scale_fill() +
#   geom_fruit(data=pipolb_data_matrix_melted, geom=geom_tile,
#              mapping=aes(y=Var1, x=Var2, alpha=value, fill=Var2),
#              color = NA, offset = 0.02,size = 0.02) + ggtree::vexpand(.1, -1)

#Add outer bars - viral genes
p <- p + new_scale_fill() + 
         geom_fruit(data=prophage_data, geom=geom_bar,
                    mapping=aes(y=piPolB_repre, x=Viral_genes_in_pipolin, fill=Assembly_gaps),
                    pwidth=0.2, 
                    orientation="y", 
                    stat="identity") +
         scale_fill_manual(values=c("#32d940","#ffea00",rep("#d94123",length(levels(prophage_data$Assembly_gaps))-2)))


svglite("pipolb_tree_800aa.svg", width = 50, height = 2000)
p + layout_rectangular() + theme(legend.position="none")  + geom_tiplab(offset = .1, hjust = .1, size = 1) 
dev.off()

write.table(get_taxa_name(p),file="pipolin_repre.txt", row.names = FALSE, col.names = FALSE, quote=FALSE) ### This is for ggenomes



### Plotting other stuff in tree (Atts, Integration site)
p <- ggtree(pipolb_tree, size=0.1, layout="fan", open.angle=5) #, layout="fan", open.angle=5
p <- p %<+% pipolin_data +
  geom_tippoint(aes(color=order), size=0.2, alpha=0.8) + 
  scale_color_manual(values=order_color_vector) + ggtree::vexpand(.1, -1) 

# add integration site
p <- p + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=integration_site),
    width=0.05,
    offset=0.1
  )

#Add att num info
p <- p + new_scale_fill() +  geom_fruit(geom=geom_bar,
             mapping=aes(y=piPolB_repre, x=atts),
             pwidth=0.2, 
             orientation="y", 
             stat="identity")

p <- p +  geom_fruit(geom=geom_bar,
                     mapping=aes(y=piPolB_repre, x=piPolB_len_sum),
                     pwidth=0.2, 
                     orientation="y", 
                     stat="identity")

p + layout_rectangular() + geom_tiplab(offset = .1, hjust = .1, size = 1) + theme(legend.position="none") 

### Plot PADLOC HMM data
padloc_hmm_pipolins <- read.table("pipolin_padloc_cdhit.txt", header=TRUE)
padloc_hmm_pipolins[nrow(padloc_hmm_pipolins) + 1,] <- c(rep(0, ncol(padloc_hmm_pipolins)))
padloc_hmm_pipolins[nrow(padloc_hmm_pipolins),]$piPolB_repre <- "tr|Q6X3W4|Q6X3W4_BP35C"
rownames(padloc_hmm_pipolins) <- padloc_hmm_pipolins$piPolB_repre
colnames(padloc_hmm_pipolins) <- gsub("Cluster_", "", colnames(padloc_hmm_pipolins))
padloc_hmm_pipolins_matrix <- as.matrix(padloc_hmm_pipolins[,2:ncol(padloc_hmm_pipolins)])
padloc_hmm_pipolins_matrix_subset <- padloc_hmm_pipolins_matrix[,colSums(padloc_hmm_pipolins_matrix)>=5]
padloc_hmm_pipolins_matrix_subset_heatmap <-  heatmap(padloc_hmm_pipolins_matrix_subset, scale="none", RowSideColors=color_row_vector) #, Rowv = NA
padloc_hmm_pipolins_matrix_subset_heatmap_ordered <- padloc_hmm_pipolins_matrix_subset[,padloc_hmm_pipolins_matrix_subset_heatmap$colInd]

p <- ggtree(pipolb_tree, size=0.1, layout="fan", open.angle=5) #, layout="fan", open.angle=5
p <- p %<+% pipolin_data +
  geom_tippoint(aes(color=order), size=0.2, alpha=0.8) + 
  scale_color_manual(values=order_color_vector) + ggtree::vexpand(.5, -1) 
p <- gheatmap(p, padloc_hmm_pipolins_matrix_subset_heatmap_ordered, color = NA, width = 5, 
              colnames_angle=90, font.size = 2, hjust=1) + scale_fill_continuous(high="#c27108",low="#f0f0f0") 

p + layout_rectangular() 

### Plot pass genes
cluster_info <- read.table("../pipolin_protein_clustering/Cluster_info_simplified.txt", header=TRUE, sep="\t")
pass_indices <- as.numeric(rownames(cluster_info[cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_pass <- pipolb_data_matrix[,pass_indices]
pipolb_data_matrix_pipolb_pass_filtered <- pipolb_data_matrix_pipolb_pass[,colSums(pipolb_data_matrix_pipolb_pass) > 10]

p <- ggtree(pipolb_tree, size=0.4, layout="fan") #, layout="fan", open.angle=5
p <- p %<+% pipolin_data +
  geom_tippoint(aes(color=order), size=0.5, alpha=0.8) + 
  scale_color_manual(values=order_color_vector) + ggtree::vexpand(.1, -1)
p_rec <- gheatmap(p, pipolb_data_matrix_pipolb_pass_filtered, color = NA, width = 5, 
                  colnames_angle=90, font.size = 8, hjust=1) + scale_fill_continuous(high="#572203",low="#e8e5e3")  # , colnames_offset_y = 0, hjust=1
svglite("pipolb_tree_800aa_with_Pass.svg", width = 50, height = 80)
p_rec + layout_rectangular() + theme(legend.position="none") # + geom_tiplab(offset = .1, hjust = .1, size = 5) 
dev.off()

