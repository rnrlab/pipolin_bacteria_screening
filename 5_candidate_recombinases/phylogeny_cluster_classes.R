library(phytools)
library(ggtreeExtra)
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

# 1) Read tree
pipolb_tree <- read.tree("pipolb_seq_filter800aa_mafft-linsi_trimAl_auto.aln.contree")
outG <- match("tr|Q6X3W4|Q6X3W4_BP35C", pipolb_tree$tip.label)
pipolb_tree <- reroot(pipolb_tree, outG, 0.5) #Root on outG (Bam35 DNAP)

# 2) Read pipolin data
pipolin_data <- read.table("table_info_pipolins_phylo.txt", header = TRUE, sep="\t")

#Color coding
pipolin_data$genus <- relevel(factor(pipolin_data$genus), "Other")
genus_other_names_vector <- levels(pipolin_data$genus)
genus_color_vector <- c("#3b3b3b",
                        "#9221fc", "#ff920d", "#fc58e4", "#9e6e4a", "#de4040", 
                        "#69d3fa", "#ffa8ab", "#86e374", "#4781e6") #Pseudosulfitobacter ("#10780b") removed

# 3) Read clustering data
pipolb_data <- read.table("piPolB_context_gene_presence.txt", header = TRUE)
pipolb_data[nrow(pipolb_data) + 1,] <- c(rep(0, ncol(pipolb_data)))
pipolb_data[nrow(pipolb_data),]$piPolB_repre <- "tr|Q6X3W4|Q6X3W4_BP35C" #add outgroup name
rownames(pipolb_data) <- pipolb_data$piPolB_repre
colnames(pipolb_data) <- gsub("Cluster_", "", colnames(pipolb_data))
pipolb_data_matrix <- as.matrix(pipolb_data[,2:ncol(pipolb_data)])

# 4) Plot  categories
cluster_info <- read.table("../pipolin_protein_clustering/Cluster_info_simplified.txt", header=TRUE, sep="\t")
pass_indices <- as.numeric(rownames(cluster_info[cluster_info$Pass=="Pass",]))

DUF_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="DUF and UPF"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_DUF <- pipolb_data_matrix[,c(DUF_indices)]
DUF_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_DUF)>=1)
names(DUF_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_DUF)>=1)

EXC_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Excisionase"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_EXC <- pipolb_data_matrix[,c(EXC_indices)]
EXC_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_EXC)>=1)
names(EXC_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_EXC)>=1)

AAA_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Helicase-AAA"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_AAA <- pipolb_data_matrix[,c(AAA_indices)]
AAA_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_AAA)>=1)
names(AAA_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_AAA)>=1)

Transp_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Membrane"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_Transp <- pipolb_data_matrix[,c(Transp_indices)]
Transp_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_Transp)>=1)
names(Transp_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_Transp)>=1)

DNAres_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="NA Cleavage"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_DNAres <- pipolb_data_matrix[,c(DNAres_indices)]
DNAres_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_DNAres)>=1)
names(DNAres_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_DNAres)>=1)

DNAmet_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="NA Modification"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_DNAmet <- pipolb_data_matrix[,c(DNAmet_indices)]
DNAmet_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_DNAmet)>=1)
names(DNAmet_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_DNAmet)>=1)

NOP_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="No PFAM"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_NOP <- pipolb_data_matrix[,c(NOP_indices)]
NOP_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_NOP)>=1)
names(NOP_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_NOP)>=1)

NUCMet_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Nucleotide metabolism"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_NUCMet <- pipolb_data_matrix[,c(NUCMet_indices)]
NUCMet_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_NUCMet)>=1)
names(NUCMet_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_NUCMet)>=1)

Other_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Other"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_Other <- pipolb_data_matrix[,c(Other_indices)]
Other_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_Other)>=1)
names(Other_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_Other)>=1)

OtherDNAbind_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Other DNA binding"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_OtherDNAbind <- pipolb_data_matrix[,c(OtherDNAbind_indices)]
OtherDNABind_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_OtherDNAbind)>=1)
names(OtherDNABind_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_OtherDNAbind)>=1)


OtherHydrol_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Other hydrolases"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_OtherHydrol <- pipolb_data_matrix[,c(OtherHydrol_indices)]
OtherHydrol_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_OtherHydrol)>=1)
names(OtherHydrol_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_OtherHydrol)>=1)

Peptidase_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Peptidase"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_Peptidase <- pipolb_data_matrix[,c(Peptidase_indices)]
Peptidase_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_Peptidase)>=1)
names(Peptidase_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_Peptidase)>=1)

Rel_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Relaxase"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_REL <- pipolb_data_matrix[,c(Rel_indices)]
Rel_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_REL)>=1)
names(Rel_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_REL)>=1)

Srec_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Srec"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_SREC <- pipolb_data_matrix[,c(Srec_indices)]
SRec_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_SREC)>=1)
names(SRec_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_SREC)>=1)

Tnp_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Transposase"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_Tnp <- pipolb_data_matrix[,c(Tnp_indices)]
Tnp_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_Tnp)>=1)
names(Tnp_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_Tnp)>=1)

Yrec_indices <- as.numeric(rownames(cluster_info[cluster_info$Custom_category=="Yrec"&cluster_info$Pass=="Pass",]))
pipolb_data_matrix_pipolb_associated_YREC <- pipolb_data_matrix[,c(Yrec_indices)]
YRec_presence <- as.integer(rowSums(pipolb_data_matrix_pipolb_associated_YREC)>=1)
names(YRec_presence) <- names(rowSums(pipolb_data_matrix_pipolb_associated_YREC)>=1)

Clases_presence <- data.frame("YRec" = YRec_presence, "Exc" = EXC_presence, "SRec" = SRec_presence, "Relaxase" = Rel_presence, "IS/Tnp" = Tnp_presence,
                              "DNA bind other" = OtherDNABind_presence, "AN modification" = DNAmet_presence, "AN cleaveage" = DNAres_presence, 
                              "Nuc Metab" = NUCMet_presence, "Peptidase" = Peptidase_presence, "Hydrolases other" = OtherHydrol_presence,
                              "Helicase-AAA" = AAA_presence, "Membrane" = Transp_presence,
                              "DUF-UPF" = DUF_presence, "No PFAM" = NOP_presence, "Other" = Other_presence)

Class_test <- data.frame("YRec" = YRec_presence)

p <- ggtree(pipolb_tree, size=0.4, layout="fan") #, layout="fan", open.angle=5
p <- p %<+% pipolin_data +
  geom_tippoint(aes(color=genus), size=0.5, alpha=0.8) + 
  scale_color_manual(values=genus_color_vector) + ggtree::vexpand(.1, -1)

p_classes <- gheatmap(p, Clases_presence, color = NA, width = 5, 
                      colnames_angle=90, font.size = 2, hjust=1) + scale_fill_continuous(high="#572203",low="#e8e5e3")  # , colnames_offset_y = 0, hjust=1
svglite("pipolb_tree_800aa_with_Clases_CIRC.svg", width = 10, height = 10)
p_classes + theme(legend.position="none")  + layout_rectangular() # + geom_tiplab(offset = .1, hjust = .1, size = 5) 
dev.off()




