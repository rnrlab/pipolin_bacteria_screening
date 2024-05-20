# libraries
library("data.table")
library("RColorBrewer")
library("tidyverse")
library("rstatix")
library("ggstatsplot")
library("ggplot2")
library(randomcoloR)
library("gridExtra")
library(cowplot)
library("svglite")

## load in pipolin, ciMGEs, plasmid and phage tbl (V_ not needed)
MGE_tbl = fread("MGE_info_table.tsv")

# load in RGs, NRGs, nhNRGs (wGRR < 0.1 and < 0.2)
RG_NRGs =  fread("RG_NRG_list.tsv") %>% mutate(MGE_id = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>% left_join(MGE_tbl, by="MGE_id")
RG_NRGs_20  =  fread("RG_NRG_list_20.tsv") %>% mutate(MGE_id = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>% left_join(MGE_tbl, by="MGE_id")
RG_NRGs_25  =  fread("RG_NRG_list_25.tsv") %>% mutate(MGE_id = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>% left_join(MGE_tbl, by="MGE_id")
RG_NRGs_30  =  fread("RG_NRG_list_30.tsv") %>% mutate(MGE_id = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>% left_join(MGE_tbl, by="MGE_id")

# load in pipolin info
pipolin_info = fread("../Gene_annotation/pipolin_summary_new.tsv")
colnames(pipolin_info)[2] <- "MGE_id" #may be needed to change 3 to 2

#Plot ratio RGs
RG_NRGs_filtered = RG_NRGs %>% filter(Recombination_Type != "nhNRG") 
RG_NRGs_filtered_20 = RG_NRGs_20 %>% filter(Recombination_Type != "nhNRG")
RG_NRGs_filtered_25 = RG_NRGs_25 %>% filter(Recombination_Type != "nhNRG")
RG_NRGs_filtered_30 = RG_NRGs_25 %>% filter(Recombination_Type != "nhNRG")

ratio_RGs <- RG_NRGs_filtered %>% group_by(MGE_id) %>%
                                  summarize(ratioRG = sum(Recombination_Type == "RG") / n(), MGE = first(MGE))
ratio_RGs_20 <- RG_NRGs_filtered_20 %>% group_by(MGE_id) %>%
                                  summarize(ratioRG = sum(Recombination_Type == "RG") / n(), MGE = first(MGE))
ratio_RGs_25 <- RG_NRGs_filtered_25 %>% group_by(MGE_id) %>%
                                  summarize(ratioRG = sum(Recombination_Type == "RG") / n(), MGE = first(MGE))
ratio_RGs_30 <- RG_NRGs_filtered_30 %>% group_by(MGE_id) %>%
                                  summarize(ratioRG = sum(Recombination_Type == "RG") / n(), MGE = first(MGE))

plt <- ggbetweenstats(
  data = ratio_RGs,
  x = MGE,
  y = ratioRG
) 
plt

rg_means <- aggregate(ratioRG ~  MGE, ratio_RGs_20, mean) %>% 
  mutate(y_val = aggregate(ratioRG ~  MGE, ratio_RGs_20, max))
plt_20 <- ggbetweenstats(
  data = ratio_RGs_20,
  x = MGE,
  y = ratioRG,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0),
  violin.args = list(width = 0, linewidth = 0)
) + ggplot2::scale_color_manual(values = c("#FFE66D","#0A014F","#FF6B6B","#a39022")) + ylim(0,1) + ggtitle("RG ratio (MGEs)") +
  labs(y= "RG/(RG+NRG)", x = "MGE") + theme_classic() +
  scale_x_discrete(limits=c("phage", "ciMGE","plasmid", "pipolin")) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust=0.5),legend.position = "None") + labs(x = NULL) +
  geom_text(data = rg_means, aes(label = paste0(round(ratioRG,2)), y = y_val$ratioRG + 0.05), color="gray40") +
  annotate(geom="text", x = 3, y = 0.95, label=paste0(round(rg_means$ratioRG[4],2)), color="gray40") +
  annotate(geom="text", x = 4, y = 0.95, label=paste0(round(rg_means$ratioRG[3],2)), color="gray40") 

plt_20

plt_25 <- ggbetweenstats(
  data = ratio_RGs_25,
  x = MGE,
  y = ratioRG
)
plt_25

grid.arrange(plt+ggtitle("wGRR < 0.1"), 
             plt_20+ggtitle("wGRR < 0.2"), 
             plt_25+ggtitle("wGRR < 0.25"), ncol=3)

### By groups (taxa = order) ###

ratio_RGs_pipolin = ratio_RGs %>% filter(MGE=="pipolin") %>% left_join(pipolin_info, by="MGE_id") %>% 
                                  group_by(order) %>% mutate(order_num = n()) %>% ungroup()
ratio_RGs_20_pipolin = ratio_RGs_20 %>% filter(MGE=="pipolin") %>% left_join(pipolin_info, by="MGE_id") %>% 
                                        group_by(order) %>% mutate(order_num = n()) %>% ungroup()
ratio_RGs_25_pipolin = ratio_RGs_25 %>% filter(MGE=="pipolin") %>% left_join(pipolin_info, by="MGE_id") %>% 
                                        group_by(order) %>% mutate(order_num = n()) %>% ungroup()

color_pipolin_groups <- c("#8a2b32", "#fa7223", "#25d9c1", "#eb3f4d", "#c782ff", 
                          "#87ffb3", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",
                          "#5bad51", "#ffa947") 

plt_pipolin <- ggbetweenstats(
  data = ratio_RGs_pipolin %>% filter(order_num >= 20),
  x = order,
  y = ratioRG,
  pairwise.display = "none"
) + ggplot2::scale_color_manual(values = color_pipolin_groups)
plt_pipolin 

rg_p_means <- aggregate(ratioRG ~  order, ratio_RGs_20_pipolin %>% filter(order_num >= 20 & order != ""), mean) %>% 
  mutate(y_val = aggregate(ratioRG ~  order, ratio_RGs_20_pipolin %>% filter(order_num >= 20 & order != ""), max))
order_n <- ratio_RGs_20_pipolin %>% filter(order_num >= 20 & order != "") %>% group_by(order) %>% summarise(n = n())
x_label_names <- paste0(order_n$order, "\n(n = ",order_n$n,")")
names(x_label_names) <- order_n$order

plt_pipolin_20 <- ggbetweenstats(
  data = ratio_RGs_20_pipolin %>% filter(order_num >= 20 & order != ""),
  x = order,
  y = ratioRG,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0),
  violin.args = list(width = 0, linewidth = 0)
) + ggplot2::scale_color_manual(values = color_pipolin_groups) + ylim(0,1) +
  ggtitle("RG ratio (pipolins)") +
  labs(y= "RG/(RG+NRG)", x = NULL) + theme_classic() +
  geom_text(data = rg_p_means, aes(label = paste0(round(ratioRG,2)), y = y_val$ratioRG + 0.05), color="gray40") +
  scale_x_discrete(limits=c("Enterobacterales", "Vibrionales", "Aeromonadales", "Alteromonadales",
                            "Hyphomicrobiales", "Rhodobacterales", 
                            "Micrococcales", "Mycobacteriales",
                            "Lachnospirales", "Eubacteriales",
                            "Lactobacillales", "Bacillales"),
                   labels=x_label_names) + 
  theme(axis.text.x = element_text(angle = 45, size=10, vjust=0.5), legend.position = "None") +
  annotate(geom="text", x = 1, y = 0.95, label=paste0(round(rg_p_means$ratioRG[4],2)), color="gray40")

plt_pipolin_20

grid.arrange(plt_20, plt_pipolin_20, ncol=2, widths = c(1/4, 3/4))
ratio_rg_plots <- plot_grid(plt_20, plt_pipolin_20, nrow=1, align = "h", rel_widths=c(1/4, 3/4))

plt_pipolin_25 <- ggbetweenstats(
  data = ratio_RGs_25_pipolin %>% filter(order_num > 10),
  x = order,
  y = ratioRG,
  pairwise.display = "none"
) + ggplot2::scale_color_manual(values = color_pipolin_groups)
plt_pipolin_25

grid.arrange(ncol=1,
             plt_pipolin+ggtitle("wGRR < 0.1"), 
             plt_pipolin_20+ggtitle("wGRR < 0.2"), 
             plt_pipolin_25+ggtitle("wGRR < 0.25"))


#RG/NRG/nhNRG barplots
RG_NRGs_num <- RG_NRGs %>% group_by(MGE, Recombination_Type) %>% summarize(rec_num = n())
RG_NRGs_num$Recombination_Type <- factor(RG_NRGs_num$Recombination_Type, levels = c("NRG", "nhNRG", "RG"))

p1 <- ggplot(RG_NRGs_num, aes(x=Recombination_Type, y=rec_num, fill=Recombination_Type)) +  
            geom_bar(stat = "identity") +
            geom_text(aes(label = rec_num), color="black", hjust = "inward")+
            scale_fill_manual(values = c("#add8e6ff", "lightgray", "#cd3333ff")) +
            facet_wrap(~MGE,nrow=1, scales = "free_x") +
            coord_flip() +
            theme_classic() +
            theme(legend.position = "None") +
            labs(y = "Class", x = "Genes") 
            

p1

RG_NRGs_20_num <- RG_NRGs_20 %>% group_by(MGE, Recombination_Type) %>% summarize(rec_num = n())
RG_NRGs_20_num$Recombination_Type <- factor(RG_NRGs_20_num$Recombination_Type, levels = c("NRG", "nhNRG", "RG"))

p2 <- ggplot(RG_NRGs_20_num, aes(x=Recombination_Type, y=rec_num, fill=Recombination_Type)) +  
            geom_bar(stat = "identity") +
            geom_text(aes(label = rec_num), size = 3, color="black", hjust = "inward")+
            scale_fill_manual(values =c("#add8e6ff", "lightgray", "#cd3333ff")) +
            facet_wrap(~MGE,nrow=1, scales = "free_x") +
            coord_flip() +
            theme_classic() +
            theme(legend.position = "None") +
            labs(y = "Class", x = NULL) 
p2

ratio_rg_plots <- plot_grid(plt_20, plt_pipolin_20, nrow=1, align = "h", rel_widths=c(1/4, 3/4))
fig6_CDE <- grid.arrange(p2+ggtitle("RG assignation at wGRR < 0.2"), ratio_rg_plots, ncol=1, heights=c(2.2/7, 4.8/7))
svglite("RG_ratio_plots.svg", width = 12, height = 6)
grid.arrange(p2+ggtitle("RG assignation at wGRR < 0.2"), ratio_rg_plots, ncol=1, heights=c(2.2/7, 4.8/7))
invisible(dev.off())

RG_NRGs_30_num <- RG_NRGs_30 %>% group_by(MGE, Recombination_Type) %>% summarize(rec_num = n())
RG_NRGs_30_num$Recombination_Type <- factor(RG_NRGs_30_num$Recombination_Type, levels = c("NRG", "nhNRG", "RG"))

p3 <- ggplot(RG_NRGs_25_num, aes(x=Recombination_Type, y=rec_num, fill=Recombination_Type)) +  
  geom_bar(stat = "identity") +
  geom_text(aes(label = rec_num), color="black", hjust = "inward")+
  scale_fill_manual(values = c("darkgray", "lightgray", "red")) +
  facet_wrap(~MGE,nrow=1, scales = "free_x") +
  coord_flip() +
  theme_minimal()
p3


grid.arrange(ncol=1,
             p1+ggtitle("wGRR < 0.1"), 
             p2+ggtitle("wGRR < 0.2"), 
             p3+ggtitle("wGRR < 0.3"))


















### DEPRECATED: 
#Search for pipolin to plasmid/ciMGE/phage RGs
# load eggnog mapper and pfam cluster annotation for pipolins
pipolin_eggnog = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/pipolin_protein_clustering/reannotation_eggnog/bacteria_2022_pipolins_proteins_partial_EggNog.tsv") %>%
  select(Gene_ID=V1, COG=V7, Indiv_pfam=V21)
pipolin_gene_cluster_annotation = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/pipolin_protein_clustering/Cluster_info_simplified.txt") %>%
  select(Representative, hhblits_pfam35_tophit, Custom_category, top_hit_rec_hmm)
pipolin_gene_pfam = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/pipolin_protein_clustering/mmseqs_clustering/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db_clu.tsv",
                          header=FALSE) %>% select(Representative=V1, Gene_ID=V2) %>% left_join(pipolin_gene_cluster_annotation, by="Representative") %>%
  mutate(hhblits_pfam35_tophit_filled = if_else(is.na(hhblits_pfam35_tophit), "", hhblits_pfam35_tophit)) %>%
  select(Gene_ID, Top_pfam=hhblits_pfam35_tophit_filled, Custom_category, top_hit_rec_hmm) %>%
  left_join(pipolin_eggnog, by="Gene_ID") %>%
  mutate(COG = if_else(is.na(COG), "", COG),
         Indiv_pfam = if_else(is.na(Indiv_pfam), "", Indiv_pfam))

#Load RG assigments
BBH_08_pident_02_wgrr = fread("Pipolins_MGEs_RefSeq_RGs_20.tsv") 

# load in RGs, NRGs and add annotations 
BBH_08_pident_02_wgrr_annot = BBH_08_pident_01_wgrr %>% 
                              mutate(Gene_ID=qseqid, MGE_id=query) %>%  
                              left_join(MGE_tbl, by=c("MGE_id")) %>%
                              filter(MGE=="pipolin") %>% 
                              select(-MGE) %>%
                              mutate(MGE_id=subject)%>%
                              left_join(MGE_tbl, by=c("MGE_id")) %>%
                              left_join(pipolin_gene_pfam, by=c("Gene_ID")) 
                              
BBH_08_pident_02_wgrr_annot %>% filter(Top_pfam=="WYL")

#RG ratio to length correlation 
ggplot(ratio_RGs_pipolin %>% filter(order_num > 10), aes(x=pipolin_length, y=ratioRG, color=order)) + 
  geom_point(size=0.5, alpha=0.2) +
  scale_color_manual(values = distinctColorPalette(length(unique(ratio_RGs_pipolin$order))))


#Exploring RGs in pipolins
RGs_pipolin = RG_NRGs %>% filter(MGE=="pipolin" & Recombination_Type=="RG") %>% left_join(pipolin_info, by="MGE_id")

RGs_pipolin_alpha = RGs_pipolin %>% filter(class=="Alphaproteobacteria")



