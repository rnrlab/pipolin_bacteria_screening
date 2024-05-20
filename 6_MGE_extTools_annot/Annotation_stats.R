library("data.table")
library("tidyverse")
library("seqinr")
library("readxl")
library(ggstatsplot)
library(RColorBrewer)
library(randomcoloR)
library("gridExtra")
library("pheatmap")
library("svglite")
library(cowplot)
setDTthreads(30)

pipolin_order_color <-c("#8a2b32", "#fa7223", "#25d9c1", "#eb3f4d", "#c782ff", 
                        "#87ffb3", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",
                        "#5bad51", "#ffa947") 


#1) Read annotations file and numeber of cds per mge
MGE_tbl <- fread("../RG_assignation/MGE_info_table.tsv")
gene_annotations <- fread("Gene_annotations.tsv")

prot_info_lst = read.fasta(seqtype = "AA", as.string = T, file="../MGE_pipolins.faa")
protein_sizes_df = data.frame(protein_id = unlist(lapply(prot_info_lst,getName)), 
                              prot_length = unlist(lapply(prot_info_lst, getLength))) %>% mutate(protein_id = as.character(protein_id))
prot_number = protein_sizes_df %>% distinct() %>% 
  mutate(replicon =  str_remove(protein_id, pattern="_[0-9]{1,6}$")) %>% 
  group_by(replicon) %>% summarise(protein = n()) %>% ungroup() %>% select(MGE_id=replicon, prot_num=protein)
rm(prot_info_lst, protein_sizes_df)

gene_annotations <- gene_annotations %>% left_join(prot_number, by=c("MGE_id"))

#2) Pipolin info
pipolin_info <- fread("pipolin_summary_new.tsv") %>% select(MGE_id = Pipolin_ID, phylum, class, order, family, genus)

#Density of gene numbers
prot_number_info <- prot_number %>% left_join(MGE_tbl, by=c("MGE_id")) %>% left_join(pipolin_info, by=c("MGE_id"))

ggplot(data=prot_number_info, aes(x=prot_num, group=MGE, fill=MGE)) + geom_density(adjust=1.5, alpha=.4)

ggplot(data=prot_number_info %>% group_by(order) %>% mutate(taxa_n = n()) %>% filter(taxa_n>=10 & order!=""), 
       aes(x=prot_num, group=order, fill=order)) + geom_density(adjust=1.5, alpha=.4) + facet_wrap(~order)


#3) Stats

#3.1) PHROG
PHROG_perc <- gene_annotations %>%
  group_by(MGE) %>%
  summarise(Percentage = sum(PHROG_tag == "Viral genes", na.rm = TRUE) / n() * 100)
print(PHROG_perc)

PHROG_perc_cat <- gene_annotations %>%
  filter(PHROG_tag == "Viral genes") %>%  # Remove rows where VFDB_hit is NA
  group_by(MGE, PHROG_category) %>%
  summarise(Count = n()) %>%
  group_by(MGE) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(MGE)
PHROG_perc_cat_plot <- ggplot(PHROG_perc_cat, aes(fill=PHROG_category, y=Percentage, x=MGE)) +  
       geom_bar(position="dodge", stat="identity") +
       scale_fill_brewer(palette = "Paired") + ggtitle("A) PHROG annotation by categories")


PHROG_perc_cat_pipolins <- gene_annotations %>% filter(MGE=="pipolin") %>% left_join(pipolin_info, by=c("MGE_id")) %>%
  filter(PHROG_tag == "Viral genes") %>%  # Remove rows where VFDB_hit is NA
  group_by(order, PHROG_category) %>%
  summarise(Count = n()) %>%
  group_by(order) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(order)

phrog_cat_pipolin <- ggplot(PHROG_perc_cat_pipolins %>% filter(Count>3 & order!=""), aes(fill=PHROG_category, y=Percentage, x=order)) +  
  geom_bar(position = position_dodge(preserve = 'single'), stat="identity") +
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 

grid.arrange(PHROG_perc_cat_plot, phrog_cat_pipolin, ncol=1)

PHROG_stats_MGE_ID <-  gene_annotations %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            Vir_count = sum(PHROG_tag == "Viral genes", na.rm = TRUE),  
            Vir_perc = (Vir_count / total_count) * 100) %>%
  left_join(MGE_tbl, by=c("MGE_id")) %>% ungroup()
PHROG_plt <- ggbetweenstats(
  data = PHROG_stats_MGE_ID,
  x = MGE,
  y = Vir_perc
)
PHROG_plt

PH_means <- aggregate(Vir_perc ~  MGE, PHROG_stats_MGE_ID, mean) %>% 
  mutate(y_val = aggregate(Vir_perc ~  MGE, PHROG_stats_MGE_ID, max))
PHROG_plt_nostat <- ggbetweenstats(
  data = PHROG_stats_MGE_ID,
  x = MGE,
  y = Vir_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0),
  violin.args = list(width = 0, linewidth = 0)
) + ggplot2::scale_color_manual(values = c("#FFE66D","#0A014F","#FF6B6B","#a39022")) + ylim(0,100) + ggtitle("A - PHROGs") +
    labs(y= "Viral genes (%) ", x = NULL) + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                                                   panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(data = PH_means, aes(label = paste0(round(Vir_perc,1),"%"), y = y_val$Vir_perc + 6), color="gray40") +
  scale_x_discrete(limits=c("phage", "plasmid", "ciMGE", "pipolin")) +
  annotate(geom="text", x = 2, y = 95, label=paste0(round(PH_means$Vir_perc[4],1),"%"), color="gray40") #add plasmid value (out of range)

PHROG_plt_nostat  

PHROG_hits_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            Vir_count = sum(PHROG_tag == "Viral genes", na.rm = TRUE),  
            Vir_perc = (Vir_count / total_count) * 100) %>%
  left_join(pipolin_info, by=c("MGE_id")) %>%
  group_by(order) %>% 
  filter(n() >= 20 & order!="")

PHROG_plt_pipolin <- ggbetweenstats(
  data = PHROG_hits_pipolin,
  x = order,
  y = Vir_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  palette = "Paired" #,point.args = list(size=2, position="jitter")
) + ggplot2::scale_color_manual(values = pipolin_order_color) + ylim(0,100) + theme_classic()
PHROG_plt_pipolin

grid.arrange(PHROG_plt_nostat, PHROG_plt_pipolin,  widths = c(1/3, 2/3), ncol=2)

#3.2) CONJScan
CONJ_perc <- gene_annotations %>%
  group_by(MGE) %>%
  summarise(Percentage = sum(CONJ_tag == "Conjugation", na.rm = TRUE) / n() * 100)
print(CONJ_perc)

CONJ_stats_MGE_ID <-  gene_annotations %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            CONJ_count = sum(CONJ_tag == "Conjugation", na.rm = TRUE),  
            CONJ_perc = (CONJ_count / total_count) * 100) %>%
  left_join(MGE_tbl, by=c("MGE_id")) %>% ungroup()
CONJ_plt <- ggbetweenstats(
  data = CONJ_stats_MGE_ID,
  x = MGE,
  y = CONJ_perc 
) 
CONJ_plt

CJ_means <- aggregate(CONJ_perc ~  MGE, CONJ_stats_MGE_ID, mean) %>% 
  mutate(y_val = aggregate(CONJ_perc ~  MGE, CONJ_stats_MGE_ID, max))
CONJ_plt_nostat <- ggbetweenstats(
  data = CONJ_stats_MGE_ID,
  x = MGE,
  y = CONJ_perc,  
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0),
  violin.args = list(width = 0, linewidth = 0)
) + ggplot2::scale_color_manual(values = c("#FFE66D","#0A014F","#FF6B6B","#a39022")) + ylim(0,100) + ggtitle("B - CONJScan") +
  labs(y= "Conjugation genes (%) ", x = NULL) + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                                                       panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(data = CJ_means, aes(label = paste0(round(CONJ_perc,1),"%"), y = y_val$CONJ_perc + 6), color="gray40") +
  scale_x_discrete(limits=c("phage", "plasmid", "ciMGE", "pipolin"))


CONJ_plt_nostat 
grid.arrange(PHROG_plt_nostat, CONJ_plt_nostat, ncol=2)

# Calculate percentage of elements with 'yes' and 'no' for each MGE class
CONJ_presence_MGE_ID <- gene_annotations %>%
  group_by(MGE_id) %>%
  mutate(conjugation = ifelse("Conjugation" %in% CONJ_tag, "yes", "no")) %>%
  ungroup() %>% select(MGE, MGE_id, conjugation) %>% distinct()
CONJ_presence_MGE_ID_perc <- CONJ_presence_MGE_ID %>%
  group_by(MGE, conjugation) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
ggplot(CONJ_presence_MGE_ID_perc, aes(x = MGE, y = percentage, fill = conjugation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Percentage of Elements with Conjugation genes",
       x = "Class", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Categories and MGEs
Conj_perc_cat <- gene_annotations %>%
  filter(CONJ_tag == "Conjugation") %>%  # Remove rows where VFDB_hit is NA
  group_by(MGE, CONJ_hmm) %>%
  summarise(Count = n()) %>%
  group_by(MGE) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(MGE)
Conj_perc_cat_plot <- ggplot(Conj_perc_cat%>%filter(Count>3), aes(y=Percentage, x=CONJ_hmm)) +  
  geom_bar(position="stack", stat="identity") + #dodge
  facet_grid(MGE~., scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7)) +
  ggtitle("B) CONJScan annotation by profiles")
Conj_perc_cat_plot
#Now for pipolins
Conj_perc_cat_pipolins <- gene_annotations %>% filter(MGE=="pipolin") %>% left_join(pipolin_info, by=c("MGE_id")) %>%
  filter(CONJ_tag == "Conjugation") %>%  # Remove rows where VFDB_hit is NA
  group_by(order, CONJ_hmm) %>%
  summarise(Count = n()) %>%
  group_by(order) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(order)

Conj_perc_cat_pip_plot <- ggplot(Conj_perc_cat_pipolins %>% filter(Count > 3), aes(fill=CONJ_hmm, y=Count, x=order)) +  
  geom_bar(position = position_dodge(preserve = 'single'), stat="identity") +
  scale_fill_manual(values = distinctColorPalette(length(unique(Conj_perc_cat$CONJ_hmm)))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

svglite("Fig_S_PHROG_CONJSCan_categories.svg", width = 15, height = 20)
grid.arrange(PHROG_perc_cat_plot, phrog_cat_pipolin, Conj_perc_cat_plot, Conj_perc_cat_pip_plot, 
             ncol=1, heights=c(2/10,2/10,4/10,2/10))
invisible(dev.off())


CONJ_hits_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            CONJ_count = sum(CONJ_tag == "Conjugation", na.rm = TRUE),  
            CONJ_perc = (CONJ_count / total_count) * 100) %>% ungroup() %>%
  left_join(pipolin_info, by=c("MGE_id")) %>%
  group_by(order) %>% mutate(taxa_n = n()) %>%
  filter(taxa_n>=10 & order!="")

CONJ_plt_pipolin <- ggbetweenstats(
  data = CONJ_hits_pipolin,
  x = order,
  y = CONJ_perc,
  palette = "Paired",
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
)+ ggplot2::scale_color_manual(values = distinctColorPalette(length(unique(CONJ_hits_pipolin$order)))) + ylim(0,100)
CONJ_plt_pipolin

grid.arrange(CONJ_plt_nostat, CONJ_plt_pipolin,  widths = c(1/3, 2/3), ncol=2)

# Calculate percentage of elements with 'yes' and 'no' for each pipolin class
CONJ_presence_MGE_ID_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>% 
  group_by(MGE_id) %>%
  mutate(conjugation = ifelse("Conjugation" %in% CONJ_tag, "yes", "no")) %>%
  ungroup() %>% left_join(pipolin_info, by=c("MGE_id")) %>% select(order, MGE_id, conjugation) %>% distinct() %>% 
  group_by(order) %>% filter(n()>=10 & order != "") 

CONJ_presence_MGE_ID_perc_pipolin <- CONJ_presence_MGE_ID_pipolin %>%
  group_by(order, conjugation) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
ggplot(CONJ_presence_MGE_ID_perc_pipolin, aes(x = order, y = percentage, fill = conjugation)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Percentage of Elements with Conjugation genes", x = "Class", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#3.3) VF
VF_perc <- gene_annotations %>%
           group_by(MGE) %>%
          summarise(Percentage_VF = sum(VFDB_hit == "Virulence Factor", na.rm = TRUE) / n() * 100)
print(VF_perc)


VFDB_stats_MGE_ID <-  gene_annotations %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            VF_count = sum(VFDB_hit == "Virulence Factor", na.rm = TRUE),  
            VF_perc = (VF_count / total_count) * 100) %>%
  left_join(MGE_tbl, by=c("MGE_id")) %>% ungroup()
VF_plt <- ggbetweenstats(
  data = VFDB_stats_MGE_ID,
  x = MGE,
  y = VF_perc
)
VF_plt 

VF_means <- aggregate(VF_perc ~  MGE, VFDB_stats_MGE_ID, mean) %>% 
  mutate(y_val = aggregate(VF_perc ~  MGE, VFDB_stats_MGE_ID, max))
VF_plt_nostat <- ggbetweenstats(
  data = VFDB_stats_MGE_ID,
  x = MGE,
  y = VF_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  violin.args = list(width = 0, linewidth = 0),
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0)
)  + ggplot2::scale_color_manual(values = c("#FFE66D","#0A014F","#FF6B6B","#a39022")) + ylim(0,100) + ggtitle("C - VFDB") +
  labs(y= "Virulence genes (%) ", x = NULL) + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                                                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(data = VF_means, aes(label = paste0(round(VF_perc,1),"%"), y = y_val$VF_perc + 6), color="gray40") +
  scale_x_discrete(limits=c("phage", "plasmid", "ciMGE", "pipolin")) +
  annotate(geom="text", x = 2, y = 96, label=paste0(round(VF_means$VF_perc[4],1),"%"), color="gray40") #add plasmid value (out of range)

VF_plt_nostat 

grid.arrange(PHROG_plt_nostat, CONJ_plt_nostat, VF_plt_nostat, nrow=1)

VFDB_hits_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>% group_by(MGE_id) %>%
                                          summarise(total_count = n(), 
                                                    VF_count = sum(VFDB_hit == "Virulence Factor", na.rm = TRUE),  
                                                    VF_perc = (VF_count / total_count) * 100) %>% ungroup() %>%
                                                    left_join(pipolin_info, by=c("MGE_id")) %>%
                                                    group_by(order) %>% mutate(taxa_n = n()) %>%
                                                    filter(taxa_n>=10 & order!="")
  
VF_plt_pipolin <- ggbetweenstats(
  data = VFDB_hits_pipolin,
  x = order,
  y = VF_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
)+ ggplot2::scale_color_manual(values = distinctColorPalette(length(unique(VFDB_hits_pipolin$order)))) + ylim(0,100)
VF_plt_pipolin

grid.arrange(VF_plt_nostat, VF_plt_pipolin,  widths = c(1/3, 2/3), ncol=2)

VF_perc_cat <- gene_annotations %>% filter(MGE=="pipolin") %>% left_join(pipolin_info, by=c("MGE_id")) %>%
  filter(VFDB_hit == "Virulence Factor") %>%  # Remove rows where VFDB_hit is NA
  group_by(genus, VF_class) %>%
  summarise(Count = n()) %>%
  group_by(genus) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(genus)

ggplot(VF_perc_cat, aes(fill=VF_class, y=Count, x=genus)) +  
  geom_bar(position = position_dodge(preserve = 'single'), stat="identity") +
  scale_fill_manual(values = distinctColorPalette(length(unique(VF_perc_cat$VF_class)))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



VFDB_hits_pipolin_tmp <- gene_annotations %>% filter(MGE=="pipolin" & VF_class == "Immune modulation (VFC0258)") %>%
                                              left_join(pipolin_info, by=c("MGE_id")) %>%
                                              filter(genus=="Escherichia")

#3.4) AMRFinder
AMR_perc <- gene_annotations %>%
  group_by(MGE) %>%
  summarise(Percentage = sum(AMR_tag == "AMR", na.rm = TRUE) / n() * 100)
print(AMR_perc)

AMR_stats_MGE_ID <-  gene_annotations %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            AMR_count = sum(AMR_tag == "AMR", na.rm = TRUE),  
            AMR_perc = (AMR_count / total_count) * 100) %>%
  left_join(MGE_tbl, by=c("MGE_id")) %>% ungroup()

amr_means <- aggregate(AMR_perc ~  MGE, AMR_stats_MGE_ID, mean) %>% 
  mutate(y_val = aggregate(AMR_perc ~  MGE, AMR_stats_MGE_ID, max))
AMR_plt <- ggbetweenstats(
  data = AMR_stats_MGE_ID,
  x = MGE,
  y = AMR_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  violin.args = list(width = 0, linewidth = 0),
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0)
) + ggplot2::scale_color_manual(values = c("#FFE66D","#0A014F","#FF6B6B","#a39022")) + ylim(0,100) + ggtitle("D - AMRFinderPlus") +
  labs(y= "AMR genes (%) ", x = NULL) + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                                               panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(data = amr_means, aes(label = paste0(round(AMR_perc,1),"%"), y = y_val$AMR_perc + 6), color="gray40") + 
  scale_x_discrete(limits=c("phage", "plasmid", "ciMGE", "pipolin")) 

AMR_plt

plot4 <- grid.arrange(PHROG_plt_nostat, CONJ_plt_nostat, VF_plt_nostat, AMR_plt, nrow=1)
ggsave(plot=plot4, file="Fig_annotation_4plots.svg", width = 12,  height = 4)



AMR_hits_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>% group_by(MGE_id) %>%
  summarise(total_count = n(), 
            AMR_count = sum(AMR_tag == "AMR", na.rm = TRUE),  
            AMR_perc = (AMR_count / total_count) * 100) %>% ungroup() %>%
  left_join(pipolin_info, by=c("MGE_id")) %>%
  group_by(order) %>% mutate(taxa_n = n()) %>%
  filter(taxa_n>=10 & order!="")


AMR_plt_pipolin <- ggbetweenstats(
  data = AMR_hits_pipolin,
  x = order,
  y = AMR_perc,  
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE) + ggplot2::scale_color_manual(values = distinctColorPalette(length(unique(AMR_hits_pipolin$order)))) + ylim(0,100)
AMR_plt_pipolin

grid.arrange(AMR_plt, AMR_plt_pipolin,  widths = c(1/3, 2/3), ncol=2)

#3.5) YR and SR
Rec_perc <- gene_annotations %>%
  group_by(MGE) %>%
  summarise(Percentage = sum(Rec_tag == "Recombinase (except IS)", na.rm = TRUE) / n() * 100)
print(Rec_perc)

#3.X) PADLOC
PADLOC_perc <- gene_annotations %>%
  group_by(MGE) %>%
  summarise(Percentage = sum(PADLOC_tag == "Defense", na.rm = TRUE) / n() * 100)
print(PADLOC_perc)

PADLOC_stats_MGE_ID <-  gene_annotations %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            Def_count = sum(PADLOC_tag == "Defense", na.rm = TRUE),  # Count of "A" values
            Def_perc = (Def_count / total_count) * 100) %>%
            left_join(MGE_tbl, by=c("MGE_id")) %>% left_join(prot_number, by=c("MGE_id")) %>% ungroup()

PADLOC_plt <- ggbetweenstats(
  data = PADLOC_stats_MGE_ID,
  x = MGE,
  y = Def_perc
)

PADLOC_plt

def_means <- aggregate(Def_perc ~  MGE, PADLOC_stats_MGE_ID, mean) %>% 
  mutate(y_val = aggregate(Def_perc ~  MGE, PADLOC_stats_MGE_ID, max))
PADLOC_plt_nostats <- ggbetweenstats(
  data = PADLOC_stats_MGE_ID,
  x = MGE,
  y = Def_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  violin.args = list(width = 0, linewidth = 0),
  results.subtitle = FALSE,
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0)
) + ggplot2::scale_color_manual(values = c("#FFE66D","#0A014F","#FF6B6B","#a39022")) + ylim(0,100) + ggtitle("E - PADLOC (All MGE)") +
  labs(y= "Defense genes (%) ", x = NULL) + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                                                   panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_text(data = def_means, aes(label = paste0(round(Def_perc,1),"%"), y = y_val$Def_perc + 6), color="gray40") +
  scale_x_discrete(limits=c("phage", "plasmid", "ciMGE", "pipolin")) 
PADLOC_plt_nostats

#Bias of low amount of genes? 
p1 <- ggplot(PADLOC_stats_MGE_ID, aes(x=prot_num, y=Def_perc, color=MGE))+ geom_point(size=0.2,alpha=0.2) #percentage
p2 <- ggplot(PADLOC_stats_MGE_ID, aes(x=prot_num, y=Def_count, color=MGE))+ geom_point(size=0.2,alpha=0.2) #abs num
grid.arrange(p1, p2, ncol=2)

#Now with facet
p1 <- ggplot(PADLOC_stats_MGE_ID, aes(x=prot_num, y=Def_perc, color=MGE))+ geom_point(size=0.2,alpha=0.2) + facet_wrap(~MGE, nrow=1)
p2 <- ggplot(PADLOC_stats_MGE_ID, aes(x=prot_num, y=Def_count, color=MGE))+ geom_point(size=0.2,alpha=0.2) + facet_wrap(~MGE, nrow=1)
grid.arrange(p1, p2, nrow=2)

#Only similar size MGE
p1 <- ggplot(PADLOC_stats_MGE_ID, aes(x=prot_num, y=Def_perc, color=MGE))+ geom_point(size=1.5,alpha=0.2) + xlim(0,100)
p2 <- ggplot(PADLOC_stats_MGE_ID, aes(x=prot_num, y=Def_count, color=MGE))+ geom_point(size=1.5,alpha=0.2) + xlim(0,100)
grid.arrange(p1, p2, ncol=2)

#2D density plot
PADLOC_stats_MGE_ID_p100 <- PADLOC_stats_MGE_ID %>% filter(prot_num<=100 & Def_count>0)
p1 <- ggplot(PADLOC_stats_MGE_ID_p100, aes(prot_num, Def_perc))  +
  stat_density2d(geom="tile", aes(fill = after_stat(density)), contour = FALSE) +
  facet_wrap(~MGE)+  geom_point(color = "white", size=0.1, alpha=0.02)
p2 <- ggplot(PADLOC_stats_MGE_ID_p100, aes(prot_num, Def_count))  +
  stat_density2d(geom="tile", aes(fill = after_stat(density)), contour = FALSE) +
  facet_wrap(~MGE)+  geom_point(color = "white", size=0.1, alpha=0.02)
grid.arrange(p1, p2, ncol=2)

#filter elements <= 50 prot
PADLOC_stats_MGE_ID_50 <-  gene_annotations %>% group_by(MGE_id) %>%
  summarise(total_count = n(),  # Total count of rows for each MGE_id
            Def_count = sum(PADLOC_tag == "Defense", na.rm = TRUE),  
            Def_perc = (Def_count / total_count) * 100) %>%
  left_join(MGE_tbl, by=c("MGE_id"))%>% left_join(prot_number, by=c("MGE_id")) %>% filter(prot_num<=50) %>% ungroup()

PADLOC_plt_50 <- ggbetweenstats(
  data = PADLOC_stats_MGE_ID_50,
  x = MGE,
  y = Def_perc
)
PADLOC_plt_50
grid.arrange(PADLOC_plt, PADLOC_plt_50, ncol=2)

# Calculate percentage of elements with 'yes' and 'no' for each MGE class
PADLOC_presence_MGE_ID <- gene_annotations %>%
  group_by(MGE_id) %>%
  mutate(def = ifelse("Defense" %in% PADLOC_tag, "yes", "no")) %>%
  ungroup() %>% select(MGE, MGE_id, def) %>% distinct()
PADLOC_presence_MGE_ID_perc <- PADLOC_presence_MGE_ID %>%
  group_by(MGE, def) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
PADLOC_presence_MGE_p1 <- ggplot(PADLOC_presence_MGE_ID_perc, aes(x = MGE, y = percentage, fill = def)) +
                          geom_bar(stat = "identity", position = "fill") +
                          scale_y_continuous(labels = scales::percent) +
                          labs(title = "Percentage of Elements with defense genes",
                               x = "Class", y = "Percentage") +
                          theme_minimal() +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1))


PADLOC_sys_types <- gene_annotations %>%
                      filter(PADLOC_system != "") %>%
                      group_by(MGE, PADLOC_system) %>%
                      summarise(sys_genes =n())
ggplot(PADLOC_sys_types, aes(fill=MGE, y=sys_genes, x=PADLOC_system)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~MGE, ncol=1)


#In pipolin groups
PADLOC_hits_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>%   
                       group_by(MGE_id) %>%
                       summarise(total_count = n(), 
                                  Def_count = sum(PADLOC_tag == "Defense", na.rm = TRUE),  
                                  Def_perc = (Def_count / total_count) * 100) %>% ungroup() %>%
                       left_join(pipolin_info, by=c("MGE_id")) %>%
                       group_by(order) %>% 
                       filter(n() >= 20 & order!="") 
                       

### plot aes ##########################################################################################
pipolin_order_color <-c("#8a2b32", "#fa7223", "#25d9c1", "#eb3f4d", "#c782ff", 
                        "#87ffb3", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",
                        "#5bad51", "#ffa947") 

def_pipolin_means <- aggregate(Def_perc ~  order, PADLOC_hits_pipolin, mean) %>% 
                     mutate(y_val = aggregate(Def_perc ~  order, PADLOC_hits_pipolin, max))

order_n <- PADLOC_hits_pipolin %>% group_by(order) %>% summarise(n = n())
x_label_names <- paste0(order_n$order, "\n(n = ",order_n$n,")")
names(x_label_names) <- order_n$order
#######################################################################################

PADLOC_plt_pipolin <- ggbetweenstats(
  data = PADLOC_hits_pipolin,
  x = order,
  y = Def_perc,
  pairwise.comparisons = FALSE,
  pairwise.display = "none",
  results.subtitle = FALSE,
  centrality.point.args = list(size = 1, color = "red"),
  centrality.label.args = list(alpha = 0),
  violin.args = list(width = 0, linewidth = 0)
  ) +
  ggplot2::scale_color_manual(values = pipolin_order_color) + ylim(0,100) + ggtitle("F - PADLOC (Pipolins by host)") +
  labs(y= "Defense genes (%) ", x = "Host order") + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                                                          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                                          axis.text.x=element_text(size=7)) +
  geom_text(data = def_pipolin_means, aes(label = paste0(round(Def_perc,1),"%"), y = y_val$Def_perc + 6), color="gray40") +
  scale_x_discrete(limits=c("Enterobacterales", "Vibrionales", "Aeromonadales", "Alteromonadales",
                            "Hyphomicrobiales", "Rhodobacterales", 
                            "Micrococcales", "Mycobacteriales",
                            "Lachnospirales", "Eubacteriales",
                            "Lactobacillales", "Bacillales"),
                   labels=x_label_names) + 
  theme(axis.text.x = element_text(angle = 45, size=7.5, vjust=0.5))
  

PADLOC_plt_pipolin 



PADLOC_sys_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>%
                                           left_join(pipolin_info, by=c("MGE_id")) %>%
                                           filter(PADLOC_system != "") %>%
                                           group_by(order) %>% 
                                           filter(order %in% c("Enterobacterales", "Vibrionales", "Aeromonadales", "Alteromonadales",
                                                               "Hyphomicrobiales", "Rhodobacterales", 
                                                               "Micrococcales", "Mycobacteriales",
                                                               "Lachnospirales", "Eubacteriales",
                                                               "Lactobacillales")) %>%
                                           group_by(order, PADLOC_system) %>%
                                           summarise(sys_genes=n(),
                                                     sys_pipolins = length(unique(MGE_id))) %>%
                                           mutate( total_sys = sum(sys_pipolins),
                                                     p_system = sys_pipolins/total_sys,
                                                     log2_p_system = log2(p_system),
                                                     p_log2_p_system = p_system*log2_p_system) 
PADLOC_sys_pipolin_entropy <- PADLOC_sys_pipolin %>%  group_by(order) %>%
                                                     summarise(sys_entropy = -sum(p_log2_p_system)) 
                                                      

#Note: We assume there are no >1 system copy per pipolin 

pipolin_order_color <-c("#8a2b32", "#fa7223", "#25d9c1", "#eb3f4d", "#c782ff", 
                        "#87ffb3", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",
                        "#5bad51", "#ffa947") 

#Entropy plot 
entropy_plot <- ggplot(data=PADLOC_sys_pipolin_entropy, aes(x=reorder(order, sys_entropy), y=sys_entropy, fill=order)) +
  geom_bar(stat="identity", width = 0.7) + 
  scale_fill_manual(values=pipolin_order_color[-3]) +
  coord_flip() + theme_classic() + theme(legend.position = "none", axis.text.y=element_blank()) +
  ggtitle("H - Entropy") + labs(x = NULL, y = "Bits")

entropy_plot

#heatpmap systems
heatmap_order <- PADLOC_sys_pipolin_entropy[order(PADLOC_sys_pipolin_entropy$sys_entropy),]$order
         
def_heatmap <- ggplot(PADLOC_sys_pipolin %>% filter(sys_pipolins >= 1), aes(x=order, y=PADLOC_system)) +
                geom_tile(aes(fill = sys_pipolins)) +
                #geom_text(aes(label = round(sys_pipolins, 0)), size = 2, color="white") + 
                scale_fill_gradient(low = "#ff9191", high = "darkred", trans = 'log') +
                theme_bw()+
                theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7.5),
                      axis.text.y = element_text(size = 8, face="italic")) + 
                coord_flip() + ggtitle("G - Defence system frequency among pipolins") +
                labs(x = "Order", y="System") + scale_x_discrete(limits=heatmap_order)# +
                theme(legend.position = "bottom")
  
def_heatmap
#grid.arrange(def_heatmap, entropy_plot, nrow=1, widths=c(2/3,1/3))
plot_grid(def_heatmap, entropy_plot, nrow=1, align = "h", rel_widths=c(6/7,1/7))

plot4 <- grid.arrange(PHROG_plt_nostat, CONJ_plt_nostat, VF_plt_nostat, AMR_plt, nrow=1)
padloc_grid <- plot_grid(PADLOC_plt_nostats, PADLOC_plt_pipolin, nrow=1, align = "h", rel_widths = c(1/4, 3/4))
padloc_variety <- plot_grid(def_heatmap, entropy_plot, nrow=1, align = "h", rel_widths=c(6/7,1/7))

svglite("Fig_5_annotation_stats.svg", width = 14, height = 15)
grid.arrange(plot4, padloc_grid, padloc_variety, nrow=3, heights = c(2/8, 2.5/8, 3.5/7))
invisible(dev.off())


#barplots systems
subplot_def <- function(order_name) {
  ggplot(PADLOC_sys_pipolin  %>% filter(order==order_name), aes(fill=order, y=sys_pipolins, x=PADLOC_system)) + 
    geom_bar(position="stack", stat="identity") +
    facet_wrap(~order, scales = "free", nrow=1) +
    ggtitle("System diveristy in pipolins") +
    labs(y= "Pipolins") + 
    scale_fill_manual(values = c("#ffa947","#eb3f4d")) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), legend.position = "none") 
}

def_row_1 <- ggplot(PADLOC_sys_pipolin  %>% filter(order=="Vibrionales" | order=="Enterobacterales" ), aes(fill=order, y=log2(sys_pipolins), x=PADLOC_system)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~order, scales = "free", space="free") +
  ggtitle("System diveristy in pipolins") +
  labs(y= "log2(Pipolins)", x=NULL) + 
  scale_fill_manual(values = c("#eb3f4d","#ffa947")) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), legend.position = "none") 
def_row_1

def_row_2 <- ggplot(PADLOC_sys_pipolin  %>% filter(order=="Hyphomicrobiales" | order=="Aeromonadales" | order=="Alteromonadales"), aes(fill=order, y=sys_pipolins, x=PADLOC_system)) + 
                    geom_bar(position="stack", stat="identity") +
                    facet_grid(~order, scales = "free", space="free")+
                    labs(y= "Pipolins", x=NULL) + 
                    scale_fill_manual(values = c("#8a2b32", "#fa7223","#87ffb3")) + 
                    theme_classic() + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), legend.position = "none") 
def_row_2

def_row_3 <- ggplot(PADLOC_sys_pipolin  %>% filter(order!="Vibrionales" & order!="Enterobacterales" & order!="Aeromonadales"  & order!="Hyphomicrobiales" & order!="Alteromonadales"), aes(fill=order, y=sys_pipolins, x=PADLOC_system)) + 
                    geom_bar(position="stack", stat="identity") +
                    facet_grid(~order, scales = "free", space="free") +
                    labs(y= "Pipolins", x = "Host order") + 
                    scale_fill_manual(values = c("#c782ff", "#4781e6","#8e2abf","#e32bac","#ff94f1",
                                                 "#5bad51")) + 
                    theme_classic() + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8), legend.position = "none") 
def_row_3

plot4 <- grid.arrange(PHROG_plt_nostat, CONJ_plt_nostat, VF_plt_nostat, AMR_plt, nrow=1)
padloc_grid <- grid.arrange(PADLOC_plt_nostats, PADLOC_plt_pipolin, nrow=1, widths = c(1/4, 3/4))

padloc_pipolin <- grid.arrange(def_row_1, def_row_2, def_row_3, nrow=3)

annot_1 <- grid.arrange(plot4, padloc_grid, nrow=2, heights = c(1/2,1/2))
ggsave(plot=annot_1, file="Fig_annotation_plots_1.pdf", width = 10,  height = 20)
ggsave(plot=padloc_pipolin, file="Fig_annotation_plots_2.pdf", width = 10,  height = 10)


grid.arrange(plot4, padloc_grid, def_row_1, def_row_2, def_row_3, nrow=5)

# Calculate percentage of elements with 'yes' and 'no' for each pipolin class
PADLOC_presence_MGE_ID_pipolin <- gene_annotations %>% filter(MGE=="pipolin") %>% 
  group_by(MGE_id) %>%
  mutate(def = ifelse("Defense" %in% PADLOC_tag, "yes", "no")) %>%
  ungroup() %>% left_join(pipolin_info, by=c("MGE_id")) %>% select(order, MGE_id, def) %>% distinct() %>% 
  group_by(order) %>% filter(n()>=5 & order != "") 
PADLOC_presence_MGE_ID_perc_pipolin <- PADLOC_presence_MGE_ID_pipolin %>%
  group_by(order, def) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
PADLOC_presence_MGE_p2 <- ggplot(PADLOC_presence_MGE_ID_perc_pipolin, aes(x = order, y = percentage, fill = def)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Percentage of pipolins with defense genes", x = "Class", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
grid.arrange(PADLOC_presence_MGE_p1, PADLOC_presence_MGE_p2, ncol=2, widths = c(1/4, 3/4))





