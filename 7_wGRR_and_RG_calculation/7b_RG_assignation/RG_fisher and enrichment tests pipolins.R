# libraries
library("data.table")
library("RColorBrewer")
library("tidyverse")
library("rstatix")
library("gridExtra")
library(cowplot)
library(svglite)

# data sets

## load in pipolin, ciMGEs, plasmid and phage tbl (V_ not needed)
MGE_tbl = fread("MGE_info_table.tsv")

# load in pipolin info
pipolin_info = fread("../Gene_annotation/pipolin_summary_new.tsv")
colnames(pipolin_info)[2] <- "MGE_id"

# load eggnog mapper and pfam cluster annotation for pipolins
pipolin_eggnog = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/pipolin_protein_clustering/reannotation_eggnog/bacteria_2022_pipolins_proteins_partial_EggNog.tsv") %>%
                      select(Gene_ID=V1, COG=V7, Indiv_pfam=V21)
pipolin_gene_cluster_annotation = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/pipolin_protein_clustering/Cluster_info_simplified.txt") %>%
  select(Representative, hhblits_pfam35_tophit, Custom_category, top_hit_rec_hmm)
pipolin_gene_pfam = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/pipolin_protein_clustering/Clustered_db_clu_fixedNames.tsv",
                          header=FALSE) %>% select(Representative=V1, Gene_ID=V2) %>% left_join(pipolin_gene_cluster_annotation, by="Representative") %>%
                    mutate(hhblits_pfam35_tophit_filled = if_else(is.na(hhblits_pfam35_tophit), "", hhblits_pfam35_tophit)) %>%
                    select(Gene_ID, Top_pfam=hhblits_pfam35_tophit_filled, Custom_category, top_hit_rec_hmm) %>%
                    left_join(pipolin_eggnog, by="Gene_ID") %>%
                    mutate(COG = if_else(is.na(COG), "", COG),
                           Indiv_pfam = if_else(is.na(Indiv_pfam), "", Indiv_pfam))

# load in RGs, NRGs and add annotations 
RG_NRGs_list = fread("RG_NRG_list_20.tsv") %>% mutate(MGE_id = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>% left_join(MGE_tbl, by="MGE_id")

Gene_annotations = fread("../Gene_annotation/Gene_annotations.tsv") %>% select(-MGE_id, -MGE) 


RG_NRGs = RG_NRGs_list %>% left_join(Gene_annotations, by=c("Gene_ID")) %>%
                           filter(MGE=="pipolin") %>%
                           left_join(pipolin_gene_pfam, by=c("Gene_ID")) %>% 
                           left_join(pipolin_info, by=c("MGE_id")) %>%
                           group_by(order) %>% 
                           mutate(num_taxa = n_distinct(MGE_id)) %>%
                           filter(num_taxa >= 20) %>% ungroup() %>% filter(order!="") %>%
                           mutate(Top_pfam_fixed = if_else(is.na(Top_pfam), "", Top_pfam), #fix NAs
                                 COG_fixed = if_else(is.na(COG), "", COG),
                                 Indiv_pfam_fixed = if_else(is.na(Indiv_pfam), "", Indiv_pfam),
                                 Custom_category = if_else(is.na(Custom_category), "", Custom_category)) %>%
                           mutate(COG_fixed = if_else(COG_fixed=="-", "", COG_fixed),
                                  Indiv_pfam_fixed = if_else(Indiv_pfam_fixed=="-", "", Indiv_pfam_fixed)) %>%
                           select(Gene_ID, MGE_id, taxa=order, Recombination_Type, PADLOC_system, 
                                  Rec_class, AMR_class, VF_class, CONJ_hmm, PHROG_category, 
                                  Top_pfam=Top_pfam_fixed, COG=COG_fixed, Indiv_pfam=Indiv_pfam_fixed, Custom_category, top_hit_rec_hmm) 
                            
#these rows are for testing
RG_pairs_20 <- fread("../RG_assignation/Pipolins_MGEs_RefSeq_RGs_20.tsv") %>% left_join(MGE_tbl, by=c("subject"="MGE_id"))
RG_NRGs_test = RG_NRGs %>% filter(taxa=="Sphingomonadales")
rows_with_na <- RG_NRGs %>% select(Gene_ID, Top_pfam, COG, Indiv_pfam) %>% filter(rowSums(is.na(.)) > 0)

rg_query <- RG_NRGs %>% filter(taxa=="Enterobacterales" & PADLOC_system=="mza")
view(RG_pairs_20 %>%  left_join(rg_query, by=c("qseqid"="Gene_ID")) %>% filter(taxa=="Enterobacterales" & PADLOC_system=="mza"))

view(RG_NRGs %>% filter(MGE_id=="G_371040_0v0"))

#stacked barplots RG/NRG/nhNRG
RG_NRGs_test_for_plot = data.frame(NULL); for (c in (RG_NRGs %>% names)[c(5,10,8,9,11)]) { #10,5,9,8
  tmp_df <- RG_NRGs %>% select(taxa, Recombination_Type, subcategory=c) %>%  
                        filter(subcategory != "" & subcategory != "-") %>%
                        group_by(taxa, Recombination_Type, subcategory) %>%
                        summarise(Genes=n(), log_genes=log2(n())) %>% ungroup() %>%
                        group_by(taxa) %>%
                        mutate(c_n = sum(Genes)) %>%
                        group_by(taxa,  subcategory) %>%
                        filter(sum(Genes)/c_n >= 0.005 & sum(Genes) >= 3) %>% #0.005, 10
                        mutate(ratio = sum(Genes)/c_n,
                               sum_genes = sum(Genes),
                               Genes = if_else(Recombination_Type != "RG", -Genes, Genes),
                               log_genes = if_else(Recombination_Type != "RG", -log_genes, log_genes),
                               category=c)
  RG_NRGs_test_for_plot = bind_rows(RG_NRGs_test_for_plot, tmp_df)
}

### ALL ###
ggplot(RG_NRGs_test_for_plot, aes(fill=Recombination_Type, y=log_genes, x=reorder(subcategory, -Genes))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
  coord_flip() +
  theme_bw(base_size = 12) +
  facet_grid(~taxa, scales = "free",   drop=TRUE) + ggtitle("wGRR < 0.2")

### Enterobacterales ###
E1 <- ggplot(RG_NRGs_test_for_plot %>%  filter(taxa=="Enterobacterales" & category != "Top_pfam")  , 
            aes(fill=Recombination_Type, y=Genes, x=reorder(subcategory, -Genes))) + 
            geom_bar(position="stack", stat="identity") +
            scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
            coord_flip() +
            theme_bw(base_size=10) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
            facet_grid(category~., scales = "free", space="free")  + ggtitle("RGs in Enterobacterales") + ylab(NULL) + xlab(NULL)
E1
E2 <- ggplot(RG_NRGs_test_for_plot %>%  filter(taxa=="Enterobacterales" & category == "Top_pfam")  , 
             aes(fill=Recombination_Type, y=Genes, x=reorder(subcategory, -Genes))) + 
             geom_bar(position="stack", stat="identity") +
             scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
             coord_flip() +
             theme_bw(base_size=10) +
             theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
             facet_grid(category~., scales = "free", space="free")  + ggtitle("") + xlab(NULL)
E2

E_join <- plot_grid(E1, E2, ncol=1, align = "v", rel_heights = c(2/5,3/5))
E_final <- grid.arrange(E_join, p_fixed, ncol=2,
             widths = c(1.5/5,3.5/5)) #p from gggenomes in another script (needs R studio)
svglite("RGs_enterobacteria.svg", width = 15, height = 7)
plot(E_final)
invisible(dev.off())


#alt: facet_grid(taxa ~ ., scales = "free", space = "free") + ggtitle("wGRR < 0.2")
#alt: facet_grid(~taxa, scales = "free",   drop=TRUE) + ggtitle("wGRR < 0.2")
#alt: facet_wrap(~taxa,  scales = "free",  drop=TRUE, ncol = 4) + ggtitle("wGRR < 0.2")
#alt: grid.arrange for each category


### REST OF ORDERS ###
Rdef <- ggplot(RG_NRGs_test_for_plot %>%  filter((taxa!="Vibrionales"&taxa!="Enterobacterales"&taxa!="Eubacteriales"&taxa!="Alteromonadales"&taxa!="Rhodobacterales") & category == "PADLOC_system")  , 
               aes(fill=Recombination_Type, y=Genes, x=reorder(subcategory, -Genes))) +  #Alt:log_genes
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
  coord_flip() +
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
  facet_grid(~taxa, scales = "free")  + labs(x=NULL, y = NULL)
Rdef

Rdef_vibrio <- ggplot(RG_NRGs_test_for_plot %>%  filter((taxa=="Vibrionales") & category == "PADLOC_system")  , 
                      aes(fill=Recombination_Type, y=Genes, x=reorder(subcategory, -Genes))) +  #Alt:log_genes
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
  coord_flip() +
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
  facet_grid(~taxa, scales = "free")  + labs(x=NULL, y = NULL)
Rdef_vibrio

Rdef_join <- grid.arrange(Rdef_vibrio, Rdef, nrow=1, widths=c(3/7,4/7))

# tmp categ other
phrog_conj_rest <- ggplot(RG_NRGs_test_for_plot %>%  filter(taxa!="Enterobacterales" & (category == "PHROG_category" | category =="CONJ_hmm"))  , 
                          aes(fill=Recombination_Type, y=Genes, x=reorder(subcategory, -Genes))) +  #Alt:log_genes
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
  coord_flip() +
  theme_bw(base_size=9.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none") + 
  facet_grid(category~factor(taxa, levels = c("Vibrionales", "Aeromonadales", "Alteromonadales",
                                              "Hyphomicrobiales", "Rhodobacterales", 
                                              "Micrococcales", "Mycobacteriales",
                                              "Lachnospirales", "Eubacteriales",
                                              "Lactobacillales", "Bacillales")), scales = "free") +
  labs(x = NULL)
phrog_conj_rest
svglite("RGs_other_nopfam2.svg", width = 12, height = 7)
grid.arrange(Rdef_join, phrog_conj_rest, ncol=1)
invisible(dev.off())


### Pfam RG plots (supplementary figure)
Rpfam <- ggplot(RG_NRGs_test_for_plot %>%  filter(category == "Top_pfam" & sum_genes > 10)  , 
             aes(fill=Recombination_Type, y=Genes, x=reorder(subcategory, -Genes))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("darkgray","lightblue","brown3")) +
  coord_flip() +
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", axis.text.x = element_text(angle = 270, hjust=0)) + 
  facet_wrap(~taxa, scales = "free", ncol=2)  + ggtitle("Top cluster PFAM RG/NRG distribution")
Rpfam


svglite("Fig_S_RGs_pfam.svg", width = 12, height = 24)
plot(Rpfam)
invisible(dev.off())




### Significance tests ###

Fischer_complete = data.frame(NULL); for (c in (RG_NRGs %>% names)[c(5,8,9,10,11)]) {  

### 1) Count NRG and RG for each taxa and function category  
count_tbl = RG_NRGs %>% select(taxa, Recombination_Type, category=c) %>%  #Get taxa , Rec type (RG/NRG/NRG-nh), and category
  filter(Recombination_Type %in% c("RG","NRG")) %>% 
  mutate(category = ifelse(nchar(category)==0,"notHit",category)) %>% 
  group_by(taxa, Recombination_Type, category) %>% dplyr::count() %>% ungroup() %>% 
  pivot_wider( names_from = Recombination_Type, values_from = n) %>%
  replace(is.na(.),0)  
 

Fisher_tbl_c = data.frame(NULL); for(i in unique(RG_NRGs$taxa)) { 

# fisher test tbl
# test the first hypothesis RG > NRG
fischers_test_RG_greater_NRG = 
  rstatix::row_wise_fisher_test(count_tbl %>% filter(taxa==i) %>% column_to_rownames("category")  %>% select(RG, NRG) 
                                ,alternative = "greater",  detailed = T, p.adjust.method	= "hochberg") %>% as.data.frame() %>%
  select(n_per_mge_per_category = n, category=group, test_greater=alternative, p_val_greater = p.adj, signifcance_greater = p.adj.signif)

# test the second hypothesis NRG > RG
fischers_test_RG_less_NRG = 
  rstatix::row_wise_fisher_test(count_tbl %>% filter(taxa==i) %>% column_to_rownames("category")  %>% select(RG, NRG) 
                                ,alternative = "less",  detailed = T, p.adjust.method	= "hochberg") %>% as.data.frame() %>% 
  select(n_per_mge_per_category = n,category=group, test_less=alternative, p_val_less = p.adj, signifcance_less = p.adj.signif)

#final fisher tbl
tmp_Fisher_tbl = fischers_test_RG_greater_NRG %>% left_join(fischers_test_RG_less_NRG, by = c("category","n_per_mge_per_category")) %>%
  left_join(count_tbl %>% filter(taxa==i), by = c("category")) %>%
  mutate(sum_RG = sum(RG), sum_NRG=sum(NRG)
         ,fraction_RG = round(RG/sum_RG * 100, 10) 
         ,fraction_NRG = round(NRG/sum_NRG * 100, 10) 
         # consider only significant differences (for plot), keep non-significant (ns) if it is "ns" on both sides
         ,sign_diff = case_when(signifcance_greater != "ns" ~ "greater", signifcance_less != "ns" ~ "less", T ~ "ns")
         # calculate the Difference/sum ratio 
         ,expected_fraction_equal = round( (RG+NRG)/n_per_mge_per_category*100,5)
         ,Diff_ratio_RG = round( (fraction_RG-expected_fraction_equal)/(fraction_RG+expected_fraction_equal),3)
         ,Diff_ratio_NRG = round( (fraction_NRG-expected_fraction_equal)/(fraction_NRG+expected_fraction_equal),3)
                  )
  

Fisher_tbl_c = Fisher_tbl_c %>% bind_rows(tmp_Fisher_tbl)  
}
Fisher_tbl_c = Fisher_tbl_c %>% mutate(Category_group = c)

### safe tbl
Fischer_complete = Fischer_complete %>% bind_rows(Fisher_tbl_c)

}
#safe tbl
#Fischer_complete %>% write_tsv("")

#Plot results
#atler:  filter(NRG>1 & RG>1) %>%  %>%  filter(NRG>10 | RG>10)
Fisher_tbl_plot = Fischer_complete %>% mutate(facet_plot_order=case_when(taxa =="Enterobacterales" ~ "01", 
                                                                         taxa =="Vibrionales" ~ "02",
                                                                         taxa =="Bacillales" ~ "03",
                                                                         taxa =="Mycobacteriales" ~ "04",
                                                                         taxa =="Hyphomicrobiales" ~ "05",
                                                                         taxa =="Lactobacillales" ~ "06",
                                                                         taxa =="Aeromonadales" ~ "07",
                                                                         taxa =="Micrococcales" ~ "08",
                                                                         taxa =="Lachnospirales" ~ "09",
                                                                         taxa =="Alteromonadales" ~ "10",
                                                                         taxa =="Rhodobacterales" ~ "11",
                                                                         taxa =="Eubacteriales" ~ "12",
                                                                         taxa =="Rhodospirillales" ~ "13",
                                                                         taxa =="Sphingomonadales" ~ "14"))
#
### 
plot_obj_tools <- Fisher_tbl_plot %>% filter(category != "notHit"  & signifcance_greater != "ns"& Category_group!="Top_pfam") %>%
    pivot_longer(cols = c(Diff_ratio_RG, Diff_ratio_NRG), names_to = "Recombining_Type", values_to = "DSR") %>%
    ggplot(aes(x=DSR, y = category, fill = Recombining_Type))+
    geom_col(position = "dodge", col = "black") +
    scale_fill_manual(values=c("lightblue","brown3"))+ylab("")+xlab("Diff_Sum_Ratio") +
    theme_bw()+
    theme(legend.position = "none", text = element_text(size = 10), strip.background = element_blank(), panel.spacing = unit(0.5, "cm"))+
    scale_y_discrete(limits=rev)+xlim(-1,1) +
    facet_grid(Category_group~facet_plot_order, scales = "free", labeller = labeller( facet_plot_order=c("01"="Enterobacterales",
                                                                                           "02"="Vibrionales",
                                                                                           "03"="Bacillales",
                                                                                           "04"="Mycobacteriales",
                                                                                           "05"="Hyphomicrobiales",
                                                                                           "06"="Lactobacillales",
                                                                                           "07"="Aeromonadales",
                                                                                           "08"="Micrococcales",
                                                                                           "09"="Lachnospirales",
                                                                                           "10"="Alteromonadales",
                                                                                           "11"="Rhodobacterales",
                                                                                           "12"="Eubacteriales",
                                                                                           "13"="Rhodospirillales",
                                                                                          "14"="Sphingomonadales")))

plot_obj_pfam <- Fisher_tbl_plot %>% filter(category != "notHit"  & Category_group=="Top_pfam") %>% #& signifcance_greater != "ns"
  pivot_longer(cols = c(Diff_ratio_RG, Diff_ratio_NRG), names_to = "Recombining_Type", values_to = "DSR") %>%
  ggplot(aes(x=DSR, y = category, fill = Recombining_Type))+
  geom_col(position = "dodge", col = "black") +
  scale_fill_manual(values=c("lightblue","brown3"))+ylab("")+xlab("Diff_Sum_Ratio") +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 10), strip.background = element_blank(), panel.spacing = unit(0.5, "cm"))+
  scale_y_discrete(limits=rev)+xlim(-1,1) +
  facet_grid(Category_group~facet_plot_order, scales = "free", labeller = labeller( facet_plot_order=c("01"="Enterobacterales",
                                                                                                       "02"="Vibrionales",
                                                                                                       "03"="Bacillales",
                                                                                                       "04"="Mycobacteriales",
                                                                                                       "05"="Hyphomicrobiales",
                                                                                                       "06"="Lactobacillales",
                                                                                                       "07"="Aeromonadales",
                                                                                                       "08"="Micrococcales",
                                                                                                       "09"="Lachnospirales",
                                                                                                       "10"="Alteromonadales",
                                                                                                       "11"="Rhodobacterales",
                                                                                                       "12"="Eubacteriales",
                                                                                                       "13"="Rhodospirillales",
                                                                                                       "14"="Sphingomonadales")))


full_plot_enr <- grid.arrange(plot_obj_tools,plot_obj_pfam, heights=c(1/15,14/15))


write.table(Fischer_complete, file="Fisher_enrichment_pipolins.tsv", sep="\t", row.names = FALSE)

svglite("Full_dif_sum_ratio.svg", width = 15, height = 90)
plot(full_plot_enr)
invisible(dev.off())



