library("data.table")
library("tidyverse")
library("readxl")
library(ggstatsplot)
setDTthreads(30)

#Read MGE info and add MGE type info
MGE_tbl <- fread("../RG_assignation/MGE_info_table.tsv")

#Read gene info
Gene_info <- fread("MGE_Gene_info.tsv") %>% left_join(MGE_tbl, by=c("MGE_id"))



### Include annotations
#1) YR and SR
YR_SR <- fread("YR_and_SR/MGE_pipolin_YR-SR_results.txt") %>% group_by(Gene_ID) %>%
                                                              slice(which.max(Score)) %>% #Keep hit with best score
                                                              ungroup() %>%
                                                              mutate(Rec_class = ifelse(Recombinase %in% c("c3_n1ser", "c2_n1ser", "c1_n1ser"), "SR", 
                                                              ifelse(Recombinase %in% c("Cyan", "Int_P2","Arch2", "Int_SXT","IntKX","Int_Des","Arch1",
                                                              "Int_BPP-1","TnpR","Int_Brujita","Myc","RitC","Integron","RitB","Candidate",
                                                              "Int_CTnDOT","TnpA","RitA","Int_Tn916","Xer"),"YR","Other")),
                                                               Rec_tag = "Recombinase (except IS)") %>% select(-Score)


#2) PADLOC 
PADLOC <- fread("PADLOC/MGE_pipolins_padloc.csv", sep =",") %>% 
          select(Gene_ID=target.name,PADLOC_system = system,Evalue=full.seq.E.value) %>%
          group_by(Gene_ID) %>%
          slice(which.max(Evalue)) %>% #Keep hit with best eval
          ungroup() %>%
          mutate(PADLOC_tag = "Defense") %>%
          select(-Evalue)


#3) PHROG
PHROG_data <- fread("PHROG/phrog_annot_v4.tsv")
PHROG <- fread("PHROG/PHROG_results_table.tsv") %>% mutate(PHROG_tag = "Viral genes",
                                                            phrog = as.numeric(str_remove(PHROG_hit, pattern="^phrog_"))) %>%
                                                     left_join(PHROG_data, by=c("phrog")) %>%
                                                     mutate(PHROG_category = category) %>%
                                                     group_by(Gene_ID) %>%
                                                     slice(which.max(Score)) %>% #Keep hit with best score
                                                     ungroup() %>%
                                                     select(Gene_ID,PHROG_tag,PHROG_category)


#4) AMRFinderPlus 
AMR <- fread("AMRFinder/MGE_pipolins_AMRfinder.tsv", 
             select=c("Protein identifier", "Class","% Identity to reference sequence"),
             col.names = c("Gene_ID", "Class","pident")) %>% 
             group_by(Gene_ID) %>%
             slice(which.max(pident)) %>% #Keep hit with best pident
             ungroup() %>%
             mutate(AMR_tag = "AMR") %>%
             select(Gene_ID, AMR_class = Class, AMR_tag)


#5) VFDB 
VFDB_len <- fread("VFDB/VF_id_info_len.tsv")

VFDB_hits <- fread("VFDB/VFDB_mmseqs/MGE_pipolins_filtered_VFDB_ali.m8", 
                   select = c(1,2,3,4,12), 
                   col.names = c("Gene_ID", "VFID", "pident", "length", "Score")) %>% 
                   filter(pident >= 0.8) %>%
                   left_join(VFDB_len, by=c("VFID")) %>%
                   filter(length/Protein_length >= 0.8)  %>%
                   group_by(Gene_ID) %>%
                   slice(which.max(Score)) %>% #Keep hit with best score
                   ungroup() %>%
                   mutate(VFDB_hit = "Virulence Factor") %>%
                   select(Gene_ID, VF_gene, VF_class, VFDB_hit) 

#6) CONJScan
CONJ <- fread("CONJScan/MGE_pipolin_CONJ_results.txt") %>% group_by(Gene_ID) %>%
                                                            slice(which.max(Score)) %>% #Keep hit with best score
                                                            ungroup() %>% mutate(CONJ_tag="Conjugation") %>% select(-Score)
#...
#END: Join
Gene_annotation = Gene_info %>% left_join(YR_SR, by=c("Gene_ID"))  %>% left_join(PADLOC, by=c("Gene_ID")) %>% left_join(AMR, by=c("Gene_ID"))%>% left_join(VFDB_hits, by=c("Gene_ID"))  %>% left_join(CONJ, by=c("Gene_ID")) %>% left_join(PHROG, by=c("Gene_ID"))

### Write tsv
write_tsv(Gene_annotation, "Gene_annotations.tsv", na="")


