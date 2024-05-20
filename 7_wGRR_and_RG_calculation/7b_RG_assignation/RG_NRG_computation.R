#0) Requirements
library(data.table)
library(seqinr)
library(tidyverse)
setDTthreads(30)

#1: Get genes with BBH (RG+NRG)
BBH_ids <- fread("../BBH.tsv", select=c("qseqid", "sseqid")) 
mge_codes <- fread("../mge_codes.tsv", col.names = c("MGE_id","code"))

BBH_ids_list <- BBH_ids %>% select(Gene_ID=qseqid) %>% distinct() %>% bind_rows(BBH_ids%>% ungroup %>% select(Gene_ID=sseqid)) %>% distinct() 

rm(BBH_ids)

BBH_ids_names = BBH_ids_list %>%  mutate(MGE_code = str_remove(Gene_ID, pattern="_[0-9]{1,6}$")) %>% 
                            left_join(mge_codes, by=c("MGE_code"="code")) %>%
                            mutate(Gene_ID_mge = paste(MGE_id, gsub("^.*?_", "", Gene_ID), sep = "_")) %>% 
                            select(Gene_ID = Gene_ID_mge) 


#2: Get RGs
RGs_ids <- fread("Pipolins_MGEs_RefSeq_RGs_20.tsv", select=c("qseqid", "sseqid"))
RGs_ids_names <- RGs_ids %>% select(Gene_ID=qseqid) %>% distinct() %>% bind_rows(RGs_ids%>% ungroup %>% select(Gene_ID=sseqid)) %>% distinct() %>% mutate(Recombination_Type="RG")

#3: Compute NRGs
NRGs = BBH_ids_names %>% anti_join(RGs_ids_names,by="Gene_ID") %>% mutate(Recombination_Type="NRG")


#5: Compute nh-NRG
Gene_info <- fread("../Gene_annotation/MGE_Gene_info_filtered.tsv", select = c("Gene_ID"))
nh_NRGs = Gene_info %>% anti_join(BBH_ids_names,by="Gene_ID")  %>%  mutate(Recombination_Type="nhNRG")


#6: Get all CDS
Gene_info_full <- bind_rows(RGs_ids_names, NRGs, nh_NRGs)
write_tsv(Gene_info_full, "RG_NRG_list_20.tsv")

