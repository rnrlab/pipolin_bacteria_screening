#0) Requirements
library(data.table)
library(seqinr)
library(tidyverse)
library(randomcoloR)
library(RColorBrewer)
setDTthreads(30)

MGE_tbl = fread("../RG_assignation/MGE_info_table.tsv")
pipolin_info = fread("../Gene_annotation/pipolin_summary_new.tsv")
colnames(pipolin_info)[2] <- "MGE_id"

#1) Import wGRR
wGRR_df_pair = fread("../wGRR_df.tsv", select = c("query", "subject", "wGRR")) %>% 
               mutate(pair_id = paste(query,subject,sep ="_")) %>% 
               select(wGRR,pair_id)
wGRR_taxa <- function(taxa){
    fread("../wGRR_df.tsv", select = c("query", "subject", "wGRR")) %>%
    mutate(MGE_id=query) %>% 
    left_join(pipolin_info,by=c("MGE_id")) %>% 
    filter(order==taxa) %>%
    select(query, subject, wGRR) %>%
    mutate(MGE_id=subject) %>% 
    left_join(pipolin_info,by=c("MGE_id")) %>% 
    filter(order==taxa) %>%
    mutate(pair_id = paste(query,subject,sep ="_")) %>%
    select(query, subject, wGRR, order, pair_id)
}

wGRR_df_pair_ent <- wGRR_taxa("Enterobacterales")
wGRR_df_pair_vibrio <- wGRR_taxa("Vibrionales")
wGRR_df_pair_aerom <- wGRR_taxa("Aeromonadales")
wGRR_df_pair_hyph <- wGRR_taxa("Hyphomicrobiales")
wGRR_df_pair_lac <- wGRR_taxa("Lactobacillales")
wGRR_df_pair_lach <- wGRR_taxa("Lachnospirales")
wGRR_df_pair_mic <- wGRR_taxa("Micrococcales")
wGRR_df_pair_myc <- wGRR_taxa("Mycobacteriales")
wGRR_df_pair_bac <- wGRR_taxa("Bacillales")
  
#2) Import pipolb id values 
pipolb_id = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                  select = c(1,2,3)) %>% 
                  mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
                  left_join(wGRR_df_pair, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_ent = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                         select = c(1,2,3)) %>% 
                   mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
                   left_join(wGRR_df_pair_ent, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_vibrio = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                   select = c(1,2,3)) %>% 
                   mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
                   left_join(wGRR_df_pair_vibrio, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_aerom = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                   select = c(1,2,3)) %>% 
                   mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
                   left_join(wGRR_df_pair_aerom, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_hyph = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                    select = c(1,2,3)) %>% 
                    mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
                    left_join(wGRR_df_pair_hyph, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_lac = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                    select = c(1,2,3)) %>% 
                    mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
                    left_join(wGRR_df_pair_lac, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_lach = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                      select = c(1,2,3)) %>% 
  mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
  left_join(wGRR_df_pair_lach, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_mic = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                       select = c(1,2,3)) %>% 
  mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
  left_join(wGRR_df_pair_mic, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_myc = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                      select = c(1,2,3)) %>% 
  mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
  left_join(wGRR_df_pair_myc, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)
pipolb_id_bac = fread("/home/victor2022/Pipolin_screening/screening_v13/results_Bacteria_genbank_2022-11-11_14.01.58/piPolB_analysis/all_vs_all_pipolb/pipolb_800_DB_ali.m8",
                      select = c(1,2,3)) %>% 
  mutate(pair_id = paste(str_remove(V1, pattern="_[0-9]{1,6}$"),str_remove(V2, pattern="_[0-9]{1,6}$"),sep ="_")) %>% 
  left_join(wGRR_df_pair_bac, by= c("pair_id")) %>% filter(wGRR>0 & V3 >= 0.5)



color_pipolin_groups <- c("#8a2b32", "#25d9c1", "#eb3f4d",
                          "#6ed494", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",   #87ffb3 -> #6ed494 for this plot  
                          "#ffa947") 

pipolb_id_comb <- rbind(pipolb_id_ent, pipolb_id_vibrio,pipolb_id_aerom, pipolb_id_hyph, pipolb_id_lac,
                        pipolb_id_lach, pipolb_id_mic, pipolb_id_myc, pipolb_id_bac)
pipolb_id_test <- rbind(pipolb_id_vibrio,pipolb_id_aerom, pipolb_id_hyph, pipolb_id_lac)

plot_div_pipolb <- ggplot(pipolb_id_comb, aes(x=V3, y=wGRR, group=order)) +
  geom_smooth(aes(color=order),method="auto", linewidth = 0.5) + ylim(0,1) + theme_classic() +
  scale_x_reverse() +
  ggplot2::scale_color_manual(values = color_pipolin_groups)
plot_div_pipolb


#All pipolins      
ggplot(pipolb_id, aes(x=V3, y=wGRR)) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis", trans="log")  +
  geom_smooth(method="auto" , color="white") + ylim(0,1) + xlim(0.5,1)  + theme_classic() +
  scale_x_reverse()





