library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)

### Parsing MSF results ###
#results_Aeromonadales_2024-07-02_13.16.18/135624
#results_Vibrionales_2024-07-02_12.56.41/135623
#results_Enterobacterales_2024-07-02_16.24.48/91347

#results_Enterobacter_2024-07-02_12.34.31/547

parse_msf_results <- function(results_file, taxa) {
  
  msf_results <- read.table(results_file, sep="\t", header=TRUE) %>%
    mutate(pipolin_id = paste(genome, sys_id, sep = "-")) %>%
    group_by(pipolin_id) %>%
    mutate(piPolB_in_pipolin = if_else(any(gene_name == "Clu1-GpiPolB"), "Yes", "No"),
           pipolin_gene_n = n_distinct(hit_id),
           pipolin_comp = paste(sort(unique(hit_gene_ref)), collapse = ","),
           pipolin_synt = paste(hit_gene_ref, collapse = ",")) %>%
    ungroup() %>%
    group_by(genome) %>%
    mutate(pipolin_in_genome = if_else(any(gene_name == "Clu1-GpiPolB"), "Yes", "No"),
           sys_n = n_distinct(sys_id),
           n_piPolB = sum(gene_name == "Clu1-GpiPolB")) %>%
    ungroup() %>%
    group_by(genome) %>%
    mutate(piPolB_genome_pipolin = paste(pipolin_in_genome, piPolB_in_pipolin, sep = ","),
           taxa = taxa)
  
  return(msf_results)
  
}

msf_results_ent <- parse_msf_results("results_Enterobacterales_2024-07-02_16.24.48/91347_pipolinMSF_search_results.tsv", "Enterobacterales")
msf_results_vib <- parse_msf_results("results_Vibrionales_2024-07-02_12.56.41/135623_pipolinMSF_search_results.tsv", "Vibrionales")
msf_results_aer <- parse_msf_results("results_Aeromonadales_2024-07-02_13.16.18/135624_pipolinMSF_search_results.tsv", "Aeromonadales")

msf_results_gamma <- rbind(msf_results_ent, msf_results_vib, msf_results_aer)


#Get genome numbers and type of systems
msf_genomes_gamma <- msf_results_gamma %>% distinct(genome, .keep_all = TRUE) %>% group_by(pipolin_in_genome, taxa) %>% summarize(count = n())

msf_system_distrib_gamma <- msf_results_gamma %>% distinct(pipolin_id, .keep_all = TRUE) %>% group_by(genome) %>%
                                                  mutate(piPolB_sys_distrib = paste0(sys_n, paste(sort(piPolB_in_pipolin), collapse = ","), sep="-"))
msf_system_distrib_gamma_n <- msf_system_distrib_gamma  %>% distinct(genome, .keep_all = TRUE) %>% group_by(piPolB_sys_distrib, taxa) %>% summarize(count = n())
#Distrib yes/no plot
msf_system_distrib_gamma_n_p <- ggplot(data=msf_system_distrib_gamma_n, aes(x=piPolB_sys_distrib, y=count, fill=taxa)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("#8a2b32","#eb3f4d","#ffa947")) +
  coord_flip() +
  facet_grid(~taxa, scales = "free")
msf_system_distrib_gamma_n_p

#line plot
msf_system_line_gamma_p <- ggplot(msf_system_distrib_gamma %>% filter(sys_n >1), aes(x = hit_pos, y = genome)) +
  geom_point(aes(color=piPolB_in_pipolin)) + 
  geom_segment(aes(x = min(hit_pos), xend = max(hit_pos), y = genome, yend = genome, color=taxa), alpha = 0.3) +
  scale_color_manual(values=c("darkred","blue4","black", "darkgreen", "red")) +
  theme_minimal() +
  facet_grid(taxa~., scales = "free", space="free") +
  labs(x = "Position", y = "Genome", title = "Genome Positions")
msf_system_line_gamma_p

#Check system composition
msf_pipolin_gamma <- msf_results_gamma %>% distinct(pipolin_id, .keep_all = TRUE)
msf_pipolins_comp_n <- msf_pipolin_gamma %>% group_by(piPolB_in_pipolin, pipolin_comp,taxa) %>% summarize(count = n())
msf_pipolins_comp_n_p <- ggplot(data=msf_pipolins_comp_n, aes(x=pipolin_comp, y=count, fill=piPolB_in_pipolin)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkred")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_pipolins_comp_n_p

#and synteny
msf_pipolins_synt_n <- msf_pipolin_gamma %>% group_by(piPolB_in_pipolin, pipolin_synt,taxa) %>% summarize(count = n())
msf_pipolins_synt_n_p <- ggplot(data=msf_pipolins_synt_n, aes(x=pipolin_synt, y=count, fill=piPolB_in_pipolin)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkred")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_pipolins_synt_n_p

#and synteny only more frequent (count > X)
msf_pipolins_synt_n <- msf_pipolin_gamma %>% group_by(piPolB_in_pipolin, pipolin_synt,taxa) %>% summarize(count = n())
msf_pipolins_synt_n_p <- ggplot(data=msf_pipolins_synt_n %>% filter(count > 1), aes(x=pipolin_synt, y=count, fill=piPolB_in_pipolin)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkred")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_pipolins_synt_n_p

#get synteny nums
#Enterobacterales: #ADD TAG
msf_pipolins_synt_n_ent_pipolin <- msf_pipolins_synt_n %>% filter(pipolin_synt == "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu1-GpiPolB,Clu11-IntP2" | 
                                                                  pipolin_synt == "Clu11-IntP2,Clu1-GpiPolB,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                  taxa == "Enterobacterales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Complete system")
msf_pipolins_synt_n_ent_onlypiPolBloss <- msf_pipolins_synt_n %>% filter(pipolin_synt == "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2" | 
                                                                         pipolin_synt == "Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                         taxa == "Enterobacterales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Hallmark deletion")
msf_pipolins_synt_n_ent_rest_pipolb <- msf_pipolins_synt_n %>% filter(pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu1-GpiPolB,Clu11-IntP2" & 
                                                               pipolin_synt != "Clu11-IntP2,Clu1-GpiPolB,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                               pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2" & 
                                                               pipolin_synt != "Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                               piPolB_in_pipolin == "Yes",
                                                               taxa == "Enterobacterales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Rearranged with hallmark")
msf_pipolins_synt_n_ent_rest_no_pipolb <- msf_pipolins_synt_n %>% filter(pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu1-GpiPolB,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu1-GpiPolB,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                        pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                        piPolB_in_pipolin == "No",
                                                                      taxa == "Enterobacterales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Rearranged without hallmark")


#Vibrionales:
msf_pipolins_synt_n_vib_pipolin <- msf_pipolins_synt_n %>% filter(pipolin_synt == "Clu1-GpiPolB,Clu53-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" | 
                                                                    pipolin_synt == "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu53-IntP2,Clu1-GpiPolB" |
                                                                    pipolin_synt == "Clu1-GpiPolB,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" | 
                                                                    pipolin_synt == "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu1-GpiPolB" &
                                                                    taxa == "Vibrionales")  %>% mutate(sys_name = "1 - Pipolins", sys_class = "Complete system")
msf_pipolins_synt_n_vib_onlypiPolBloss <- msf_pipolins_synt_n %>% filter(pipolin_synt == "Clu53-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" | 
                                                                    pipolin_synt == "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu53-IntP2" |
                                                                    pipolin_synt == "Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" | 
                                                                    pipolin_synt == "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787" &
                                                                    taxa == "Vibrionales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Hallmark deletion")

msf_pipolins_synt_n_vib_rest_pipolb <- msf_pipolins_synt_n %>% filter(pipolin_synt != "Clu1-GpiPolB,Clu53-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" & 
                                                                    pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu53-IntP2,Clu1-GpiPolB" &
                                                                    pipolin_synt != "Clu1-GpiPolB,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" & 
                                                                    pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu1-GpiPolB" &
                                                                    pipolin_synt != "Clu53-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" &
                                                                    pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu53-IntP2" &
                                                                    pipolin_synt != "Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" & 
                                                                    pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787" &
                                                                      piPolB_in_pipolin == "Yes",
                                                                           taxa == "Vibrionales")  %>% mutate(sys_name = "1 - Pipolins", sys_class = "Rearranged with hallmark")
msf_pipolins_synt_n_vib_rest_no_pipolb <- msf_pipolins_synt_n %>% filter(pipolin_synt != "Clu1-GpiPolB,Clu53-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" & 
                                                                        pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu53-IntP2,Clu1-GpiPolB" &
                                                                        pipolin_synt != "Clu1-GpiPolB,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" & 
                                                                        pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu1-GpiPolB" &
                                                                        pipolin_synt != "Clu53-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" &
                                                                        pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu53-IntP2" &
                                                                        pipolin_synt != "Clu3-DUF2787,Clu29-Hyp2,Clu4-Hyp1,Clu32-IntSXT" & 
                                                                        pipolin_synt != "Clu32-IntSXT,Clu4-Hyp1,Clu29-Hyp2,Clu3-DUF2787" &
                                                                        piPolB_in_pipolin == "No",
                                                                        taxa == "Vibrionales")  %>% mutate(sys_name = "1 - Pipolins", sys_class = "Rearranged without hallmark")


#Aeromonadales
msf_pipolins_synt_n_aer_pipolin <- msf_pipolins_synt_n %>% filter(pipolin_synt == "Clu2-IntSXT,Clu435-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu1-GpiPolB,Clu11-IntP2" | 
                                                                    pipolin_synt == "Clu11-IntP2,Clu1-GpiPolB,Clu3-DUF2787,Clu29-Hyp2,Clu435-Hyp1,Clu2-IntSXT" |
                                                                    pipolin_synt == "Clu1-GpiPolB,Clu32-IntSXT,Clu3-DUF2787,Clu29-Hyp2,Clu763-HypAer," | 
                                                                    pipolin_synt == "Clu763-HypAer,Clu29-Hyp2,Clu3-DUF2787,Clu32-IntSXT,Clu1-GpiPolB" &
                                                                    taxa == "Aeromonadales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Complete system")
msf_pipolins_synt_n_aer_onlypiPolBloss <- msf_pipolins_synt_n %>% filter(pipolin_synt == "Clu2-IntSXT,Clu435-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu11-IntP2" | 
                                                                           pipolin_synt == "Clu11-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu435-Hyp1,Clu2-IntSXT" |
                                                                           pipolin_synt == "Clu32-IntSXT,Clu3-DUF2787,Clu29-Hyp2,Clu763-HypAer," | 
                                                                           pipolin_synt == "Clu763-HypAer,Clu29-Hyp2,Clu3-DUF2787,Clu32-IntSXT" &
                                                                           taxa == "Aeromonadales") %>% mutate(sys_name = "1 - Pipolins", sys_class = "Hallmark deletion")

msf_pipolins_synt_n_aer_rest_pipolb <- msf_pipolins_synt_n %>% filter(pipolin_synt != "Clu2-IntSXT,Clu435-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu1-GpiPolB,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu1-GpiPolB,Clu3-DUF2787,Clu29-Hyp2,Clu435-Hyp1,Clu2-IntSXT" &
                                                                        pipolin_synt != "Clu1-GpiPolB,Clu32-IntSXT,Clu3-DUF2787,Clu29-Hyp2,Clu763-HypAer," & 
                                                                        pipolin_synt != "Clu763-HypAer,Clu29-Hyp2,Clu3-DUF2787,Clu32-IntSXT,Clu1-GpiPolB" &
                                                                        pipolin_synt != "Clu2-IntSXT,Clu435-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu435-Hyp1,Clu2-IntSXT" &
                                                                        pipolin_synt != "Clu32-IntSXT,Clu3-DUF2787,Clu29-Hyp2,Clu763-HypAer," & 
                                                                        pipolin_synt != "Clu763-HypAer,Clu29-Hyp2,Clu3-DUF2787,Clu32-IntSXT" &
                                                                        piPolB_in_pipolin == "Yes",
                                                                      taxa == "Aeromonadales")  %>% mutate(sys_name = "1 - Pipolins", sys_class = "Rearranged with hallmark")
msf_pipolins_synt_n_aer_rest_no_pipolb <- msf_pipolins_synt_n %>% filter(pipolin_synt != "Clu2-IntSXT,Clu435-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu1-GpiPolB,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu1-GpiPolB,Clu3-DUF2787,Clu29-Hyp2,Clu435-Hyp1,Clu2-IntSXT" &
                                                                        pipolin_synt != "Clu1-GpiPolB,Clu32-IntSXT,Clu3-DUF2787,Clu29-Hyp2,Clu763-HypAer," & 
                                                                        pipolin_synt != "Clu763-HypAer,Clu29-Hyp2,Clu3-DUF2787,Clu32-IntSXT,Clu1-GpiPolB" &
                                                                        pipolin_synt != "Clu2-IntSXT,Clu435-Hyp1,Clu29-Hyp2,Clu3-DUF2787,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu3-DUF2787,Clu29-Hyp2,Clu435-Hyp1,Clu2-IntSXT" &
                                                                        pipolin_synt != "Clu32-IntSXT,Clu3-DUF2787,Clu29-Hyp2,Clu763-HypAer," & 
                                                                        pipolin_synt != "Clu763-HypAer,Clu29-Hyp2,Clu3-DUF2787,Clu32-IntSXT" &
                                                                        piPolB_in_pipolin == "No",
                                                                      taxa == "Aeromonadales")  %>% mutate(sys_name = "1 - Pipolins", sys_class = "Rearranged without hallmark")

#Distribution of E-values depending on the presence of piPolB
msf_results_gamma_eval <- msf_results_gamma %>% mutate(hit_i_eval = ifelse(hit_i_eval == 0.0, 1.0e-300, hit_i_eval))
pipolin_d_distrib <- ggplot(data=msf_results_gamma_eval, aes(x=log10(hit_i_eval), fill=piPolB_in_pipolin)) +
  geom_histogram( color=NA, alpha=0.5, position = 'identity', bins=20) +
  scale_fill_manual(values=c("cornflowerblue","brown3")) +
  theme_bw() + 
  facet_grid(taxa~gene_name, scales = "free")
pipolin_d_distrib

pipolin_d_dens <- ggplot(data=msf_results_gamma_eval, aes(x=log10(hit_i_eval), fill=piPolB_in_pipolin)) +
  #geom_histogram(aes(y=..density..), color=NA, alpha=0.2, position = 'identity', bins=10) +
  geom_density(alpha=.4, aes(color=piPolB_in_pipolin), linewidth = 0) +
  scale_fill_manual(values=c("gray30","darkred")) +
  scale_color_manual(values=c(NA,NA)) +
  theme_bw() + 
  facet_grid(taxa~gene_name, scales = "free")
pipolin_d_dens


pipolin_d_distrib_len <- ggplot(data=msf_results_gamma_eval, aes(x=log10(hit_i_eval)/hit_seq_len, fill=piPolB_in_pipolin)) +
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity', bins=20) +
  scale_fill_manual(values=c("gray30","darkred")) +
  facet_grid(taxa~gene_name, scales = "free")
pipolin_d_distrib_len



### NEGATIVE CONTROL MZA ###
msf_results_ent_mza <- parse_msf_results("results_Enterobacterales_2024-07-02_16.24.48/91347_pipolinMSF_mza_search_results.tsv", "Enterobacterales")

#Get genome numbers and type of systems
msf_genomes_ent_mza <- msf_results_ent_mza %>% distinct(genome, .keep_all = TRUE) %>% group_by(pipolin_in_genome, taxa) %>% summarize(count = n())

msf_system_distrib_ent_mza <- msf_results_ent_mza %>% distinct(pipolin_id, .keep_all = TRUE) %>% group_by(genome) %>%
  mutate(piPolB_sys_distrib = paste0(sys_n, paste(sort(piPolB_in_pipolin), collapse = ","), sep="-"))
msf_system_distrib_ent_mza_n <- msf_system_distrib_ent_mza  %>% distinct(genome, .keep_all = TRUE) %>% group_by(piPolB_sys_distrib, taxa) %>% summarize(count = n())

#Distrib yes/no plot
msf_system_distrib_ent_mza_n_p <- ggplot(data=msf_system_distrib_ent_mza_n, aes(x=piPolB_sys_distrib, y=count, fill=taxa)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("darkred")) +
  coord_flip() +
  facet_grid(~taxa, scales = "free")
msf_system_distrib_ent_mza_n_p

#line plot
msf_system_line_ent_mza_p <- ggplot(msf_system_distrib_ent_mza , aes(x = hit_pos, y = genome)) +
  geom_point(aes(color=piPolB_in_pipolin)) + 
  geom_segment(aes(x = min(hit_pos), xend = max(hit_pos), y = genome, yend = genome, color=taxa), alpha = 0.3) +
  theme_minimal() + 
  scale_color_manual(values=c("darkred","black", "red")) +
  labs(x = "Position", y = "Genome", title = "Genome Positions")
msf_system_line_ent_mza_p

#Check system composition
msf_pipolin_ent_mza <- msf_results_ent_mza %>% distinct(pipolin_id, .keep_all = TRUE)
msf_results_ent_mza_comp_n <- msf_pipolin_ent_mza %>% group_by(piPolB_in_pipolin, pipolin_comp,taxa) %>% summarize(count = n())
msf_results_ent_mza_comp_n_p <- ggplot(data=msf_results_ent_mza_comp_n, aes(x=pipolin_comp, y=count, fill=piPolB_in_pipolin)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkred")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_results_ent_mza_comp_n_p

#and synteny
msf_pipolins_ent_mza_synt_n <- msf_pipolin_ent_mza %>% group_by(piPolB_in_pipolin, pipolin_synt,taxa) %>% summarize(count = n())
msf_pipolins_ent_mza_synt_n_p <- ggplot(data=msf_pipolins_ent_mza_synt_n, aes(x=pipolin_synt, y=count, fill=piPolB_in_pipolin)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkred")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_pipolins_ent_mza_synt_n_p

#In numbers:
msf_pipolins_ent_mza_synt_n_ent_pipolin <- msf_pipolins_ent_mza_synt_n %>% filter(pipolin_synt == "Clu1-GpiPolB,Clu94-mza5,Clu110-mza4,Clu90-mza3,Clu91-mza2,Clu87-mza1" | 
                                                                          pipolin_synt == "Clu87-mza1,Clu91-mza2,Clu90-mza3,Clu110-mza4,Clu94-mza5,Clu1-GpiPolB") %>%
                                                                          mutate(sys_name = "2 - mza+piPolB", sys_class = "Complete system")
msf_pipolins_ent_mza_synt_n_ent_pipolbLoss <- msf_pipolins_ent_mza_synt_n %>% filter(pipolin_synt == "Clu94-mza5,Clu110-mza4,Clu90-mza3,Clu91-mza2,Clu87-mza1" | 
                                                                            pipolin_synt == "Clu87-mza1,Clu91-mza2,Clu90-mza3,Clu110-mza4,Clu94-mza5")%>%
                                                                            mutate(sys_name = "2 - mza+piPolB", sys_class = "Hallmark deletion")
msf_pipolins_ent_mza_synt_n_ent_other_piPolB <- msf_pipolins_ent_mza_synt_n %>% filter(pipolin_synt != "Clu1-GpiPolB,Clu94-mza5,Clu110-mza4,Clu90-mza3,Clu91-mza2,Clu87-mza1" & 
                                                                            pipolin_synt != "Clu87-mza1,Clu91-mza2,Clu90-mza3,Clu110-mza4,Clu94-mza5,Clu1-GpiPolB" &
                                                                            pipolin_synt != "Clu94-mza5,Clu110-mza4,Clu90-mza3,Clu91-mza2,Clu87-mza1" & 
                                                                            pipolin_synt != "Clu87-mza1,Clu91-mza2,Clu90-mza3,Clu110-mza4,Clu94-mza5" &
                                                                            piPolB_in_pipolin == "Yes") %>%
                                                                            mutate(sys_name = "2 - mza+piPolB", sys_class = "Rearranged with hallmark")
msf_pipolins_ent_mza_synt_n_ent_other_nopiPolB <- msf_pipolins_ent_mza_synt_n %>% filter(pipolin_synt != "Clu1-GpiPolB,Clu94-mza5,Clu110-mza4,Clu90-mza3,Clu91-mza2,Clu87-mza1" & 
                                                                            pipolin_synt != "Clu87-mza1,Clu91-mza2,Clu90-mza3,Clu110-mza4,Clu94-mza5,Clu1-GpiPolB" &
                                                                            pipolin_synt != "Clu94-mza5,Clu110-mza4,Clu90-mza3,Clu91-mza2,Clu87-mza1" & 
                                                                            pipolin_synt != "Clu87-mza1,Clu91-mza2,Clu90-mza3,Clu110-mza4,Clu94-mza5" &
                                                                            piPolB_in_pipolin == "No") %>%
                                                                            mutate(sys_name = "2 - mza+piPolB", sys_class = "Rearranged without hallmark")


#E-values
msf_results_ent_mza_eval <- msf_results_ent_mza %>% mutate(hit_i_eval = ifelse(hit_i_eval == 0.0, 1.0e-300, hit_i_eval))
pipolb_mza_d_distrib <- ggplot(data=msf_results_ent_mza_eval, aes(x=log10(hit_i_eval), fill=piPolB_in_pipolin)) +
  geom_histogram( color=NA, alpha=0.5, position = 'identity', bins=20) +
  scale_fill_manual(values=c("cornflowerblue","brown3")) +
  theme_bw() + 
  facet_grid(taxa~gene_name, scales = "free")
pipolb_mza_d_distrib

### NEGATIVE CONTROL HHE ###
msf_results_ent_hhe <- parse_msf_results("results_Enterobacterales_2024-07-02_16.24.48/91347_pipolinMSF_hhe_search_results.tsv", "Enterobacterales") 
                 
msf_results_ent_hhe <- msf_results_ent_hhe %>% group_by(pipolin_id) %>%
  mutate(hhe_in_system = if_else(any(gene_name == "Clu28-hhe"), "Yes", "No")) %>%
  ungroup() 

#Get genome numbers and type of systems
msf_genomes_ent_hhe <- msf_results_ent_hhe %>% distinct(genome, .keep_all = TRUE) %>% group_by(hhe_in_system, taxa) %>% summarize(count = n())

msf_system_distrib_ent_hhe <- msf_results_ent_hhe %>% distinct(pipolin_id, .keep_all = TRUE) %>% group_by(genome) %>%
  mutate(hhe_sys_distrib = paste0(sys_n, paste(sort(hhe_in_system), collapse = ","), sep="-"))
msf_system_distrib_ent_hhe_n <- msf_system_distrib_ent_hhe  %>% distinct(genome, .keep_all = TRUE) %>% group_by(hhe_sys_distrib, taxa) %>% summarize(count = n())

#Distrib yes/no plot
msf_system_distrib_ent_hhe_n_p <- ggplot(data=msf_system_distrib_ent_hhe_n, aes(x=hhe_sys_distrib, y=count, fill=taxa)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("darkorange")) +
  coord_flip() +
  facet_grid(~taxa, scales = "free")
msf_system_distrib_ent_hhe_n_p
grid.arrange(msf_system_distrib_gamma_n_p, msf_system_distrib_ent_hhe_n_p)

#line plot
msf_system_line_ent_hhe_p <- ggplot(msf_system_distrib_ent_hhe %>% filter(sys_n >1) , aes(x = hit_pos, y = genome)) +
  geom_point(aes(color=hhe_in_system)) + 
  geom_segment(aes(x = min(hit_pos), xend = max(hit_pos), y = genome, yend = genome, color=taxa), alpha = 0.3) +
  theme_minimal() + 
  scale_color_manual(values=c("darkred","black", "darkorange")) +
  labs(x = "Position", y = "Genome", title = "Genome Positions")
msf_system_line_ent_hhe_p


#Check system composition
msf_pipolin_ent_hhe <- msf_results_ent_hhe %>% distinct(pipolin_id, .keep_all = TRUE)
msf_results_ent_hhe_comp_n <- msf_pipolin_ent_hhe %>% group_by(hhe_in_system, pipolin_comp,taxa) %>% summarize(count = n())
msf_results_ent_hhe_comp_n_p <- ggplot(data=msf_results_ent_hhe_comp_n, aes(x=pipolin_comp, y=count, fill=hhe_in_system)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkorange")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_results_ent_hhe_comp_n_p

#and synteny
msf_pipolins_ent_hhe_synt_n <- msf_pipolin_ent_hhe %>% group_by(hhe_in_system, pipolin_synt,taxa) %>% summarize(count = n())
msf_pipolins_ent_hee_synt_n_p <- ggplot(data=msf_pipolins_ent_hhe_synt_n, aes(x=pipolin_synt, y=count, fill=hhe_in_system)) +
  geom_bar(stat="identity", color="black", alpha = 0.7, position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=c("gray30","darkorange")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  facet_grid(taxa~., scales = "free", space = "free")
msf_pipolins_ent_hee_synt_n_p

#synteny numbers
msf_hhe_synt_n_ent_pipolin <- msf_pipolins_ent_hhe_synt_n %>% filter(pipolin_synt == "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2,Clu28-hhe" | 
                                                                    pipolin_synt == "Clu28-hhe,Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT")%>%
                                                                    mutate(sys_name = "3 - core+hhe", sys_class = "Complete system")
msf_hhe_synt_n_ent_onlyhheloss <- msf_pipolins_ent_hhe_synt_n %>% filter(pipolin_synt == "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2" | 
                                                                           pipolin_synt == "Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT") %>%
                                                                    mutate(sys_name = "3 - core+hhe", sys_class = "Hallmark deletion")
msf_pipolins_synt_n_ent_rest_hhe <- msf_pipolins_ent_hhe_synt_n %>% filter(pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2,Clu28-hhe" & 
                                                                        pipolin_synt != "Clu28-hhe,Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                        pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2" & 
                                                                        pipolin_synt != "Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                          hhe_in_system == "Yes") %>%
                                                                        mutate(sys_name = "3 - core+hhe", sys_class = "Rearranged with hallmark")
msf_pipolins_synt_n_ent_rest_no_hhe <- msf_pipolins_ent_hhe_synt_n %>% filter(pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2,Clu28-hhe" & 
                                                                           pipolin_synt != "Clu28-hhe,Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                           pipolin_synt != "Clu2-IntSXT,Clu4-Hyp1,Clu14-Hyp2,Clu3-DUF2787,Clu11-IntP2" & 
                                                                           pipolin_synt != "Clu11-IntP2,Clu3-DUF2787,Clu14-Hyp2,Clu4-Hyp1,Clu2-IntSXT" &
                                                                             hhe_in_system == "No") %>%
                                                                           mutate(sys_name = "3 - core+hhe", sys_class = "Rearranged without hallmark")



#E-values
msf_results_ent_hhe_eval <- msf_results_ent_hhe %>% mutate(hit_i_eval = ifelse(hit_i_eval == 0.0, 1.0e-300, hit_i_eval))
pipolin_d_distrib_hhe <- ggplot(data=msf_results_ent_hhe_eval, aes(x=log10(hit_i_eval), fill=hhe_in_system)) +
  geom_histogram( color=NA, alpha=0.5, position = 'identity', bins=20) +
  scale_fill_manual(values=c("cornflowerblue","brown3")) +
  theme_bw() + 
  facet_grid(taxa~gene_name, scales = "free")
pipolin_d_distrib_hhe

### PLOTS FOR SHOWING ###
#System prevalence - pipolins
label_order <- c(rep("Aeromonadales",2),rep("Enterobacterales",2),rep("Vibrionales",2))
sys_type <- rep("1 - Pipolins",6)
result_col <- rep(c("Yes","No"),3)
N_genomes <- c(307,307,9496,9496,448,448)
N_systems <- c(82,225,383,9113,52,396)
prev_df <- data.frame(order=label_order,System_found=result_col, N_genomes=N_genomes, N_systems=N_systems,
                 prev=N_systems/N_genomes, sys_type=sys_type)
ggplot(prev_df, aes(fill=System_found, y=prev, x=order)) + 
  geom_bar(position="fill", stat="identity") +
  geom_text(aes(x = order, y = prev+0.05, label = N_systems), hjust = 0,color="#f2bfbf")+
  theme_classic() +
  xlab("Order")+
  ylab("N_systems") +
  scale_fill_manual(values=c("gray30","#d92929"), name="System found\nin genome") + coord_flip()

#System prevalence - mza and hhe 
label_order_c <- c(rep("Enterobacterales",4))
sys_type_c <- c(rep("2 - mza+piPolB",2),rep("3 - core+hhe",2))
result_col_c <- rep(c("Yes","No"),1)
N_genomes_c <- rep(9496,4)
N_systems_c <- c(95,9401,381,9115)
prev_df_c <- data.frame(order=label_order_c,System_found=result_col_c, N_genomes=N_genomes_c, N_systems=N_systems_c,
                      prev=N_systems_c/N_genomes_c,sys_type=sys_type_c)
p_prev <-ggplot(rbind(prev_df,prev_df_c), aes(fill=System_found, y=prev, x=order)) + 
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  geom_text(aes(x = order, y = prev+0.05, label = N_systems), hjust = 0, size = 3, color="#f2bfbf")+
  xlab("Order")+
  ylab("Prevalence") +
  scale_fill_manual(values=c("gray30","#d92929"), name="System found\nin genome") + coord_flip()+
  facet_wrap(~sys_type, ncol = 1) +
  ggtitle("System prevalence in GenBank assemblies")
p_prev

#Sys num per genome
msf_system_num_gamma <- msf_results_gamma %>% distinct(pipolin_id, .keep_all = TRUE) %>% 
  group_by(taxa,sys_n) %>% summarise(sys_n_abs = n()) %>% mutate(sys_type="1 - Pipolins")
msf_system_num_mza <- msf_results_ent_mza %>% distinct(pipolin_id, .keep_all = TRUE) %>% 
  group_by(taxa,sys_n) %>% summarise(sys_n_abs = n()) %>% mutate(sys_type="2 - mza+piPolB")
msf_system_num_hhe <- msf_results_ent_hhe %>% distinct(pipolin_id, .keep_all = TRUE) %>% 
  group_by(taxa,sys_n) %>% summarise(sys_n_abs = n()) %>% mutate(sys_type="3 - core+hhe")

p_nG <- ggplot(rbind(msf_system_num_gamma, msf_system_num_mza, msf_system_num_hhe), aes(x=sys_n, y=sys_n_abs, fill=taxa)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic()  +
  geom_text(aes(x = sys_n, y = sys_n_abs+5, label = sys_n_abs), hjust = 0, size = 3, color="gray20")+
  xlab("Systems per genome")+
  ylab("Number of systems") +
  scale_fill_manual(values=c("#8a2b32","#eb3f4d","#ffa947"), name = "Order") + coord_flip() +
  facet_wrap(~sys_type, ncol = 1) +
  ggtitle("Number of systems in a single assembly")

p_sub_1 <- grid.arrange(p_prev, p_nG, nrow=1)
#Synteny
synt_stats_df <- rbind(msf_pipolins_synt_n_ent_pipolin,msf_pipolins_synt_n_ent_onlypiPolBloss,msf_pipolins_synt_n_ent_rest_pipolb,msf_pipolins_synt_n_ent_rest_no_pipolb,
                       msf_pipolins_synt_n_vib_pipolin, msf_pipolins_synt_n_vib_onlypiPolBloss, msf_pipolins_synt_n_vib_rest_pipolb, msf_pipolins_synt_n_vib_rest_no_pipolb,
                       msf_pipolins_synt_n_aer_pipolin, msf_pipolins_synt_n_aer_onlypiPolBloss, msf_pipolins_synt_n_aer_rest_pipolb, msf_pipolins_synt_n_aer_rest_no_pipolb,
                       msf_pipolins_ent_mza_synt_n_ent_pipolin, msf_pipolins_ent_mza_synt_n_ent_pipolbLoss, msf_pipolins_ent_mza_synt_n_ent_other_piPolB, msf_pipolins_ent_mza_synt_n_ent_other_nopiPolB,
                       msf_hhe_synt_n_ent_pipolin, msf_hhe_synt_n_ent_onlyhheloss, msf_pipolins_synt_n_ent_rest_hhe, msf_pipolins_synt_n_ent_rest_no_hhe)

synt_stats_n_grouped <- synt_stats_df %>% group_by(sys_name,sys_class) %>% summarise(sys_class_n = sum(count))
p_synt_bars <- ggplot(synt_stats_n_grouped, aes(x=sys_class, y=sys_class_n, fill=sys_class)) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(x = sys_class, y = sys_class_n+5, label = sys_class_n), hjust = 0, size = 3, color="gray20")+
  theme_classic() +
  theme(legend.position = "none") +
  xlab("System classification")+
  ylab("Number of systems") +
  scale_fill_manual(values=c("brown1","cadetblue2","brown4","cadetblue4")) + coord_flip() +
  facet_wrap(~sys_name, ncol = 1) +
  ggtitle("Distribution of system subcategories")

ggplot(synt_stats_n_grouped, aes(fill=sys_class, y=sys_class_n, x=sys_name)) + 
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  scale_fill_manual(values=c("brown1","cadetblue2","brown4","cadetblue4")) + coord_flip() 




#Stats
synt_stats_n_grouped_wide <- synt_stats_n_grouped %>% pivot_wider(names_from = sys_class, values_from = sys_class_n)
synt_stats_n_grouped_matrix <- as.matrix(synt_stats_n_grouped_wide[,-1])
rownames(synt_stats_n_grouped_matrix) <- synt_stats_n_grouped_wide$sys_name

sys_chi2 <- chisq.test(synt_stats_n_grouped_matrix)
residuals_df <- as.data.frame(sys_chi2$residuals)
residuals_df$sys_name <- rownames(sys_chi2$residuals) 
residuals_df_long <- residuals_df%>%pivot_longer(names_to = "sys_class", values_to = "sys_class_n",cols = -sys_name)
p_resi <- ggplot(residuals_df_long, aes(sys_name, sys_class, fill= sys_class_n)) + 
  geom_tile() +
  xlab("System")+
  ylab("Subcategory") +
  scale_fill_gradient2(low = "cadetblue3", mid = "white", high = "brown3", midpoint = 0, name = "") +
  geom_text(aes(label = round(sys_class_n,3)), color = "black", size = 3) +
  theme_classic() +
  coord_flip() +
  ggtitle("Chi-squared residuals")


p_sub_2 <-grid.arrange(p_synt_bars, p_resi, nrow=1)
p_sub_3 <-grid.arrange(pipolin_d_distrib, pipolb_mza_d_distrib,pipolin_d_distrib_hhe,ncol=1,heights=c(3/6, 1.5/6,1.5/6))
full_plot <- grid.arrange(p_sub_1, p_sub_2, p_sub_3, nrow=3, heights=c(1.7/6, 1.2/6,3/6))
ggsave("MacSyFinder_plots.svg", full_plot, width = 14, height = 20)
