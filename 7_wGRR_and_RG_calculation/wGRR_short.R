library(data.table)
library(seqinr)
library(tidyverse)
setDTthreads(30)

#Load protein fasta from prodigal output
prot_info_lst = read.fasta(seqtype = "AA", as.string = T,
                           file="MGE_pipolins_filtered_encoded.faa")
protein_sizes_df = data.frame(protein_id = unlist(lapply(prot_info_lst,getName)), 
                              prot_length = unlist(lapply(prot_info_lst, getLength))) %>% mutate(protein_id = as.character(protein_id))
# get number of proteins per genome
prot_number = protein_sizes_df %>% distinct() %>% 
  mutate(replicon =  str_remove(protein_id, pattern="_[0-9]{1,6}$")) %>% 
  group_by(replicon) %>% summarise(protein = n()) %>% ungroup()
rm(prot_info_lst, protein_sizes_df)

# load the BBH.tsv shortened file and calculate for each element size the p_ident to size ratio
path_data = "BBH.tsv"
score_df = fread(path_data, select = c(1,2,3), col.names=c("qseqid","sseqid","pident")) %>% 
                    mutate(query = str_remove(qseqid, pattern="_[0-9]{1,6}$")) %>%
                    select(-qseqid) %>%
                    mutate(subject = str_remove(sseqid, pattern="_[0-9]{1,6}$")) %>%
                    select(-sseqid) %>%  
                    left_join(y = prot_number, by= c("query"="replicon"),copy = T) %>% rename(prot_query = protein) %>%
                    left_join(y = prot_number, by= c("subject"="replicon"),copy = T) %>% rename(prot_subject = protein) %>%
                    arrange(subject) %>% 
                    mutate(ident_s = round(pident/prot_subject,5), 
                           ident_q = round(pident/prot_query,5))


#print(nrow(score_df))
#print(score_df)
score_df_g = score_df %>% group_by(query,subject) %>% summarize(score_s =  round(sum(ident_s, na.rm = TRUE),3)
                                                               ,score_q = round(sum(ident_q, na.rm = TRUE),3)
                                                               ,n_similar_proteins =n()) %>% ungroup() 

#write_tsv(score_df_g,"Score_df.tsv")
rm(score_df)
#Removed from mutate:
#edge_id_s = paste(subject,query, sep ="_")
#edge_id_q = paste(query,subject, sep ="_")

# attach protein numbers again and take only the score of the element with fewest genes/proteins
wGRR_df  = score_df_g  %>% 
  left_join(y = prot_number, by= c("query"="replicon"),copy = T) %>% 
  rename(prot_query = protein) %>%
  left_join(y = prot_number, by= c("subject"="replicon"),copy = T) %>% rename(prot_subject = protein) %>%
  mutate(wGRR = ifelse(prot_subject < prot_query, score_s , score_q) ) %>% 
  select(query, subject, wGRR,n_similar_proteins,  prot_query,prot_subject)
rm(score_df_g)
write_tsv(wGRR_df, "wGRR_df_encoded.tsv")
#revert encoding

mge_codes <- fread("mge_codes.tsv", col.names = c("MGE_id","code"))
wGRR_df_names <- wGRR_df %>% left_join(mge_codes, by=c("query"="code")) %>%
                             mutate(query_id = MGE_id) %>% select(-query,-MGE_id) %>% 
                             left_join(mge_codes, by=c("subject"="code")) %>%
                             mutate(subject_id = MGE_id) %>% select(-subject,-MGE_id) %>%
                             select("query_id","subject_id","wGRR","n_similar_proteins","prot_query","prot_subject")
colnames(wGRR_df_names)[1] <- c("query")
colnames(wGRR_df_names)[2] <- c("subject")

write_tsv(wGRR_df_names, "wGRR_df.tsv")
