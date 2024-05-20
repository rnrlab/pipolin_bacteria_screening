#0) Requirements
library(data.table)
library(seqinr)
library(tidyverse)
setDTthreads(30)

#1) Import wGRR
wGRR = fread("../wGRR_df_encoded.tsv", select = c("query", "subject", "wGRR", "prot_query", "prot_subject")) %>% 
       mutate(pair_id = paste(query,subject,sep ="_")) %>% select(wGRR,pair_id)


#3) Get BBH between MGEs with wGRR < 0.2
BBH_id80_10 = fread("../BBH.tsv", select =  c("qseqid", "sseqid", "pident", "cov_q", "cov_s")) %>% 
                  filter(pident >= 0.8 & cov_q >= 0.8 & cov_s >= 0.8) %>%
                  mutate(pair_id = paste(str_remove(qseqid, pattern="_[0-9]{1,6}$"),str_remove(sseqid, pattern="_[0-9]{1,6}$"),sep ="_")) %>%
                  left_join(y = wGRR, by= c("pair_id")) %>% 
                  filter(wGRR < 0.3)

#4) Apply filter n_RG
RGs = BBH_id80_10 %>%   arrange(pair_id) %>%
                        group_by(pair_id) %>% 
                        mutate(n_RG=dplyr::n()) %>% 
                        ungroup() %>%
                        filter(n_RG < 25) #Filter pairs with more than 25 hits

rm(BBH_id80_10)
#5: Rever to original name
mge_codes <- fread("../mge_codes.tsv", col.names = c("MGE_id","code"))
RGs_out = RGs %>% mutate(query_code = str_remove(qseqid, pattern="_[0-9]{1,6}$"),
                                   subject_code = str_remove(sseqid, pattern="_[0-9]{1,6}$")) %>% 
                                  left_join(mge_codes, by=c("query_code"="code")) %>%
                                  mutate(query= MGE_id) %>% select(-query_code,-MGE_id) %>% 
                                  left_join(mge_codes, by=c("subject_code"="code")) %>%
                                  mutate(subject = MGE_id) %>% select(-subject_code,-MGE_id) %>%
                                  mutate(qseqid_mge = paste(query, gsub("^.*?_", "", qseqid), sep = "_"),
                                         sseqid_mge = paste(subject, gsub("^.*?_", "", sseqid), sep = "_")) %>%
                                  select(qseqid= qseqid_mge, sseqid = sseqid_mge, query, subject, wGRR)

write_tsv(RGs_out, "Pipolins_MGEs_RefSeq_RGs_30.tsv")




