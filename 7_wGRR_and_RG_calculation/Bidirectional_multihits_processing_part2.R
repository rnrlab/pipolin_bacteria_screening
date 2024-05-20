library(data.table)
library(seqinr)
library(tidyverse)
setDTthreads(30)


names_lst = c("qseqid", "sseqid", "pident", "cov_q", "cov_s", "bitscore")

solved_df =  fread("BBH_multi_hit_ids_sorted_qquery_solvedLC.tsv", col.names = names_lst) %>% mutate( 
                         query = str_remove(qseqid, pattern="_[0-9]{1,6}$"),
                         subject = str_remove(sseqid, pattern="_[0-9]{1,6}$"),
                         pair_id = paste(query, subject, sep = "_")) %>% 
                         group_by(qseqid,subject,pair_id) %>% 
                         arrange(desc(bitscore),desc(pident),desc(cov_q)) %>% 
                         mutate(n_q=1:n()) %>%
                         group_by(sseqid,query,pair_id) %>% 
                         arrange(desc(bitscore),desc(pident),desc(cov_s)) %>%
                         mutate(n_s=1:n())


solved_df_out = solved_df %>% ungroup() %>% select(-pair_id, -bitscore, -query, -subject) %>% filter(n_q ==1 & n_s == 1) %>% 
                                            select(-n_q, -n_s)

write_tsv(solved_df_out, "BBH_multi_hit_solved.tsv")

