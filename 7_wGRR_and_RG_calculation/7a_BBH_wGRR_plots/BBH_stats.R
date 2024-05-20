#0) Requirements
library(data.table)
library(seqinr)
library(tidyverse)
library(svglite)
setDTthreads(30)

#1) Import wGRR
wGRR = fread("../wGRR_df_encoded.tsv", select = c("query", "subject", "wGRR", "prot_query", "prot_subject")) %>% 
       mutate(pair_id = paste(query,subject,sep ="_")) %>% select(wGRR,pair_id)


#2) Get BBH between MGEs 
BBH_id = fread("../BBH.tsv", select =  c("qseqid", "sseqid", "pident")) %>% 
                  mutate(pair_id = paste(str_remove(qseqid, pattern="_[0-9]{1,6}$"),str_remove(sseqid, pattern="_[0-9]{1,6}$"),sep ="_")) %>%
                  left_join(y = wGRR, by= c("pair_id")) %>% select(-qseqid, -sseqid)

gc()
#3) 2D histogram
#data_test <- head(BBH_id, n=1000000)
heatmap_bbh <- ggplot(BBH_id, aes(x=pident, y=wGRR)) +
  geom_bin2d(bins = 75) +
  scale_fill_gradient2('n', low = "white", mid = "#ff9e9e", high = "darkred", midpoint = 10000000) +
  theme_classic()
heatmap_bbh

#4) Density
BBH_id80 = fread("../BBH.tsv", select =  c("qseqid", "sseqid", "pident")) %>% 
  mutate(pair_id = paste(str_remove(qseqid, pattern="_[0-9]{1,6}$"),str_remove(sseqid, pattern="_[0-9]{1,6}$"),sep ="_")) %>%
  filter(pident>=0.8) %>% left_join(y = wGRR, by= c("pair_id")) %>% select(wGRR)

BBH_id80_cov80 = fread("../BBH.tsv", select =  c("qseqid", "sseqid", "pident","cov_q","cov_s")) %>% 
  mutate(pair_id = paste(str_remove(qseqid, pattern="_[0-9]{1,6}$"),str_remove(sseqid, pattern="_[0-9]{1,6}$"),sep ="_")) %>%
   filter(pident>=0.8, cov_q>=0.8, cov_s>=0.8) %>% left_join(y = wGRR, by= c("pair_id")) %>% select(wGRR)

BBH_id80_comp <- rbind(BBH_id80 %>% mutate(tag="80% id"),
                       BBH_id80_cov80 %>% mutate(tag="80% id + 80% cov"))

rm(BBH_id80, BBH_id80_cov80)
gc()

hist_bbh <- ggplot(BBH_id80_comp,  aes(x=wGRR, fill=tag)) +
            geom_histogram(color=NA, bins=100, position = 'identity') +
            scale_fill_manual(values=c("#ff9e9e", "darkred")) +
            labs(fill="") +
            theme_classic() +
            theme(legend.position.inside=c(0,5000000))

svglite("histogram.svg", width = 12, height = 6)
hist_bbh
invisible(dev.off())
