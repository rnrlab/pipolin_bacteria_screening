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
wGRR_df = fread("../wGRR_df.tsv", select = c("query", "subject", "wGRR"))

wGRR_df_test = wGRR_df %>% filter(grepl("G_",query) & grepl("G_",subject))

#2) transform to matrix
wGRR_matrix <- spread(wGRR_df_test, key = subject, value = as.numeric(wGRR))
rownames(wGRR_matrix) <- wGRR_matrix$query
wGRR_matrix$query <- NULL
wGRR_matrix[is.na(wGRR_matrix)] <- 0 # Replace NA values with 0

#3)prepare colSide
wGRR_df_pipolins_row = data.frame("MGE_id"=rownames(wGRR_matrix)) %>% left_join(pipolin_info,by=c("MGE_id"))
order_row = wGRR_df_pipolins_row$order
colors_nested <-wGRR_df_pipolins_row %>% select(order) %>% mutate(color = 
                 ifelse(order=="Enterobacterales","#fa5943",
                 ifelse(order=="Vibrionales", "#db713b",
                 ifelse(order=="Aeromonadales", "#c93849",
                 ifelse(order=="Alteromonadales","#783038",
                 ifelse(order=="Hyphomicrobiales","#55e66d",
                 ifelse(order=="Rhodobacterales","#79a847",
                 ifelse(order=="Rhodospirillales","#cee889",
                 ifelse(order=="Sphingomonadales","#237d2a",
                 ifelse(order=="Mycobacteriales","#792abf",
                 ifelse(order=="Micrococcales","#e655c9",
                 ifelse(order=="Lachnospirales","#d9d75f",
                 ifelse(order=="Eubacteriales","#e6bf55",
                 ifelse(order=="Lactobacillales","#5580e6",
                 ifelse(order=="Bacillales","#48c5d9","gray")))))))))))))))


heatmap(as.matrix(wGRR_matrix), scale="none", labRow=NA, labCol=NA, RowSideColors=colors_nested$color)

