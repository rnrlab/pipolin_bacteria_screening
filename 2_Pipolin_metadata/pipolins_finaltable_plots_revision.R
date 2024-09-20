#####################################################################################################
############     Step 0: Load packages & data  ######################################################
#####################################################################################################

paquetes <- c("ggplot2","data.table","shadowtext", "dplyr","tidyverse","patchwork","webr","ggridges","stringr")
unavailable <- setdiff(paquetes, rownames(installed.packages()))
install.packages(unavailable)
lapply(paquetes, library, character.only = TRUE)

# load data
pipolins_full <- fread("tables/pipolin_summary_postprocessed_new.tsv")
pipolins <- fread("tables/pipolin_bacteria_stats.tsv")
#prepare data
pipolins$family[is.na(pipolins$family)] <- "Unknown"
pipolins$order[is.na(pipolins$order)] <- "Unknown"


#####################################################################################################
############     Figure 1: Distribution of pipolin in bacterial genome assemblies  ##################
#####################################################################################################

#aggregate by order
O <- aggregate(pipolins~order, data=pipolins, sum, na.action = 'na.pass')
O$text <- paste0(O$order," ",round(O$pipolins*100/sum(O$pipolins),1)," %")

for (i in 1:nrow(pipolins)){
  pipolins$order2[i] <-as.character(O[grepl(pipolins$order[i],O$text),3])
}
#limit to 20 pipolins per order
O <- O[O$pipolins >20,]
summary <- subset(pipolins, pipolins > 0 & order %in% as.vector(O$order) )

#count by genus
summary$count <- paste0(summary$genus, " (", summary$pipolins, ")")
#color list
color_pipolin_orders <- c("#8a2b32", "#fa7223", "#25d9c1", "#eb3f4d", "#c782ff",
                          "#7FE072", "#8e2abf", "#4781e6", "#e32bac","#ff94f1",
                          "#5bad51", "#ffa947","grey40")
names(color_pipolin_orders) <- levels(as.factor(summary$order))

#prevalence plot 
summary <- summary %>% arrange(class,order,family,genus)

ggplot(summary[summary$pipolins>3,], aes(x = factor(genus,levels=genus), y = prevalence100, group=order)) +
  geom_bar(aes(fill = order), stat = "identity", col = "black") + 
  scale_y_continuous(expand = c(0,0), limits=c(0,100), breaks=c(20,40,60,80,100)) +
  geom_text(aes(label=paste0(pipolins," (",prevalence100,"%)")),size = 3, position = position_dodge(width=0.9), hjust=-0.05,  angle=90, color="grey40", fontface = "bold") + 
  facet_grid(.~fct_reorder(order, class), scale="free_x",space="free")+
  scale_fill_manual(values=color_pipolin_orders) + theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "italic")) + theme(plot.margin = margin(t = 0.5,  r = 0, b = 0, l = 0.4, unit = "cm"),legend.text=element_text(face="italic"),legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.text.x.top =   element_blank(),  panel.grid.major = element_line(size = 0.5, linetype = 'solid',  colour = "white"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "white")) +
  xlab("Bacteria Genus")+ ylab("Pipolins prevalence %")+labs(fill = "Order")
ggsave("figure_1_dual.pdf",width=12,height=6)

#####################################################################################################
###################       Figure 2: Pipolins reconstruction and structure   #########################
#####################################################################################################

#Figure 2A
#prepare the data
num_pipolins <- as.data.frame(table(pipolins_full$Genome_ID))
names(num_pipolins) <- c("Genome_ID","num_pipolins")
pipolins_full <- merge(pipolins_full,num_pipolins)

pipolins_order <- as.data.frame(table(pipolins_full$class, pipolins_full$order,pipolins_full$num_pipolins))
pipolins_order <- pipolins_order[pipolins_order$Freq>20 & pipolins_order$Var1 !="",]
pipolins_order$Var2 <- fct_reorder(factor(pipolins_order$Var2),as.numeric(factor(pipolins_order$Var1)), .desc=TRUE)

n_o_pipolins <- ggplot(pipolins_order[pipolins_order$Freq>20 & pipolins_order$Var1 !="",], aes(fill=Var2, y=Freq, x=Var3)) +   geom_bar(position="stack", stat="identity") +  scale_fill_manual(values=color_pipolin_orders, name = "Order") +
  coord_flip() +  scale_y_log10(labels=c()) +
  geom_shadowtext(aes(label=ifelse(Freq>0,Freq,"")),size = 3,angle=90, position= position_stack(vjust = 0.5), color="white", fontface = "bold") +
  xlab("") + ylab("") + ggtitle("Pipolins per assembly")+
  theme_linedraw(base_size = 12) + theme(legend.text=element_text(face="italic"), axis.ticks.x=element_blank())

pipolb_order <- as.data.frame(table(pipolins_full$class,pipolins_full$order,pipolins_full$piPolB_CDS_num))
pipolb_order <- pipolb_order[pipolb_order$Freq>20 & pipolb_order$Var1 !="" & pipolb_order$Var2 !="",]

pipolb_order$Var2 <- fct_reorder(pipolb_order$Var2,as.numeric(pipolb_order$Var1), .desc=TRUE)
n_o_piPolBs <- ggplot(pipolb_order, aes(fill=Var2, y=Freq, x=Var3)) +   geom_bar(position="stack", stat="identity") +  scale_fill_manual(values=color_pipolin_orders, name = "Order") +
  coord_flip() +  scale_y_log10(labels=c()) +
  geom_shadowtext(aes(label=ifelse(Freq>0,Freq,"")),size = 3, angle=90,position = position_stack(vjust = 0.5), color="white", fontface = "bold") +
  xlab("") + ylab("") + ggtitle("piPolBs per pipolin")+
  theme_linedraw(base_size = 12) + theme(legend.text=element_text(face="italic"), axis.ticks.x=element_blank())
#plot
n_o_pipolins / n_o_piPolBs +  plot_layout(ncol = 1, nrow = 2, guides = "collect")
ggsave("figure_2a.pdf",width = 6, height=4)

#Figure 2B
tmp <- as.data.frame(table(pipolins_full$class,pipolins_full$order,pipolins_full$Number_atts,pipolins_full$Att_type),exclude=0)
tmp <- tmp[tmp$Var3!=0 & tmp$Freq!=0,]
tmp$Var3 <- as.character(tmp$Var3)
tmp$Var3[as.integer(tmp$Var3) > 3] <- ">3"

tmp$Var3 <- factor(tmp$Var3, levels=c("1","2","3",">3"))
tmp <- tmp[tmp$Var2 %in% O$order & tmp$Var4=="['de novo']" | tmp$Var2 %in% O$order & tmp$Var4=="['pipolin conserved']",]

tmp <- tmp %>% 
  group_by(Var1,Var2,Var3,Var4) %>% 
  summarize(Freq = sum(Freq),.groups="keep")

ggplot(tmp, aes(fill=fct_reorder(factor(tmp$Var2),as.numeric(tmp$Var1), .desc = TRUE), y=log10(Freq), x=Var3)) +   geom_bar(position="stack", stat="identity") +  scale_fill_manual(values=color_pipolin_orders, name = "Order") +
  coord_flip() +  facet_wrap(~Var4,labeller=as_labeller(c(`['de novo']`="De novo DRs",`['pipolin conserved']`="E. coli-like DRs")),scales="free_x") +    geom_shadowtext(aes(label=ifelse(Freq>2,Freq,"")),size = 3, position = position_stack(vjust = 0.5), color="white", fontface = "bold") +
  xlab("Number of DRs") + ylab("Number of pipolins") + 
  theme_linedraw(base_size = 12) + theme(legend.text=element_text(face="italic"),axis.ticks.x=element_blank(),axis.text.x = element_blank())
ggsave("figure_2b.pdf",width = 6, height=3.7)

#Figure 2C
pipolins_full$Number_atts2 <- ifelse(pipolins_full$Number_atts>2, "≥3",pipolins_full$Number_atts)
cairo_pdf(file="figure_2c.pdf",width = 4, height = 4 )  
PieDonut(setNames(as.data.frame(table(pipolins_full$Number_atts2,pipolins_full$Assembly_gaps_pipolin_reconstruction)),c("DRs","Gaps","Freq")), aes(Gaps, DRs, count=Freq), showDonutName=TRUE, ratioByGroup = FALSE,r0=0,pieLabelSize = 5,
              donutLabelSize = 3) 
dev.off()

#Figure 2D and final table with delimited pipolins
delim <- read.table("data/Pipolin_delimited_new_int_site_length.tsv", sep="\t",header=TRUE)
#remove inconsistencies (very small or anomalous)
delim <- delim[delim$delim_length!=1220,]
delim <- delim[-which(delim$Pipolin_id=="G_1243401_1v0"),]
#keep only "Alt_site"
delim <- delim[delim$tdelim_type=="Alt_site",]
colnames(delim)[1] <- names(pipolins_full)[3]


#import BLAST results and extract coordinates
for (i in 1:nrow(delim)){
  blast <-  read.delim(paste0("data/delim_filtered_pipolins_fna/blast_results/",delim$Pipolin_ID[i],".txt"), comment.char = '#',sep ="\t", header=FALSE)
  blast[,13] <- apply(blast[,9:10], 1, min)   
  blast[,10] <- apply(blast[,9:10], 1, max) 
  blast[,9] <- blast[,13]
  blast <- blast[,-13]
  blast <- blast[order(blast[,3],blast[,8], decreasing = TRUE),] #in case multiple hits sort by identity and 
  if (nrow(blast)==1){
    blast <- sapply(blast, as.character)
    delim$EP_contig[i] <-  blast[2]
    delim$coords[i] <- paste0("['",blast[2]," ",blast[9],":",blast[10],"']")
    delim$coord_dif[i] <- as.numeric(blast[10]) - as.numeric(blast[9])
    
  } else if (as.numeric(blast[1,4])==as.numeric(blast[1,8]) & as.numeric(blast[2,7]<=as.numeric(blast[1,8]))){
    blast <- blast[!duplicated(blast[,1]),] #remove duplicates
    blast <- sapply(blast, as.character)
    delim$EP_contig[i] <-  blast[2]
    delim$coords[i] <- paste0("['",blast[2]," ",blast[9],":",blast[10],"']")
    delim$coord_dif[i] <- as.numeric(blast[10]) - as.numeric(blast[9])
    
  } else{
    blast <- sapply(blast, as.character)
    delim$EP_contig[i] <-  blast[1,2]
    delim$coords[i] <- paste0("['",blast[1,2]," ",min(as.numeric(blast[,9])),":",max(as.numeric(blast[,10])),"']")
    delim$coord_dif[i] <- abs(max(as.numeric(blast[,10])) - min(as.numeric(blast[,9])))
   
  }
}
#when there are many DRs mapping is better by hand/eye from Blastn output:
weird <- delim$Pipolin_ID[which(delim$delim_length!=delim$coord_dif)]
delim[which(delim$Pipolin_ID==weird[1]),7:8] <- c(984033,1002837)
delim[which(delim$Pipolin_ID==weird[1]),5] <- paste0("['",delim[which(delim$Pipolin_ID==weird[1]),4]," ",delim[which(delim$Pipolin_ID==weird[1]),7],":",delim[which(delim$Pipolin_ID==weird[1]),8],"']")
delim[which(delim$Pipolin_ID==weird[2]),7:8] <- c(75120,101524)
delim[which(delim$Pipolin_ID==weird[2]),5] <- paste0("['",delim[which(delim$Pipolin_ID==weird[2]),4]," ",delim[which(delim$Pipolin_ID==weird[2]),7],":",delim[which(delim$Pipolin_ID==weird[2]),8],"']")
delim[which(delim$Pipolin_ID==weird[3]),7:8] <- c(91423,120147)
delim[which(delim$Pipolin_ID==weird[3]),5] <- paste0("['",delim[which(delim$Pipolin_ID==weird[3]),4]," ",delim[which(delim$Pipolin_ID==weird[3]),7],":",delim[which(delim$Pipolin_ID==weird[3]),8],"']")


write.table(delim$EP_contig,"data/alt_contigs_EP.txt",col.names = FALSE,row.names = FALSE, quote=FALSE)
#grep in bash to obtain full contig headers
#add the genbank contig and coords to 'delim' table
alt_headers <- readLines("data/headers_alt_pipolins_gb.txt")
alt_headers <- gsub(">","",alt_headers)
alt_headers <- stringr::str_extract(alt_headers, "[^ ]* [^ ]*")
alt_headers <- read.table(text=alt_headers,sep=" ")
delim$GB_contig <- alt_headers$V2[alt_headers$V1==delim$EP_contig]

for (i in 1:nrow(delim)) {
    delim$GB_coords[i] <- gsub(pattern=delim$EP_contig[i],replacement=delim$GB_contig[i],delim$coords[i])
}

  
  
#combine table
pipolins_delim <- merge(pipolins_full,delim[c(1:4,5,10)],by="Pipolin_ID",all.x=TRUE)
pipolins_delim$Pipolin_length <- ifelse(is.na(pipolins_delim$delim_length),pipolins_delim$Pipolin_length,pipolins_delim$delim_length)
pipolins_delim$Atts <-  ifelse(is.na(pipolins_delim$tdelim_type),pipolins_delim$Atts,pipolins_delim$tdelim_type)
pipolins_delim$Atts[pipolins_delim$order=="Bacillales"] <- "Plasmid"
pipolins_delim$Contig_pipolin_coordinates_ExplorePipolin <- ifelse(is.na(pipolins_delim$GB_coords),pipolins_delim$Contig_pipolin_coordinates_ExplorePipolin,pipolins_delim$GB_coords)

#combine with selected orders & calculate stats
tmp <- merge(pipolins_delim[,c(7,20,39,40)],O[,1:2])
tmp <- tmp[tmp$Atts != "FALSE",]
tmp <- tmp[tmp$pipolins>20, ]

number <- as.data.frame(tmp %>% count(order) )
tmp$order <- fct_reorder(tmp$order,as.numeric(as.factor(tmp$class)),.desc=TRUE)
#plot
ggplot(tmp)+geom_density_ridges(aes(x=Pipolin_length, y=order,fill=order),alpha=0.8) + scale_fill_manual(values=color_pipolin_orders)+scale_x_continuous(labels = scales::scientific,limits=c(0,100000))+ylab("")+xlab("Pipolin length")+
  geom_text(data=number,aes(y=order, x=99000,label=paste0("n= ",n)),color="grey40",vjust=-0.5,hjust=1,inherit.aes = FALSE)+theme_classic(base_size = 14) + theme(axis.text.y = element_text(face = "italic"),legend.position = "none")
ggsave("figure_2d.pdf",width = 6, height=4)


#####################################################################################################
###############   Supplementary Figure 2: piPolB protein length distribution   ######################
#####################################################################################################
#remove the square brackets
pipolins_full$piPolB_CDS_length <- gsub("[?\\]\\[*]", "", pipolins_full$piPolB_CDS_length, perl=T) 
#duplicate rows with more than one piPolB fragment
kk <- pipolins_full %>% separate_rows(piPolB_CDS_length, sep=",")
#combine with selected orders & calculate stats
tmp <- merge(kk[,c(9,40,41)],O[,1:2])
tmp <- tmp[tmp$pipolins>20, ]
number <- data.frame(tmp %>% count(order) ,tmp[as.numeric(tmp$piPolB_CDS_length) > 800,] %>% count(order))[,c(1,2,4)]
number$perc <- round(number$n.1 * 100 / number$n,1)
tmp$order <- fct_reorder(tmp$order,as.numeric(as.factor(tmp$class)),.desc=TRUE)

#plot
ggplot(tmp, aes(x=as.numeric(piPolB_CDS_length), y=order,fill=order,color=order))+
  geom_density_ridges(jittered_points = TRUE, point_shape = 21, point_size = 0.75, point_alpha = 0.4,   alpha=0.5) + 
  scale_fill_manual(values=color_pipolin_orders) + scale_color_manual(values=color_pipolin_orders) +
  scale_x_continuous(limits=c(0,1400),breaks = c(200,400,600,800,1000,1200,1400))+ylab("")+xlab("piPolB length") +
  geom_text(data=number,aes(y=order, x=1250,label=paste0("n= ",n, " (", perc, "%)")),color="grey40",vjust=-0.5,hjust=0,inherit.aes = FALSE) + 
  theme_classic(base_size = 14) + theme(axis.text.y = element_text(face = "italic"),legend.position = "none")+
  geom_vline(xintercept=800, linetype = "dashed") + 
  theme(plot.margin = margin(0.5,1.75,0,0, "cm")) + coord_cartesian(clip = "off")
ggsave("pipolb_length.pdf",,width=6,height=8)




#####################################################################################################
#########################    ADDING OriTs and WRITTING FINAL TABLE     ##############################  
#####################################################################################################
#OriTs
oriT <- fread("~/Data/oriT/test_oriT.tsv")
kk <- oriT %>% group_by(V7, V8, V1,V11) %>% reframe(text = str_c(V2, collapse = ", "))
#keep only one hit of each oriT per pipolin
kk <- as.data.frame(kk[!duplicated(kk[,c(3,5)]),c(3,5)])
names(kk) <- c("Pipolin_ID","oriT")
pipolins_delim <- merge(pipolins_delim,kk,all.x=TRUE)


#write table
write.table(unique(pipolins_delim[,c(2,1,7,18,8,9,10,20,21,22,23,24,25,27,28,29,30,32,34,35,43,42,41,40,39,38)]),"data/Supp_dataset2_TBS.tsv",quote=FALSE,row.names = FALSE,sep="\t")

