paquetes <- c("ggplot2","data.table","shadowtext", "dplyr","tidyverse","patchwork","webr","ggridges")
unavailable <- setdiff(paquetes, rownames(installed.packages()))
install.packages(unavailable)
lapply(paquetes, library, character.only = TRUE)


#1: load full data and prepare for plots
pipolins_full <- fread("tables/pipolin_summary_new.tsv")

#prepare data
pipolins$family[is.na(pipolins$family)] <- "Unknown"
pipolins$order[is.na(pipolins$order)] <- "Unknown"
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

#2: prevalence plot (Figure 1)
summary <- summary %>% arrange(class,order,family,genus)

ggplot(summary[summary$pipolins>3,], aes(x = factor(genus,levels=genus), y = prevalence100, group=order)) +
  geom_bar(aes(fill = order), stat = "identity", col = "black") + 
  scale_y_continuous(expand = c(0,0), limits=c(0,100), breaks=c(20,40,60,80,100)) +
  geom_text(aes(label=pipolins),size = 3, position = position_dodge(width=0.9), hjust=-0.2, angle=90, color="grey20", fontface = "bold") + 
  facet_grid(.~fct_reorder(order, class), scale="free_x",space="free")+
  scale_fill_manual(values=color_pipolin_orders) + theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "italic")) + theme(plot.margin = margin(t = 0.5,  r = 0, b = 0, l = 0.4, unit = "cm"),legend.text=element_text(face="italic"),legend.position = "bottom") +
  theme(strip.background = element_blank(), strip.text.x.top =   element_blank(),  panel.grid.major = element_line(size = 0.5, linetype = 'solid',  colour = "white"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                                                                                                                                       colour = "white")) +
  xlab("Bacteria Genus")+ ylab("Pipolins prevalence %")+labs(fill = "Order")

#3: Figure 2 plots

#Figure 2A
#extract the data
num_pipolins <- as.data.frame(table(pipolins_full$Genome_ID))
names(num_pipolins) <- c("Genome_ID","num_pipolins")
pipolins_full <- merge(pipolins_full,num_pipolins)

pipolins_order <- as.data.frame(table(pipolins_full$class, pipolins_full$order,pipolins_full$num_pipolins))
pipolins_order <- pipolins_order[pipolins_order$Freq>20 & pipolins_order$Var1 !="",]
pipolins_order$Var2 <- fct_reorder(factor(pipolins_order$Var2),as.numeric(factor(pipolins_order$Var1)), .desc=TRUE)

n_o_pipolins <- ggplot(pipolins_order[pipolins_order$Freq>20 & pipolins_order$Var1 !="",], aes(fill=Var2, y=Freq, x=Var3)) +   geom_bar(position="stack", stat="identity") +  scale_fill_manual(values=color_pipolin_orders, name = "Order") +
  coord_flip() +  scale_y_log10(labels=c()) +
  geom_shadowtext(aes(label=ifelse(Freq>0,Freq,"")),size = 3,angle=90, position= position_stack(vjust = 0.5), color="white", fontface = "bold") +
  xlab("") + ylab("Number of Assemblies (log10)") + ggtitle("Pipolins per assembly")+
  theme_linedraw(base_size = 12) + theme(legend.text=element_text(face="italic"), axis.ticks.x=element_blank())

pipolb_order <- as.data.frame(table(pipolins_full$class,pipolins_full$order,pipolins_full$piPolB_num))
pipolb_order <- pipolb_order[pipolb_order$Freq>20 & pipolb_order$Var1 !="" & pipolb_order$Var2 !="",]

pipolb_order$Var2 <- fct_reorder(pipolb_order$Var2,as.numeric(pipolb_order$Var1), .desc=TRUE)
n_o_piPolBs <- ggplot(pipolb_order, aes(fill=Var2, y=Freq, x=Var3)) +   geom_bar(position="stack", stat="identity") +  scale_fill_manual(values=color_pipolin_orders, name = "Order") +
  coord_flip() +  scale_y_log10(labels=c()) +
  geom_shadowtext(aes(label=ifelse(Freq>0,Freq,"")),size = 3, angle=90,position = position_stack(vjust = 0.5), color="white", fontface = "bold") +
  xlab("") + ylab("Number of pipolins (log10)") + ggtitle("piPolBs per pipolin")+
  theme_linedraw(base_size = 12) + theme(legend.text=element_text(face="italic"), axis.ticks.x=element_blank())
#plot
n_o_pipolins / n_o_piPolBs +  plot_layout(ncol = 1, nrow = 2, guides = "collect")

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
  xlab("Number of DRs") + ylab("Number of pipolins (log10)") + 
  theme_linedraw(base_size = 12) + theme(legend.text=element_text(face="italic"),axis.ticks.x=element_blank())

#Figure 2C
pipolins_full$Number_atts2 <- ifelse(pipolins_full$Number_atts>2, "≥3",pipolins_full$Number_atts)
PieDonut(setNames(as.data.frame(table(pipolins_full$Number_atts2,pipolins_full$Assembly_gaps_pipolin_reconstruction)),c("DRs","Gaps","Freq")), aes(DRs, Gaps, count=Freq), showDonutName=TRUE, ratioByGroup = FALSE,r0=0,pieLabelSize = 5,
              donutLabelSize = 3) 

#Figure 2D and final table with delimited pipolins
delim <- read.table("data/Pipolin_delimited_new_int_site_length.tsv", sep="\t",header=TRUE)
delim <- delim[delim$delim_length!=1220,]
pipolins_delim <- pipolins_full
for (i in 1:nrow(pipolins_delim)){
  if (pipolins_delim$Pipolin_ID[i] %in% delim$Pipolin_id[delim$tdelim_type=="Alt_site"]) {
    pipolins_delim$pipolin_length[i] <- delim$delim_length[delim$Pipolin_id==pipolins_delim$Pipolin_ID[i]]
    pipolins_delim$Atts[i] <- "Alt_site"
  }
}
pipolins_delim$Atts[pipolins_delim$order=="Bacillales"] <- "Plasmid"
pipolins_delim <- pipolins_delim[pipolins_delim$Atts!="FALSE",]
#write table
write.table(pipolins_delim,"Supp_dataset1.tsv",quote=FALSE,row.names = FALSE,sep="\t")

tmp <- merge(pipolins_delim[,c(7,40,41)],O[,1:2])
tmp <- tmp[tmp$pipolins>20, ]

number <- as.data.frame(tmp %>% count(order) )
tmp$order <- fct_reorder(tmp$order,as.numeric(as.factor(tmp$class)),.desc=TRUE)

ggplot(tmp)+geom_density_ridges(aes(x=pipolin_length, y=order,fill=order),alpha=0.8) + scale_fill_manual(values=color_pipolin_orders)+scale_x_continuous(labels = scales::scientific,limits=c(0,100000))+ylab("")+xlab("Pipolin length")+
  geom_text(data=number,aes(y=order, x=99000,label=paste0("n= ",n)),color="grey40",vjust=-0.5,hjust=1,inherit.aes = FALSE)+theme_classic(base_size = 14) + theme(axis.text.y = element_text(face = "italic"),legend.position = "none")

