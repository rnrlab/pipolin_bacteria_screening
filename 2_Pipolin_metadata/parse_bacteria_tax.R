#1: Load & prepare the data
library(data.table)
bacteria <- fread("data/datasets_bacteria.txt") #created for dataformat: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/
names(bacteria) <- c("Accession","Organism","TaxId")
pipolins_assemblies <- fread("data/genome_metadata_pipolins_updatedCont.txt")

#All in one table
bacteria$pipolin <- bacteria$Accession %in% pipolins_assemblies$AN

#fix missing TaxIds
bacteria$TaxId[bacteria$Accession=="GCA_000723365.1"] <- 572419
bacteria$TaxId[bacteria$Accession=="GCA_000626115.2"] <- 1412458
bacteria$TaxId[bacteria$Accession=="GCA_000831025.2"] <- 1412465
bacteria$TaxId[bacteria$Accession=="GCA_000624415.2"] <- 1412467
bacteria$TaxId[bacteria$Accession=="GCA_000624435.2"] <- 1412469
bacteria$TaxId[bacteria$Accession=="GCA_000625555.2"] <- 1412607
bacteria$TaxId[bacteria$Accession=="GCA_001278745.1"] <- 1280
bacteria$TaxId[bacteria$Accession=="GCA_000364805.1"] <- 1265601
bacteria$TaxId[bacteria$Accession=="GCA_000831045.2"] <- 1412464

#2: Get full taxonomy
library(taxonomizr)
#only once: prepareDatabase(getAccessions=FALSE)
taxonomy <- getTaxonomy(unique(bacteria$TaxId),'accessionTaxa.sql')
bac_tax <- as.data.frame(taxonomy)
bac_tax$TaxId <- as.numeric(row.names(bac_tax$TaxId))
bacteria_full <- merge(bacteria,bac_tax, all.x=TRUE)
write.table(bacteria_full,"tables/bacteria_tax.tsv", sep="\t", row.names=FALSE)

#3: pipolin stats by genus using xtabs
bacteria_full <- fread("tables/bacteria_tax.tsv", sep="\t")
stats_pipolin <- xtabs(pipolin~genus, bacteria_full)
stats_full <- xtabs(~genus, bacteria_full)
stats_pipolin <- as.data.frame(stats_pipolin)
names(stats_pipolin)[2] <- "pipolins"
stats_full <- as.data.frame(stats_full)
names(stats_full)[2] <- "genomes"
stats <- merge(stats_pipolin,stats_full)
stats$prevalence100 <- round(stats$pipolins *100 / stats$genomes,2) 
stats$ratio100 <- round(stats$pipolins * 100 / 11431 ,2) 
#now we merge everything
stats2 <- merge(stats,unique(bacteria_full[,c(10,9,8,7,6)]), all.y=FALSE)
write.table(stats2,"tables/pipolin_bacteria_stats.tsv", sep="\t", row.names=FALSE)

#4: full pipolins table
bacteria_pipolins <- bacteria_full[bacteria_full$pipolin==TRUE]
structure <- fread("data/pipolin_summary_postprocessed.txt")
names(structure)[32] <- "Accession"
pipolin_full <- merge(structure,bacteria_full)
fwrite(pipolin_full,"tables/pipolin_summary_postprocessed_new.tsv", sep="\t")
