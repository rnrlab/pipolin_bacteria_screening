library(ggplot2)
library(tibble)
library(scales)
library(ggrepel)
library(forcats)
library(RColorBrewer)
library(randomcoloR)


### Pipolb stats ###

#1. read data
pipolb_data <- read.table("pipolb_information.txt", header = TRUE, sep="\t", quote = "\"", 
                          comment.char = "$", na.strings = "N@")

#2. basic stats
max(pipolb_data$length)
min(pipolb_data$length)

ggplot(pipolb_data, aes(x=length)) + 
  geom_histogram(binwidth = 10) + 
  geom_vline(xintercept = 700)

sum(pipolb_data$length > 700) 
sum(pipolb_data$length <= 700) 


#3. add genome info data 
genome_data <- read.table("genome_metadata_pipolins_updatedCont.txt", header = TRUE, sep="\t", quote = "\"", comment.char = "$", na.strings = "N@")
rev_genus <- c("Escherichia", "Vibrio", "Enterobacter", "Pseudosulfitobacter", "Aeromonas", "Staphylococcus", "Salmonella", "Citrobacter", "Limosilactobacillus", "Corynebacterium")
genome_data_rev_genus <- genome_data[genome_data$genus_cont %in% rev_genus,] #Relevant
genome_data_rev_genus$Genus_filtered <- genome_data_rev_genus$genus_cont
genome_data_other <- genome_data[!genome_data$genus_cont %in% rev_genus,] #Other
genome_data_other$Genus_filtered <- rep("Other",length(genome_data_other$genus_cont))
genome_data_genus_filtered <- rbind(genome_data_rev_genus, genome_data_other)
genome_data_genus_filtered$genome <- genome_data_genus_filtered$G

pipolb_data_g <- merge(genome_data_genus_filtered, pipolb_data, by = "genome")
ggplot(data=pipolb_data_g, aes(x=length, fill=Genus_filtered)) + geom_histogram(binwidth=10)


length(pipolb_data_g[pipolb_data_g$genus_cont=="Pseudosulfitobacter",]$length)
