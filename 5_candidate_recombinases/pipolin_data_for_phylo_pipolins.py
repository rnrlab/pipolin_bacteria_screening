### Read piPolB list
piPolB_list_phylo = []
with open("../piPolB_analysis/pipolb_seq_filter800aa.faa", "r") as f:
    for line in f:
        if ">G" in line:
            piPolB_list_phylo.append(line.replace(">","").replace("\n",""))


### Read pipolin info
pipolin_info = {}
with open("../pipolin_summary_new.tsv", "r") as f :
    for line in f:
        fields = line.replace("\n", "").split("\t")
        pipolin_info[fields[1]] = fields

### Compute table with requested info
table_info_pipolins_phylo = "piPolB_repre\tpipolin\torder\tatts\tintegration_site\tpiPolB_len_sum\n"
for piPolB_repre in piPolB_list_phylo:

    pipolin_id = "_".join(piPolB_repre.split("_")[0:3])
    table_info_pipolins_phylo += piPolB_repre+"\t"+pipolin_id+"\t"+pipolin_info[pipolin_id][-4]+"\t"+pipolin_info[pipolin_id][18]+"\t"+pipolin_info[pipolin_id][24].replace("hypothetical protein, ","").replace("assembly gap, ","")+"\t"+pipolin_info[pipolin_id][9]+"\n"


with open("table_info_pipolins_phylo.txt", "w") as f:
    f.write(table_info_pipolins_phylo)

