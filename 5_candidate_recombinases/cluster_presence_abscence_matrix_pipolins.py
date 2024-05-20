import os

### Read pipolin cds cluster annotation info (includes order funct etc)
cluster_info = {}
with open("../pipolin_protein_clustering/Cluster_information.txt", "r") as f:
    for line in f:
        if "Cluster_" in line:
            fields = line.replace("\n", "").split("\t")
            cluster_info[fields[0]] = fields


### Read pipolin cds cluster components
cluster_components = {}
with open("../pipolin_protein_clustering/mmseqs_clustering/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db_clu_wide.tsv", "r") as f:
    for line in f:
        fields = line.replace("\n", "").split("\t")
        cluster_components[fields[0].split("|")[0]] = fields[1:]


### Read piPolB list
piPolB_list = []
with open("../pipolin_protein_clustering/hmmsearch_pipolb/pipolb_seq_filter800aa.faa", "r") as f:
    for line in f:
        if ">G" in line:
            piPolB_list.append(line.replace(">","").replace("\n",""))


### I need a table in format: piPolB_ID \t has its pipolin X gene from cluster A? \t has its pipolin Y gene from cluster B \t ...
piPolB_context_gene_presence = {}

for piPolB_id in piPolB_list:
    piPolB_context_gene_presence[piPolB_id] = []
    pipolin_id = "_".join(piPolB_id.split("_")[0:3])
    for mmseqscluster in list(cluster_info.keys()): #50 tmp
        if pipolin_id in "-".join(cluster_components[mmseqscluster]):
            piPolB_context_gene_presence[piPolB_id].append("1")
        else:
            piPolB_context_gene_presence[piPolB_id].append("0")



piPolB_context_gene_presence_output_table = "piPolB_repre\t"+"\t".join(list(v[1] for v in cluster_info.values()))+"\n"
piPolB_context_gene_presence_output_table_only_matrix = ""
for piPolB_id in piPolB_context_gene_presence:
    piPolB_context_gene_presence_output_table += piPolB_id+"\t"+"\t".join(piPolB_context_gene_presence[piPolB_id])+"\n"
    piPolB_context_gene_presence_output_table_only_matrix += "\t".join(piPolB_context_gene_presence[piPolB_id])+"\n"


with open("piPolB_context_gene_presence.txt", "w") as f:
    f.write(piPolB_context_gene_presence_output_table)

with open("piPolB_context_gene_presence_only_matrix.txt", "w") as f:
    f.write(piPolB_context_gene_presence_output_table_only_matrix)
