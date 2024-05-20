from Bio import SeqIO


cds_info = "MGE_id\tGene_ID\tProtein_length\n"

for record in SeqIO.parse("../MGE_pipolins_filtered.faa", "fasta"):
    mge_id = "_".join(record.id.split("_")[0:-1])
    cds_info += mge_id+"\t"+record.id+"\t"+str(len(record.seq)-1)+"\n" #-1 to remove final "*"

with open("MGE_Gene_info_filtered.tsv", "w") as f:
    f.write(cds_info)