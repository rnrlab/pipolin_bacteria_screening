
import os
import subprocess
from Bio import SeqIO

cds_info = "VFID\tProtein_length\tVF_gene\tVF_class\n"

for record in SeqIO.parse("VFDB_setB_pro.fas", "fasta"):
    print(record.description)
    VF_gene = record.description.split(") - ")[0].split(" [")[-1]+")"
    VF_class = record.description.split(") - ")[1].split("] ")[0]
    cds_info += record.id+"\t"+str(len(record.seq)-1)+"\t"+VF_gene+"\t"+VF_class+"\n" #-1 to remove final "*"

with open("VF_id_info_len.tsv", "w") as f:
    f.write(cds_info)

subprocess.run("mkdir VFDB_mmseqs", shell = True)
subprocess.run("mmseqs createdb VFDB_setB_pro.fas VFDB_mmseqs/VFDB_setB_pro_DB", shell = True)
subprocess.run("mmseqs search ../../linclust_all_mge/MGE_pipolins_db VFDB_mmseqs/VFDB_setB_pro_DB VFDB_mmseqs/MGE_pipolins_filtered_VFDB_ali VFDB_mmseqs_tmp") #-s 7.5
subprocess.run("mmseqs convertalis ../../linclust_all_mge/MGE_pipolins_db VFDB_mmseqs/VFDB_setB_pro_DB VFDB_mmseqs/MGE_pipolins_filtered_VFDB_ali VFDB_mmseqs/MGE_pipolins_filtered_VFDB_ali.m8")
