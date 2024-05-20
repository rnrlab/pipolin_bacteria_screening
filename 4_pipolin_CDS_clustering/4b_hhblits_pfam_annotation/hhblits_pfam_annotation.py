import os
import subprocess
import math
from Bio import SeqIO

### Switches
filter_task = False
mafft_task = False
hhblits_task = True
results_task = True

### Filter clusters with cd hit 
if filter_task:
    subprocess.run("mkdir -p filtered_clusters",shell=True)
    subprocess.run("rm -r filtered_clusters/*",shell=True)

    for cluster in os.listdir("../mmseqs_clustering/clusters"):
        seqs_list = [record for record in SeqIO.parse("../mmseqs_clustering/clusters/"+cluster, "fasta")]
        if len(seqs_list) >= 5:
            subprocess.run("cdhit -i ../mmseqs_clustering/clusters/"+cluster+" -o filtered_clusters/"+cluster.replace(".fa", "_90id.fa")+" -c 0.9 -T 30 -d 0 -M 0", shell=True)

### align with mafft
if mafft_task:
    subprocess.run("mkdir -p filtered_clusters_mafft",shell=True)
    subprocess.run("rm -r filtered_clusters_mafft/*",shell=True)

    subprocess.run("ls filtered_clusters/* | grep -v '.clstr' | parallel --verbose 'mafft --auto {} > filtered_clusters_mafft/{/.}_mafft.aln'", shell=True)

### hhblits task

if hhblits_task:
    subprocess.run("mkdir -p filtered_clusters_mafft_hhblitsPfam",shell=True)
    subprocess.run("rm -r filtered_clusters_mafft_hhblitsPfam/*",shell=True)
    
    subprocess.run("ls filtered_clusters_mafft/* | parallel --verbose 'hhblits -i {} -o filtered_clusters_mafft_hhblitsPfam/{/.}.hhr -d $HOME/DB_Victor/Pfam_35/pfam -n 5 -M first' ", shell=True)

#Parse results
if results_task:
    cluster_hhblits_pfam = {}
    for hhr in os.listdir("filtered_clusters_mafft_hhblitsPfam"):
        with open("filtered_clusters_mafft_hhblitsPfam/"+hhr, "r") as f:
            n_hits = 3
            E_val_limit = 0.001
            cluster_N = hhr.split("_")[1]
            cluster_name = "_".join(hhr.split("_")[2:6]).replace(".hhr","")
            cluster_key = cluster_N+"\t"+cluster_name
            cluster_hhblits_pfam[cluster_key] = []
            hit = False
            for line in f:
                if n_hits != 0:
                    if ">" in line:
                        hit = True
                        desc = line.replace("\n","").replace(">", "")
                    else:
                        if hit:
                            hit = False
                            n_hits -= 1
                            E_val = eval(line.split()[1].split("=")[1])
                            if E_val < E_val_limit:
                                cluster_hhblits_pfam[cluster_key].append(desc+" ; "+str(E_val))

    cluster_hhblits_pfam_table = "#N\tCluster_name\tAnnotations\n"
    for cluster in cluster_hhblits_pfam:
        cluster_hhblits_pfam_table += cluster.replace("_90id_mafft", "")+"\t"+str(cluster_hhblits_pfam[cluster])+"\n"

    with open("Cluster_hhblits_pfam.txt", "w") as f:
        f.write(cluster_hhblits_pfam_table)




                    
