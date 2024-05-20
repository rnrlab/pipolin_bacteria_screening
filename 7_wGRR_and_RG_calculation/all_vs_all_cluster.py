import os
import subprocess
from Bio import SeqIO




cluster_size = {}
for cluster_fa in os.listdir("mmseqs_all_vs_all/clusters"):
    n = 0
    for record in SeqIO.parse("mmseqs_all_vs_all/clusters/"+cluster_fa, "fasta"):
        n += 1
    cluster_size[cluster_fa] = n

print("Clusters: ", len(cluster_size))


subprocess.run("mkdir mmseqs_all_vs_all/clusters_all_vs_all", shell=True)
for cluster_fa in cluster_size:
    #run all_vs_all
    if cluster_size[cluster_fa] > 1:
        cluster = cluster_fa.replace(".fa","")
        subprocess.run("mkdir mmseqs_all_vs_all/clusters_all_vs_all/"+cluster, shell=True)
        subprocess.run("mmseqs createdb mmseqs_all_vs_all/clusters/"+cluster_fa+" mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db", shell=True)
        subprocess.run("mmseqs prefilter mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db_pre -c 0.5 -s 7.5 --max-seqs "+str(cluster_size[cluster_fa]),shell=True)
        subprocess.run("mmseqs align mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db_pre  mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_allvsall_ali -a -c 0.5 -e 0.0001 --min-seq-id 0.35", shell=True)
        subprocess.run("mmseqs convertalis mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_db mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_allvsall_ali  mmseqs_all_vs_all/clusters_all_vs_all/"+cluster+"/"+cluster+"_allvsall_ali.m8 --format-output query,target,fident,qcov,tcov,bits", shell=True)



