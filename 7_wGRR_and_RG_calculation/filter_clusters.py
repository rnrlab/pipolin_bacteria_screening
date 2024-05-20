import sys
from Bio import SeqIO
clusters = {}
#First run, read cluster names (aka representative)
with open(sys.argv[1], "r") as f:
    for line in f:
        clusters[line.split()[0]] = []

#Second run, attatch components 
with open(sys.argv[1], "r") as f:
    for line in f:
        clusters[line.split()[0]].append(line.replace("\n","").split()[1])


cluster_pipolin_ratio = "Cluster\tn_pipolin\tn_mge\tratio\n"
clusters_pipolins_list = ""
mge_list = set() #unique list of mge candidates
for c in clusters:

    pipolin = False
    n_pipolin = 0
    n_mge = 0
    for cds_id in clusters[c]:
        n_mge += 1
        if "G_" == cds_id[:2] and "v0_" in cds_id:
            n_pipolin += 1
    if n_pipolin >= 3:
        cluster_pipolin_ratio += c+"\t"+str(n_pipolin)+"\t"+str(n_mge)+"\t"+str(n_pipolin/n_mge)+"\n"
    if n_pipolin >= 3 and n_pipolin/n_mge >= 0.03 and (n_mge-n_pipolin) < 2000:
        for cds_id in clusters[c]:
            mge_id = "_".join(cds_id.split("_")[:-1])
            if "G_" != cds_id[:2] and "v0_" not in cds_id:
                mge_list.add(mge_id)
            clusters_pipolins_list += c+"\t"+cds_id+"\n"

print(len(mge_list))
with open("Selected_MGEs.txt", "w") as f:
    for mge in mge_list:
        f.write(mge+"\n")

with open("linclust_all_mge/Clusters_pipolins.txt", "w") as f:
    f.write(clusters_pipolins_list)

with open("linclust_all_mge/Clusters_pipolin_ratio.txt", "w") as f:
    f.write(cluster_pipolin_ratio)


records_filtered = []
#Keep CDS from the selected MGE
for record in SeqIO.parse(sys.argv[2], "fasta"):
    mge_id = "_".join(record.id.split("_")[:-1])
    if mge_id in mge_list:
        records_filtered.append(record)
    else:
        if "G_" == mge_id[:2] and "v0" in mge_id:
            records_filtered.append(record)

SeqIO.write(records_filtered, sys.argv[2].replace(".faa","_filtered.faa"), "fasta")
