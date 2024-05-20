piPolB_clu_names = ""
with open("Cluster_info_simplified.txt", "r") as f:
    for line in f:
        if line.split("\t")[-3] == "piPolB":
            piPolB_clu_names +=line.split("\t")[0]+"\n"

with open("Clusters_with_piPolBs.txt", "w") as f:
    f.write(piPolB_clu_names)


