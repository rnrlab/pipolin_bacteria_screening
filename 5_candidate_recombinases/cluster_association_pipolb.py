import os
import numpy as np
np.random.seed(1)

### Get clusters with >= 800% and 75% ali cov piPolBs
pipolb_clusters = []
with open("../pipolin_protein_clustering/Clusters_with_piPolBs.txt", "r") as f:
    for line in f:
        pipolb_clusters.append(line.replace("\n",""))

### Read cluster info for output and select clusters
cluster_index = {} #index to cluster
rec_clusters_index = {} #index to rec cluster
col_list = [] #pipolb clusters
index = 0
n_test = 0
with open("../pipolin_protein_clustering/Cluster_info_simplified.txt", "r") as f:
    for line in f:
        if "Cluster_" in line:
            fields = line.replace("\n", "").split("\t")
            cluster_index[index] = {"Name":fields[0],"Number":fields[1]}
            if fields[0] in pipolb_clusters: #First, I include all clusters with piPolBs. Then I filter by cluster size
                col_list.append(index)

            #check recombinase
            if line.split("\t")[-2]=="Yrec" or line.split("\t")[-2]=="Srec":
                rec_clusters_index[index] = {"Name":fields[0],"Number":fields[1]}
                n_test += 1
                
            elif line.split("\t")[-2]=="No PFAM" and "Int_" in line.split("\t")[-6]:
                rec_clusters_index[index] = {"Name":fields[0],"Number":fields[1]}
                n_test += 1

            index += 1

### Find integrase candidates
pa_M = np.loadtxt("../pipolin_phylogeny_correlations/piPolB_context_gene_presence_only_matrix.txt", dtype="int32")
n_col = 0
selected_col = []
selected_col_rec = []
selected_col_rec_list = "Cluster_n\tCluster_name\tIndex\tNear_piPolB\tdiff_freq\tEst_p-val\n"
graph_output = "Node_1\tNode_2\tdiff_freq\n" 
diff_freqs_matrix = "Ini\t"+"\t".join([str(n+1) for n in range(pa_M.shape[1])])+"\n"

for col in pa_M.T:
    if n_col in col_list: 
        if pa_M[:, n_col].sum() >= 40: #Only piPolB clusters where at least 40 piPolBs (>800 aa and 75% cov) are present
            print("Col index:", n_col)
            print("Num pipolb:", pa_M[:, n_col].sum())

            #2 matrix subsets: with cluster and without cluster
            paM_with = np.delete(pa_M, np.where(col== 0), axis=0)
            paM_withOut = np.delete(pa_M, np.where(col== 1), axis=0)

            
            #Calculate observed freq: Sum of columns in each subset and divide by row num
            paM_with_freqs = paM_with.sum(axis = 0)/paM_with.shape[0]
            paM_withOut_freqs = paM_withOut.sum(axis = 0)/paM_withOut.shape[0]


            #Diff of freqs to get candidates
            diff_freqs = paM_with_freqs-paM_withOut_freqs
            diff_freqs_matrix += str(n_col+1)+"\t"+"\t".join([str(val) for val in diff_freqs])+"\n"
            freq_clusters = np.where(diff_freqs > 0.05)
            tmp_freq_rec_clu = []


            #Output graph and list selected cols
            for index in freq_clusters[0]:
                #Keep only recombinase clusters
                if index in list(rec_clusters_index.keys()):
                    #Register recombinase
                    if index not in tmp_freq_rec_clu:
                        tmp_freq_rec_clu.append(index)

                        tmp_diff_freq_list = []
                        #Calculate empirical p-value
                        for i in range(0,1000):
                            
                            piPolB_col = pa_M[:,n_col]
                            random_col = np.random.randint(2, size=pa_M.shape[0])
                            pairM_piPolB_rec = np.column_stack((piPolB_col, random_col))

                            #pairM_piPolB_rec = pa_M[:,[n_col,index]]
                            #np.random.shuffle(pairM_piPolB_rec[:,1]) #change in place

                            pairM_with = np.delete(pairM_piPolB_rec, np.where(col== 0), axis=0)
                            pairM_withOut = np.delete(pairM_piPolB_rec, np.where(col== 1), axis=0)

                            pairM_with_freqs = pairM_with.sum(axis = 0)/pairM_with.shape[0]
                            pairM_withOut_freqs = pairM_withOut.sum(axis = 0)/pairM_withOut.shape[0]

                            pair_diff_freqs = pairM_with_freqs-pairM_withOut_freqs
                            tmp_diff_freq_list.append(pair_diff_freqs[1])
                        
                        emp_eval = np.sum(tmp_diff_freq_list >= diff_freqs[index])/len(tmp_diff_freq_list)

                        print("piPolB:", n_col, "rec:", index, "Obs diff freq:", diff_freqs[index], "emp_eval:", emp_eval)
                        
                        #Reporting in tsv
                        selected_col_rec_list += rec_clusters_index[index]["Number"]+"\t"+rec_clusters_index[index]["Name"]+"\t"+str(index+1)+"\t"+cluster_index[n_col]["Number"]+"\t"+str(round(diff_freqs[index],4))+"\t"+str(round(emp_eval,4))+"\n" #index+1 bc pyhton to R indez
                    
                    #Printing purposes
                    if index not in selected_col_rec:
                        selected_col_rec.append(index)
                """
                #For full heatmap (deprecated)
                if index not in selected_col:
                    selected_col.append(index)

                #For graph (deprecated)
                if cluster_index[n_col] != cluster_index[index]: #to avoid self-edges
                    graph_output += cluster_index[n_col]["Number"]+"\t"+cluster_index[index]["Number"]+"\t"+str(diff_freqs[index])+"\n"
                else:
                    graph_output += cluster_index[index]["Number"] + "\n"
                """

            print(n_col, cluster_index[n_col], tmp_freq_rec_clu)                 

    n_col += 1

print("All_YR-SR_clu: ",len(list(rec_clusters_index.keys())))
print("Selected_YR-SR_clu: ",len(selected_col_rec))
print(selected_col_rec)


with open("list_pipolb_assoc_rec_clusters_list.txt", "w") as f:
    f.write(selected_col_rec_list)


with open("diff_freq_matrix.txt", "w") as f:
    f.write(diff_freqs_matrix)

""" Deprecated output

with open("pipolb_associated_clusters_pipolins.txt", "w") as f:
    f.write(graph_output)

with open("list_pipolb_assoc_clusters.txt", "w") as f:
    f.write("Clusters\n"+'\n'.join([str(n+1) for n in selected_col])+"\n")

with open("list_pipolb_assoc_clusters_index_list.txt", "w") as f:
    f.write('pipolb_assoc_clu\n'+'\n'.join([str(n) for n in selected_col])+"\n")
"""
