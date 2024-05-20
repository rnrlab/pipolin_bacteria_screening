import time
import os
from Bio import SeqIO
start_time = time.time()

### WARNING: SLOW PROGRAM

#Obtain list of all protein ids
protein_ids_list = []
with open("bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG-db_70_filtered_ids.txt", "r") as f:
    for line in f:
        protein_ids_list.append(line.replace("\n", ""))

#Obtain cluster components
cluster_components = {}
with open("mmseqs_clustering/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db_clu_wide.tsv", "r") as f:
    for line in f:
        fields = line.replace("\n", "").split("\t")
        cluster_components[fields[0]] = fields[1:]

#Obtain cluster num (order), name and mean orf length
cluster_list_num_name = {}
cluster_mean_cds_length = {}
for file in os.listdir("mmseqs_clustering/clusters"):
    cluster_num = int(file.split("_")[1])
    cluster_name = file.replace(".fa","")
    cluster_list_num_name[cluster_name] = cluster_num
    CDSlength_list = []
    for CDS in SeqIO.parse("mmseqs_clustering/clusters/"+file, "fasta"):
        CDSlength_list.append(len(CDS))
    cluster_mean_cds_length[cluster_num] = round(sum(CDSlength_list)/len(CDSlength_list))

cluster_list_num_name_ordered = sorted(cluster_list_num_name, key=lambda k: cluster_list_num_name[k], reverse=False)


#Apply order to cluster components dictionary

cluster_size_filtered_ordered = {}
for key in cluster_list_num_name_ordered:
    if "_G_" in key:
        key_clean = "G_" + key.split("_G_")[-1]
        cluster_size_filtered_ordered[key_clean] = cluster_components[key_clean]
    elif "_mobileOG_":
        mobile_OG_num = key.split("_")[-1]
        for mobile_OG_key in list(cluster_components.keys()):
            if mobile_OG_num in mobile_OG_key:
                cluster_size_filtered_ordered[mobile_OG_key] = cluster_components[mobile_OG_key]

#Extract mobileOG informaton
mobile_OG_info = {}
for cluster in cluster_size_filtered_ordered:

    if "mobileOG" in cluster or "mobileOG" in str(cluster_size_filtered_ordered[cluster]):

        for member in cluster_size_filtered_ordered[cluster]:

            if "|" in member:
                mobile_OG_desc = member.split("|")[1] +"|"+member.split("|")[3]+"|"+member.split("|")[4]

                if cluster not in list(mobile_OG_info.keys()):
                    mobile_OG_info[cluster] = {}
                    mobile_OG_info[cluster][mobile_OG_desc] = 1
                else:
                    if mobile_OG_desc not in list(mobile_OG_info[cluster].keys()):
                        mobile_OG_info[cluster][mobile_OG_desc] = 1
                    else:
                        mobile_OG_info[cluster][mobile_OG_desc] += 1
    else:
        mobile_OG_info[cluster] = {"-":0}


#Load piPolB hmmsearch results
piPolB_hmmresults = {}
with open("hmmsearch_pipolb/bacteria_dec2022_pipolin_cds_partial_filtered_hmmsearch_pipolbs.domtbl", "r") as f:
    for line in f:
        if line[0] != "#":
            fields = line.split()

            if fields[0] not in list(piPolB_hmmresults.keys()):
                piPolB_hmmresults[fields[0]] = [[fields[17]], [fields[18]], int(fields[18])-int(fields[17])]
            else:
                piPolB_hmmresults[fields[0]][0].append(fields[17])
                piPolB_hmmresults[fields[0]][1].append(fields[18])
                piPolB_hmmresults[fields[0]][2] += (int(fields[18])-int(fields[17]))


piPolB_list = []
for piPolB in piPolB_hmmresults:
        if piPolB_hmmresults[piPolB][2] > 400:
            piPolB_list.append(piPolB)


#Load MOBscan results
relaxase_MOBscan = {}
with open("hmmsearch_MOBrelaxase/results_60.csv", "r") as f:
    for line in f:
        relaxase_MOBscan[line.split()[0]] = "Relaxase_"+line.split()[1]


#Load extra rec results
extra_rec = {}
with open("hmmsearch_extra_recombinase/hmmsearch_extra_rec_results.txt", "r") as f:
    for line in f:
        id = line.split()[0]
        rec = line.split()[1]
        extra_rec[id] = rec



#Load EggNog mapper info #{Sequence}=[EggNog_Class, Description, Preferred_name, PFAM]

EggNog = {}
for id in protein_ids_list:
    EggNog[id] = ["-","-","-","-"]
with open("reannotation_eggnog/bacteria_2022_pipolins_proteins_partial_EggNog.tsv", "r") as f:
    for line in f:
        if line[0] != "#":
            fields = line.replace("\n", "").split("\t")
            member = fields[0]
            EggNog[member] = [fields[6], fields[7], fields[8], fields[20]]


# HHblits information
hhblits_pfam_cluster = {}
hhblits_pfam_cluster_tophit = {}
with open("cluster_annotation_hhblits_pfam/Cluster_hhblits_pfam.txt", "r") as f:
    for line in f:
        if "#" not in line:
            cluster = line.split()[1]
            hhblits_pfam_cluster[cluster] = line.split("\t")[-1].replace("\n", "")
            if line.split("\t")[-1].replace("\n", "") != "[]":
                hhblits_pfam_cluster_tophit[cluster] = hhblits_pfam_cluster[str(cluster).split("|")[0]].split(" ; ")[1]
            else:
                hhblits_pfam_cluster_tophit[cluster] = ""


### ANNOTATION ###
cluster_annotations = {}
for i in range(len(cluster_size_filtered_ordered)):
    cluster = list(cluster_size_filtered_ordered.keys())[i]
    cluster_members = cluster_size_filtered_ordered[cluster]
    mean_cds_length = cluster_mean_cds_length[i+1]

    size_pipolins = 0
    ### Calculate size (only pipolins)
    for member in cluster_members:
        if "v0_" in member:
            size_pipolins += 1

    if size_pipolins >= 5: 
        print("#######################")
        print("Gathering anntoations for cluster", i+1)

        piPolB_hmm = {"piPolB": 0, "Other": 0}
        piPolB_num = 0
        MOB_hmm = {}
        rec_hmm = {}
        clu_EggNog_Class = {}
        clu_desc = {}
        clu_pref_name = {}
        clu_PFAM = {}

        for member in cluster_members:
            #Check pipolb
            if member in piPolB_list:
                piPolB_hmm["piPolB"] += 1
                piPolB_num += 1
            else:
                piPolB_hmm["Other"] += 1
            
            #Check MOBrel
            if member in list(relaxase_MOBscan.keys()):
                tmp_MOBrel = relaxase_MOBscan[member]
                if tmp_MOBrel not in list(MOB_hmm.keys()):
                    MOB_hmm[tmp_MOBrel] = 1
                else:
                    MOB_hmm[tmp_MOBrel] += 1

            #Check recs
                    
            if member in list(extra_rec.keys()):
                tmp_rec = extra_rec[member]
                if tmp_rec not in list(rec_hmm.keys()):
                    rec_hmm[tmp_rec] = 1
                else:
                    rec_hmm[tmp_rec] += 1
            

            

            #Add Eggnog
            #{Sequence}=[seed_ortholog, Description, Preferred_name, PFAM]
            if EggNog[member][0] not in list(clu_EggNog_Class.keys()):
                clu_EggNog_Class[EggNog[member][0]] = 1
            else:
                clu_EggNog_Class[EggNog[member][0]] += 1

            if EggNog[member][1] not in list(clu_desc.keys()):
                clu_desc[EggNog[member][1]] = 1
            else:
                clu_desc[EggNog[member][1]] += 1


            if EggNog[member][3] not in list(clu_PFAM.keys()):
                clu_PFAM[EggNog[member][3]] = 1
            else:
                clu_PFAM[EggNog[member][3]] += 1

        #Add mobile OG
        mobile_OG_stats = "-"
        if cluster in list(mobile_OG_info.keys()):
            mobile_OG_stats = mobile_OG_info[cluster]
        
        print(rec_hmm)
    
        cluster_annotations[str(cluster).split("|")[0]] = ["Cluster_"+str(i+1), str(len(cluster_members)), str(size_pipolins), str(mean_cds_length), str(clu_EggNog_Class), str(clu_desc), str(clu_PFAM), 
                                                            str(piPolB_hmm), str(piPolB_num), str(MOB_hmm), str(rec_hmm), str(mobile_OG_stats),
                                                            hhblits_pfam_cluster_tophit[str(cluster).split("|")[0]], hhblits_pfam_cluster[str(cluster).split("|")[0]]]
        
        print("End annotating "+str(i+1)+"--- %s seconds ---" % (time.time() - start_time))
        print("#######################")

### Writing file
clusters_info_output = "Representative\tCluster\tTrue_Size\tPipolin_degree\tMean_CDS_Length\tSeed_ortholog\tDescription\tPFAM\tpipolb_hmm\tpipolb_num_clu\tMOB_scan\ttop_hit_rec_hmm\tmobileOG\thhblits_pfam35_tophit\thhblits_pfam35_3hit\n"
for cluster in cluster_annotations:
    clusters_info_output += cluster + "\t" + "\t".join(cluster_annotations[cluster]) + "\n"

with open("Cluster_information.txt", "w") as f:
    f.write(clusters_info_output) 


print("--- %s seconds ---" % (time.time() - start_time))