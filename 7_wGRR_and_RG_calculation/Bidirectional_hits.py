import sys

two_hits = 0
BH = ""
tmp_id = "NA"
tmp_alis = {"Dir1":[0,0,0,0,0,0,0],"Dir2":[0,0,0,0,0,0,0]} #mmseqs_all_vs_all/All_MGE_5-500kbp_RefSeq_allvsall_ali_nored_id35-eval4_cov50_hitID-ordered_sortb.m8
with open(sys.argv[1],"r") as f:
    for line in f: 
        if two_hits%1000000 == 0 and two_hits != 0:
            print(two_hits)
        hit_id = line.split("\t")[-1].replace("\n","")
        if hit_id != tmp_id:
            #1º) Close previous case
            if tmp_alis["Dir2"][2]!=0:
                #If True, save in appropiate format for R script
                prot_query_id = "_".join(tmp_alis["Dir1"][0].split("_")[:-1])
                prot_subject_id = "_".join(tmp_alis["Dir2"][0].split("_")[:-1])
                if prot_query_id != prot_subject_id: #Make sure we don't add cross hits within the same MGE
                    pair_id = prot_query_id+"_"+prot_subject_id
                    pair_id_rev = prot_subject_id+"_"+prot_query_id
                    BH += "\t".join(tmp_alis["Dir1"][0:6]+[pair_id])+"\n"
                    BH += "\t".join(tmp_alis["Dir2"][0:6]+[pair_id_rev])+"\n"
                    two_hits += 1
            #2º) Open new case
            tmp_id = hit_id
            tmp_alis["Dir1"] = line.replace("\n","").split("\t") #qseqid,sseqid,pident,qcov,scov,bitscore,hit_id
            tmp_alis["Dir2"] = [line.split("\t")[1], line.split("\t")[0], 0, 0, 0, 0,line.split("\t")[-1].replace("\n","")]  #inverse hit optimizable?

        else:
            if line.split("\t")[0]==tmp_alis["Dir1"][0]: #if same direction hit
                if int(line.split("\t")[5])> int(tmp_alis["Dir1"][5]): #see which has better bitscore
                    tmp_alis["Dir1"] = line.replace("\n","").split("\t")    #replace if that is the case
            elif line.split("\t")[0]==tmp_alis["Dir2"][0]:
                if int(line.split("\t")[5]) > int(tmp_alis["Dir2"][5]): #see which has better bitscore
                    tmp_alis["Dir2"] = line.replace("\n","").split("\t")     #replace if that is the case


#Add last case
if tmp_alis["Dir2"][2]!=0:
    prot_query_id = "_".join(tmp_alis["Dir1"][0].split("_")[:-1])
    prot_subject_id = "_".join(tmp_alis["Dir2"][0].split("_")[:-1])
    if prot_query_id != prot_subject_id: #Make sure we don't add cross hits within the same MGE
        pair_id = prot_query_id+"_"+prot_subject_id
        pair_id_rev = prot_subject_id+"_"+prot_query_id
        BH += "\t".join(tmp_alis["Dir1"][0:6]+[pair_id])+"\n"
        BH += "\t".join(tmp_alis["Dir2"][0:6]+[pair_id_rev])+"\n"
        two_hits += 1

with open(sys.argv[1]+"_BH", "w") as f:
    print("Writing in:", sys.argv[1]+"_BH")
    f.write(BH)
