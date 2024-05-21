from Bio import SeqIO
import os
import glob
import subprocess

n_drs = 0
n_alt_site = 0
n_nodelim = 0
n_no_full_pipolb = 0

new_length_table = "Pipolin_id\tdelim_length\ttdelim_type\n"

def write_proteins(gb_record, left_boundary, right_boundary, output_name):
    proteins_output = "" #This will store the proteins for output
    n_cds_final = 0
    for feature in gb_record.features:
            if feature.type == "CDS" and feature.location.start.position+10 > left_boundary and feature.location.end.position < right_boundary+10:  #+10 to avoid including boundary genes
                feat_header = " # ".join([">"+feature.qualifiers["locus_tag"][0], str(feature.location.start.position), str(feature.location.end.position), str(feature.location.strand), "ID="+feature.qualifiers["locus_tag"][0]])
                proteins_output += feat_header+"\n"+feature.qualifiers["translation"][0]+"\n"
                n_cds_final += 1
    if n_cds_final > 8:
        with open(output_name, "w") as tmp_fo:
            tmp_fo.write(proteins_output)

def delimit_pipolin(piPolB_coord, A_site, B_site, gb_record, tag):
    limit_check = True
    for coord in piPolB_coord:
        if coord not in range(min(A_site,B_site), max(A_site,B_site)):
            limit_check = False
        
    if limit_check:
        write_proteins(gb_record, min(A_site,B_site), max(A_site,B_site), "pipolins_5k_delimited/"+pipolin_id.replace("gbk",tag+".faa"))

    else: #In rare cases pipolin where boundaries have changed
        write_proteins(gb_record, min(piPolB_coord)-10000, max(piPolB_coord)+10000, "pipolins_5k_delimited/"+pipolin_id.replace("gbk","No_delim.faa"))
                    

if "pipolins_5k_delimited" not in os.listdir():
    subprocess.run("mkdir pipolins_5k_delimited", shell=True)
else:
    subprocess.run("rm -r pipolins_5k_delimited/*", shell=True)



for pipolin in os.listdir("../pipolin_representations/list_all_pipolins_5k/gbk_files"):
    pipolin_id = "G_"+pipolin.split("_G_")[-1]
    for record in SeqIO.parse("../pipolin_representations/list_all_pipolins_5k/gbk_files/"+pipolin,"genbank"):
        #1) Set varaibles for delimiting

        #piPolB len
        piPolB_len = 0

        #Vibrio
        clu_55_coords = "NA"
        clu_65_coords = "NA"

        #Aeromonas
        clu_2000_coords = "NA"
        clu_2353_coords = "NA"

        #Alphaprot
        clu_567_coords = "NA"
        clu_428_coords = "NA"
        clu_2075_coords = "NA"
        clu_1230_coords = "NA"
        clu_1596_coords = "NA"

        #Actinobacteria
        clu_706_coords = "NA"
        clu_882_coords = "NA"
        clu_1892_coords = "NA"
        clu_1825_coords = "NA"

        #Clostridia
        clu_1447_coords = "NA"
        clu_1138_coords = "NA"
        clu_3691_coords = "NA"
        clu_364_coords = "NA"
        clu_391_coords = "NA"

        #Lactobacillaceae
        clu_2265_coords = "NA"
        clu_3823_coords = "NA"


        no_gap = True
        piPolB_clu = "NA"
        piPolB_coord = []
        DR_coords = []
    
        for feature in record.features:
            if feature.type == "CDS":
                if "Cluster" in list(feature.qualifiers.keys()):
                    #Vibrio
                    if feature.qualifiers["Cluster"][0] == "Cluster_55":
                        clu_55_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_65":
                        clu_65_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    
                    #Aeromonas
                    if feature.qualifiers["Cluster"][0] == "Cluster_2000":
                        clu_2000_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_2353":
                        clu_2353coords = round((int(feature.location.end)+int(feature.location.start))/2)                                

                    #Alpharoteobacteria
                    if feature.qualifiers["Cluster"][0] == "Cluster_567":
                        clu_567_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_428":
                        clu_428_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_2075":
                        clu_2075_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_1230":
                        clu_1230_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_1596":
                        clu_1596_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    
                    #Actinobacteria
                    if feature.qualifiers["Cluster"][0] == "Cluster_706":
                        clu_706_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_882":
                        clu_882_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_1892":
                        clu_1892_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_1825":
                        clu_1825_coords = round((int(feature.location.end)+int(feature.location.start))/2)                               

                    #Clostridia
                    if feature.qualifiers["Cluster"][0] == "Cluster_1447":
                        clu_1447_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_1138":
                        clu_1138_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_3691":
                        clu_3691_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_364":
                        clu_364_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_391":
                        clu_391_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    
                    #Lactobacillaceae
                    if feature.qualifiers["Cluster"][0] == "Cluster_2265":
                        clu_2265_coords = round((int(feature.location.end)+int(feature.location.start))/2)
                    if feature.qualifiers["Cluster"][0] == "Cluster_3823":
                        clu_3823_coords = round((int(feature.location.end)+int(feature.location.start))/2)   

                if "Custom_category" in list(feature.qualifiers.keys()):
                    if feature.qualifiers["Custom_category"][0] == "piPolB":
                        piPolB_coord += [int(feature.location.start), int(feature.location.end)]
                        piPolB_len += len(feature.qualifiers["translation"][0])

            if feature.type == "repeat_region" and len(feature) >= 20:
                DR_coords.append(int(feature.location.start))

            if feature.type == "assembly_gap":
                    no_gap = False
            

        #Get DR boundaries (closest to piPolB)
        left_drs = []
        rigth_drs = []
        if len(piPolB_coord) != 0 and piPolB_len >= 800:
            for DR in DR_coords:
                if DR < min(piPolB_coord):
                    left_drs.append(DR)
                elif DR > max(piPolB_coord):
                    rigth_drs.append(DR)
                #else (DRs between pipolbs? we don't expect that)
        else: #piPolB not annotated (partial fragments, etc.)
            n_no_full_pipolb += 1
            #Do not include pipolin in this case
            #match_prodigal_out = glob.glob("../pipolin_protein_clustering/pipolin_prodigal_cds_faa/"+pipolin_id.replace(".gbk","")+"*")[0]
            #subprocess.run("cp "+match_prodigal_out+" pipolins_5k_delimited/"+pipolin_id.replace("gbk","faa"), shell=True)
            continue


        #First check DRs presence
        if len(left_drs)!=0 and len(rigth_drs)!=0:
            delimit_pipolin(piPolB_coord, max(left_drs), min(rigth_drs), record, "DR")
            new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(min(rigth_drs)-max(left_drs))+"\tEP_DRs\n"
            n_drs += 1

        
        #If not, we search for alternative/new boundaries
        #verify vibrio
        elif clu_55_coords != "NA" and clu_65_coords != "NA":
            n_alt_site += 1  
            delimit_pipolin(piPolB_coord, clu_55_coords, clu_65_coords, record, "Alt")
            new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_55_coords, clu_65_coords)-min(clu_55_coords, clu_65_coords))+"\tAlt_site\n"
                    
        #verify aeromonas
        elif clu_2000_coords != "NA" and clu_2353_coords != "NA":   
            delimit_pipolin(piPolB_coord, clu_2000_coords, clu_2353_coords, record, "Alt")
            new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_2000_coords, clu_2353_coords)-min(clu_2000_coords, clu_2353_coords))+"\tAlt_site\n"

        #Verify alphaprot
        elif clu_428_coords != "NA" or clu_1230_coords != "NA" or clu_1596_coords != "NA" or  clu_2075_coords != "NA": 
            
            if clu_1596_coords != "NA" and clu_2075_coords != "NA":
                n_alt_site += 1                              
                delimit_pipolin(piPolB_coord, clu_1596_coords, clu_2075_coords, record, "Alt")
                new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_1596_coords, clu_2075_coords)-min(clu_1596_coords, clu_2075_coords))+"\tAlt_site\n"

            elif clu_428_coords != "NA" and clu_1230_coords != "NA":
                n_alt_site += 1                       
                delimit_pipolin(piPolB_coord, clu_428_coords, clu_1230_coords, record, "Alt")
                new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_428_coords, clu_1230_coords)-min(clu_428_coords, clu_1230_coords))+"\tAlt_site\n"

            elif clu_428_coords != "NA" and clu_1596_coords != "NA":
                n_alt_site += 1                       
                delimit_pipolin(piPolB_coord, clu_428_coords, clu_1596_coords, record, "Alt")         
                new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_428_coords, clu_1596_coords)-min(clu_428_coords, clu_1596_coords))+"\tAlt_site\n" 
            elif clu_567_coords != "NA" and clu_428_coords != "NA":
                n_alt_site += 1                       
                delimit_pipolin(piPolB_coord, clu_428_coords, clu_567_coords, record, "Alt")   
                new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_428_coords, clu_567_coords)-min(clu_428_coords, clu_567_coords))+"\tAlt_site\n"
            else:
                n_nodelim += 1
                delimit_pipolin(piPolB_coord, min(piPolB_coord)-10000, max(piPolB_coord)+10000, record, "No_delim")

        #Verify actinobacteria
        elif clu_1892_coords != "NA" and clu_1825_coords != "NA":
            n_alt_site += 1                       
            delimit_pipolin(piPolB_coord, clu_1892_coords, clu_1825_coords, record,"Alt") 
            new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_1892_coords, clu_1825_coords)-min(clu_1892_coords, clu_1825_coords))+"\tAlt_site\n"
        elif clu_882_coords != "NA" and clu_706_coords != "NA":
            n_alt_site += 1                       
            delimit_pipolin(piPolB_coord, clu_882_coords, clu_706_coords, record, "Alt")
            new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_882_coords, clu_706_coords)-min(clu_882_coords, clu_706_coords))+"\tAlt_site\n"

        #Verify clostridia
        elif clu_1447_coords != "NA" or clu_1138_coords != "NA" or clu_3691_coords != "NA":
            if clu_1447_coords != "NA":
                if clu_364_coords != "NA":
                    n_alt_site += 1                       
                    delimit_pipolin(piPolB_coord, clu_1447_coords, clu_364_coords, record, "Alt")
                    new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_1447_coords, clu_364_coords)-min(clu_1447_coords, clu_364_coords))+"\tAlt_site\n"
                elif clu_391_coords != "NA":
                    n_alt_site += 1
                    delimit_pipolin(piPolB_coord, clu_1447_coords, clu_391_coords, record, "Alt")
                    new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_1447_coords, clu_391_coords)-min(clu_1447_coords, clu_391_coords))+"\tAlt_site\n"
                else:
                    n_nodelim += 1
                    delimit_pipolin(piPolB_coord, min(piPolB_coord)-10000, max(piPolB_coord)+10000, record, "No_delim")
            elif clu_3691_coords != "NA":
                if clu_364_coords != "NA":
                    n_alt_site += 1                       
                    delimit_pipolin(piPolB_coord, clu_3691_coords, clu_364_coords, record, "Alt")
                    new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_3691_coords, clu_364_coords)-min(clu_3691_coords, clu_364_coords))+"\tAlt_site\n"
                elif clu_391_coords != "NA":
                    n_alt_site += 1
                    delimit_pipolin(piPolB_coord, clu_3691_coords, clu_391_coords, record, "Alt")
                    new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_3691_coords, clu_391_coords)-min(clu_3691_coords, clu_391_coords))+"\tAlt_site\n"
                else:
                    n_nodelim += 1
                    delimit_pipolin(piPolB_coord, min(piPolB_coord)-10000, max(piPolB_coord)+10000, record, "No_delim")
            elif clu_1138_coords != "NA":
                if clu_364_coords != "NA":
                    n_alt_site += 1                    
                    delimit_pipolin(piPolB_coord, clu_1138_coords, clu_364_coords, record, "Alt")
                    new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_1138_coords, clu_364_coords)-min(clu_1138_coords, clu_364_coords))+"\tAlt_site\n"
                elif clu_391_coords != "NA":
                    n_alt_site += 1
                    delimit_pipolin(piPolB_coord, clu_1138_coords, clu_391_coords, record, "Alt")
                    new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_1138_coords, clu_391_coords)-min(clu_1138_coords, clu_391_coords))+"\tAlt_site\n"
                else:
                    n_nodelim += 1
                    delimit_pipolin(piPolB_coord, min(piPolB_coord)-10000, max(piPolB_coord)+10000, record, "No_delim")
            else:
                n_nodelim += 1
                delimit_pipolin(piPolB_coord, min(piPolB_coord)-10000, max(piPolB_coord)+10000, record, "No_delim")

        #Verify lactobacillaceae
        elif clu_3823_coords != "NA" and clu_2265_coords != "NA":
            n_alt_site += 1                       
            delimit_pipolin(piPolB_coord, clu_3823_coords, clu_2265_coords, record, "Alt")
            new_length_table += pipolin_id.replace(".gbk","")+"\t"+str(max(clu_3823_coords, clu_2265_coords)-min(clu_3823_coords, clu_2265_coords))+"\tAlt_site\n"

        #Rest
        else:
            n_nodelim += 1
            delimit_pipolin(piPolB_coord, min(piPolB_coord)-10000, max(piPolB_coord)+10000, record, "No_delim")

print(n_drs, n_alt_site, n_nodelim, n_no_full_pipolb)
print(n_drs + n_alt_site + n_nodelim + n_no_full_pipolb)

with open("Pipolin_delimited_new_int_site_length.tsv", "w") as f:
    f.write(new_length_table)