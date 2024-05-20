import os
import glob
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Entrez
#import warnings
#warnings.filterwarnings("error")

###Get pipolin information
pipolin_info = {}
with open("pipolin_summary.txt", "r") as f:
    for line in f:
        fields = line.replace("\n", "").split("\t")
        pipolin_info[fields[0].split("v")[0]] = fields

pass_pipolins = []
variable_pipolins = []


#MAX IS v_2 and v3 so we have 12 combinations (24 including single records)
table = "Pipolin_ID\t"
pipolin_length_dict = {}
pipolin_piPolB_dict = {}
pipolin_att_dict = {}

run_check = False #Set True for execution of the code below
### checking reconstruction lengths
if run_check:
    for genome_folder in glob.glob("Results_EP/*"):
        genome_id = genome_folder.replace("Results_EP/", "")
        for pipolin in os.listdir("Results_EP/"+genome_id+"/pipolins"):
            if ".gbk" in pipolin:
                pipolin_id = pipolin.split("v")[0]
                pipolin_len = 0
                piPolB_num = 0
                att_num = 0

                for record in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin, "genbank"):
                    pipolin_len += len(record.seq)
                    for pipolin_cds in record.features:
                        if pipolin_cds.type == "CDS":
                            if pipolin_cds.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                                piPolB_num += 1
                        if pipolin_cds.type == "repeat_region":
                            if pipolin_cds.qualifiers["rpt_family"][0] == "Att":
                                att_num += 1
                
                if pipolin_id not in list(pipolin_length_dict.keys()):
                    pipolin_length_dict[pipolin_id] = {}
                pipolin_length_dict[pipolin_id][pipolin] = pipolin_len

                if pipolin_id not in list(pipolin_piPolB_dict.keys()):
                    pipolin_piPolB_dict[pipolin_id] = {}
                pipolin_piPolB_dict[pipolin_id][pipolin] = piPolB_num

                if pipolin_id not in list(pipolin_att_dict.keys()):
                    pipolin_att_dict[pipolin_id] = {}
                pipolin_att_dict[pipolin_id][pipolin] = att_num

    print("count completed. Tabulating")

    pipolin_feature = ["length_", "att_num_", "pipolb_"]
    version_list = ["v0", "v1", "v2","v3"]
    record_type = [".multirecord", ".single_record"]

    for f in pipolin_feature:
        for v in version_list:
            for r in record_type:    
                table += f+"p"+v+r+"\t"

    table += "mean_len\tdiff_to_first\tmore_than_350\tmean_att\tdiff_att\tatt_cte\tmean_pipolb\tdiff_pipolb\tpipolb_cte\tfragments"
    table += "\n"

    for pipolin in pipolin_length_dict:
        table += pipolin+"\t"
        for f in pipolin_feature:
            for v in version_list:
                for r in record_type:

                    file_exists = False
                    if r == ".multirecord":
                        for file in pipolin_length_dict[pipolin]:
                            if v in file and ".single_record" not in file:
                                file_exists = True
                                if f == "length_":
                                    table += str(pipolin_length_dict[pipolin][file])+"\t"
                                elif f == "att_num_":
                                    table += str(pipolin_att_dict[pipolin][file])+"\t"
                                elif f == "pipolb_":
                                    table += str(pipolin_piPolB_dict[pipolin][file])+"\t"
                    else:
                        for file in pipolin_length_dict[pipolin]:
                            if  v in file and ".single_record" in file:
                                file_exists = True
                                if f == "length_":
                                    table += str(pipolin_length_dict[pipolin][file])+"\t"
                                elif f == "att_num_":
                                    table += str(pipolin_att_dict[pipolin][file])+"\t"
                                elif f == "pipolb_":
                                    table += str(pipolin_piPolB_dict[pipolin][file])+"\t"
                    
                    if not file_exists:
                        table += " \t"

        #add checks
        #1) check len
        mean_length = sum(pipolin_length_dict[pipolin].values())/len(pipolin_length_dict[pipolin])
        diff_to_first = abs(mean_length-list(pipolin_length_dict[pipolin].values())[0])
        greater_than_350 = diff_to_first < 350

        table += str(mean_length)+"\t"+str(diff_to_first)+"\t"+str(greater_than_350)+"\t"

        #2) check atts
        mean_atts = sum(pipolin_att_dict[pipolin].values())/len(pipolin_att_dict[pipolin])
        diff_to_first_atts = abs(mean_atts-list(pipolin_att_dict[pipolin].values())[0])
        greater_than_0_atts = diff_to_first_atts == 0

        table += str(mean_atts)+"\t"+str(diff_to_first_atts)+"\t"+str(greater_than_0_atts)+"\t"

        #3) check pipolb
        mean_pipolb = sum(pipolin_piPolB_dict[pipolin].values())/len(pipolin_piPolB_dict[pipolin])
        diff_to_first_pipolb = abs(mean_pipolb-list(pipolin_piPolB_dict[pipolin].values())[0])
        greater_than_0_pipolb = diff_to_first_pipolb == 0

        table += str(mean_pipolb)+"\t"+str(diff_to_first_pipolb)+"\t"+str(greater_than_0_pipolb)

        ### store pipolb id in pass or fail
        if greater_than_350 and greater_than_0_atts and greater_than_0_pipolb:
            pass_pipolins.append(pipolin)
        else:
            variable_pipolins.append(pipolin)

        table += "\n"

    with open("pipolin_versions_check.txt", "w") as f:
        f.write(table)

    with open("variable_pipolins.txt", "w") as fo:
        output_list = ""
        for pipolin_id in variable_pipolins:
            output_list += pipolin_id + "\n"
        fo.write(output_list)



"""
#This is for adding new info to the table
new_table_len = ""
with open("pipolin_length_versions.txt", "r") as f:
    for line in f:
        fields = line.replace("\n", "").split("\t")
        if "G_" in line:
            new_table_len += line.replace("\n",pipolin_info[fields[0]][6]+"\t"+pipolin_info[fields[0]][5]+"\n")
        else:
            new_table_len += line.replace("\n","fragments\tpipolin_length\n")

with open("pipolin_length_versions_calc_updated.txt", "w") as f:
    f.write(new_table_len)
"""

### checking reconstruction lengths and number contigs with pipolb. Filter according to that info 

#nums for stats
n_check = 0
n_invariable = 0
n_invariable_unchanged = 0
n_invariable_trimmed = 0
n_variable = 0
n_variable_oneContig = 0
n_variable_oneContig_trimmed = 0
n_variable_discarded = 0

run_filter = True
if run_filter:

    if "Filtered_pipolins" not in os.listdir():
        subprocess.run("mkdir Filtered_pipolins", shell=True)
    else:    
        subprocess.run("rm -r Filtered_pipolins", shell=True)
        subprocess.run("mkdir Filtered_pipolins", shell=True)

    variable_pipolins = []
    with open("variable_pipolins.txt", "r") as f:
        for line in f:
            variable_pipolins.append(line.replace("\n",""))

    for genome_folder in glob.glob("Results_EP/*"):
        genome_id = genome_folder.replace("Results_EP/", "")
        for pipolin in os.listdir("Results_EP/"+genome_id+"/pipolins"):
            pipolin_id = pipolin.split("v")[0]
            
            n_record_pipolb = 0
            pipolbs = {}
            pipolbs_id_to_coord = {}
            
            #Case A) Pipolin variable -> Requires selecting contig that has the piPolB and cut 30 kbp at each side
            if pipolin_id in variable_pipolins and ".gbk" in pipolin and "v0" in pipolin and "single_record" not in pipolin:
                n_check += 1
                n_variable += 1
                #First we count the number of piPolBs and its length
                for record in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin, "genbank"):
                    #if len(record.seq) > 5000: #We only consider contigs bigger than 5 kbp
                        pipolbs[record.name] = 0
                        record_with_pipolb = False
                        for pipolin_cds in record.features:
                            if pipolin_cds.type == "CDS":
                                if pipolin_cds.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                                    if not record_with_pipolb:
                                        n_record_pipolb += 1
                                        min_coord_pipolb = int(pipolin_cds.location.start)
                                        record_with_pipolb = True
                                    pipolbs[record.name] += len(pipolin_cds)
                                    max_coord_pipolb = int(pipolin_cds.location.end)
                        pipolbs_id_to_coord[record.name] = [min_coord_pipolb, max_coord_pipolb]

                #Case A.1) 1 or more piPolbs 
                if n_record_pipolb > 0:
                    n_variable_oneContig += 1
                    max_pipolb_record = max(pipolbs, key=pipolbs.get) #pick contig with biggest piPolB
                    for record in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin, "genbank"):
                        if max_pipolb_record == record.name:
                            
                            #left cut adjust
                            if pipolbs_id_to_coord[record.name][0] - 30000 < 0:
                                left_cut = 0
                            else: 
                                left_cut = pipolbs_id_to_coord[record.name][0] - 30000

                            #right_cut adjust
                            if pipolbs_id_to_coord[record.name][1] + 30000 > len(record.seq):
                                right_cut = len(record.seq)+1
                            else:
                                right_cut = pipolbs_id_to_coord[record.name][1] + 30000

                            if left_cut != 0 or right_cut != len(record.seq)+1:
                                n_variable_oneContig_trimmed += 1
                            
                            with open("Filtered_pipolins/"+pipolin.split(".")[0]+".trimmed.fa", "w") as fasta_out:
                                fasta_out.write(">"+pipolin.split(".")[0]+"\t["+str(left_cut)+":"+str(right_cut)+"]\n"+str(record.seq)[left_cut:right_cut]+"\n")
                            
                            SeqIO.write(record[left_cut:right_cut], "Filtered_pipolins/"+pipolin.split(".")[0]+".trimmed.gbk", "genbank")
                #Case A.2) 0 pipolbs in v0 -> discarded (2 pipolins)
                else:
                    n_variable_discarded += 1
                    print("No piPolB in ", pipolin)
            
            #Case B) Pipolin length invariable -> Extract sequece from single_record
            if pipolin_id not in variable_pipolins and "v0" in pipolin and "single_record.gbk" in pipolin:
                n_check += 1
                n_invariable += 1
                for record in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin, "genbank"):
                    
                    #obtain pipolb coordinates
                    pipolbs[record.name] = 0
                    record_with_pipolb = False
                    for pipolin_cds in record.features:
                        if pipolin_cds.type == "CDS":
                            if pipolin_cds.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                                if not record_with_pipolb:
                                    n_record_pipolb += 1
                                    min_coord_pipolb = int(pipolin_cds.location.start)
                                    record_with_pipolb = True
                                pipolbs[record.name] += len(pipolin_cds)
                                max_coord_pipolb = int(pipolin_cds.location.end)
                    pipolbs_id_to_coord[record.name] = [min_coord_pipolb, max_coord_pipolb]

                    if len(record.seq) > 100000:
                        n_invariable_trimmed += 1
                        #left cut adjust
                        if pipolbs_id_to_coord[record.name][0] - 30000 < 0:
                            left_cut = 0
                        else: 
                            left_cut = pipolbs_id_to_coord[record.name][0] - 30000

                        #right_cut adjust
                        if pipolbs_id_to_coord[record.name][1] + 30000 > len(record.seq):
                            right_cut = len(record.seq)+1
                        else:
                            right_cut = pipolbs_id_to_coord[record.name][1] + 30000
                        

                        with open("Filtered_pipolins/"+pipolin.split(".")[0]+".trimmed.fa", "w") as fasta_out:
                            fasta_out.write(">"+pipolin.split(".")[0]+"\t["+str(left_cut)+":"+str(right_cut)+"]\n"+str(record.seq)[left_cut:right_cut]+"\n")
                            
                        SeqIO.write(record[left_cut:right_cut], "Filtered_pipolins/"+pipolin.split(".")[0]+".trimmed.gbk", "genbank")

                    else:
                        n_invariable_unchanged += 1 
                        with open("Filtered_pipolins/"+pipolin.split(".")[0]+".original.fa", "w") as fasta_out:
                            fasta_out.write(">"+pipolin.split(".")[0]+"\n"+str(record.seq)+"\n")
                        SeqIO.write(record, "Filtered_pipolins/"+pipolin.split(".")[0]+".original.gbk", "genbank")               

print("processed:", n_check)
print("N_invariable: ", n_invariable)
print("n_invariable_unchanged: ", n_invariable_unchanged)
print("n_invariable_trimmed: ", n_invariable_trimmed)
print("n_variable: ", n_variable)
print("n_variable_oneContig: ", n_variable_oneContig)
print("n_variable_oneContig_trimmed: ", n_variable_oneContig_trimmed)
print("n_variable_discarded: ", n_variable_discarded)


