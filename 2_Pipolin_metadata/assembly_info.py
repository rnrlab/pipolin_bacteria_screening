import json
import os
from Bio import SeqIO

### Extra data
#a) G to assembly AN conversor
G_list = os.listdir("Results_EP")
A_list = []
G_to_assembly = {}

with open("conversion_table.txt", "r") as f:
    for line in f:
        if line.split()[1].replace(".fna", "") in G_list:
            G_to_assembly[line.split()[2]] = line.split()[1].replace(".fna", "")
            A_list.append(line.split()[2])


#b) Count number of pipolins per genome
G_to_pipolin_num = {}
for G in G_list:
    pipolin_num = 0
    for pipolin_file in os.listdir("Results_EP/"+G+"/pipolins"):
        if "single_record" in pipolin_file and "v0" in pipolin_file:
            pipolin_num += 1
    G_to_pipolin_num[G] = pipolin_num

#c) Count piPolB per genome
G_to_piPolB_num = {}
G_to_piPolB_sum = {}
for G in G_list:
    G_to_piPolB_num[G] = 0
    G_to_piPolB_sum[G] = 0
    for record in SeqIO.parse("Results_EP/"+G+"/"+G+"_piPolBs.faa", "fasta"):
        G_to_piPolB_num[G] += 1
        G_to_piPolB_sum[G] += len(record.seq)

"""
### PART A) FILTER METADATA


filtered_json = ""

with open('bacteria_dataset/ncbi_dataset/data/assembly_data_report.jsonl', 'r') as f:
    for line in f:
        line_dict = json.loads(line)
        if str(line_dict["accession"]) in A_list:
            filtered_json += line

with open("genome_metadata_pipolins.jsonl", "w") as fo2:
    fo2.write(filtered_json)


"""
### PART B) PROCESS METADATA

metadata_table = "G\tAN\torganism_full_name\torganism_short_name\tgenus\t"
metadata_table += "num_of_pipolins\t"
metadata_table += "isolation_source\thost\tgeo_loc_name\t"
metadata_table += "assembly_level\tassembly_name\tassembly_status\tassembly_type\t"
metadata_table += "biosample_accesion\t"
#metadata_table += "genome_length\t"
metadata_table += "strain\ttax_id\t"
metadata_table += "piPolB_num\tpiPolB_CDS_len_sum"
metadata_table += "\n"

n_fix = 0
n_fix_spec = 0
with open('genome_metadata_pipolins.jsonl', 'r') as f:
    for line in f:
        line_dict = json.loads(line)

        if str(line_dict["accession"]) in A_list:
            if len(line_dict["organism"]["organismName"].split(" ")) > 1:
                if "uncultured" in line_dict["organism"]["organismName"] or "Candidatus" in line_dict["organism"]["organismName"] or "bacterium" in line_dict["organism"]["organismName"].split(" ")[1]:
                    short_name = "NA"
                    genus = "NA"
                    n_fix += 1
                    #print(n_fix, line_dict["organism"]["organismName"])
                else:
                    short_name = " ".join(line_dict["organism"]["organismName"].split()[:2]).replace("[","").replace("]","")
                    genus = line_dict["organism"]["organismName"].split()[0].replace("[","").replace("]","")
            else:
                short_name = line_dict["organism"]["organismName"]
                genus = "NA"

            """
            if "bacterium" in line_dict["organism"]["organismName"].lower():
                n_fix_spec += 1
                print(n_fix_spec, line_dict["organism"]["organismName"])
            """
        
        isolation_source = "NA"
        host = "NA"
        geo_loc_name = "NA"

        for att_dict in line_dict["assemblyInfo"]["biosample"]["attributes"]:
            if "isolation_source" in list(att_dict.values()):
                isolation_source = att_dict["value"].lower()
            if "host" in list(att_dict.values()):
                host = att_dict["value"].lower()
            if "geo_loc_name" in list(att_dict.values()):
                geo_loc_name = att_dict["value"].split(":")[0].lower()
            
            
        strain = "NA"
        if "infraspecificNames" in list(line_dict["organism"].keys()):
            if "strin" in list(line_dict["organism"]["infraspecificNames"].keys()):
                strain = str(line_dict["organism"]["infraspecificNames"]["strain"])

        metadata_table += "\t".join([G_to_assembly[line_dict["accession"]], line_dict["accession"], line_dict["organism"]["organismName"], short_name, genus, str(G_to_pipolin_num[G_to_assembly[line_dict["accession"]]]),
                            isolation_source, host, geo_loc_name,
                            line_dict["assemblyInfo"]["assemblyLevel"], line_dict["assemblyInfo"]["assemblyName"], line_dict["assemblyInfo"]["assemblyStatus"], line_dict["assemblyInfo"]["assemblyType"],
                            #line_dict["assemblyInfo"]["bioprojectLineage"][0]["bioprojects"][0]["accession"], line_dict["assemblyInfo"]["bioprojectLineage"][0]["bioprojects"][0]["title"],
                            line_dict["assemblyInfo"]["biosample"]["accession"],
                            #line_dict["assemblyInfo"]["assemblyStats"]["totalSequenceLength"],
                            strain,str(line_dict["organism"]["taxId"]),
                            str(G_to_piPolB_num[G_to_assembly[line_dict["accession"]]]), str(G_to_piPolB_sum[G_to_assembly[line_dict["accession"]]])
                            ])
        metadata_table += "\n"

        

with open("genome_metadata_pipolins.txt", "w") as fo:
    fo.write(metadata_table)

