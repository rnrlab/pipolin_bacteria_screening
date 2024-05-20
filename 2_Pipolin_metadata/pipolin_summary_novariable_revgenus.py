import os
import glob
import re
from subprocess import check_call
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Entrez




### Read gbks from Results_EP/genome_id/pipolins/*single_record.gbk

pipolin_dictionary = {} #format: Genome_id : [] 


for pipolin in os.listdir("Filtered_pipolins/"):
        if ".gbk" in pipolin:

            genome_id = "G_"+pipolin.split("_")[1]
            pipolin_id = pipolin.split(".")[0].replace("_v0", "")

             #Record [pipolin_id](key-Field 0), [genome_id], [genome_num], [pipolin file], [filtering type],
            pipolin_dictionary[pipolin_id] = [genome_id]
            pipolin_dictionary[pipolin_id].append(genome_id.replace("G_",""))
            pipolin_dictionary[pipolin_id].append(pipolin)
            pipolin_dictionary[pipolin_id].append(pipolin.split(".")[1])
            
    
            #Record [pipolin length], pipolb_num, pipolb_fragments, pipolb_length_sum
            pipolin_length = 0
            num_pipolb = 0
            length_pipolb_list = []
            length_pipolb_sum = 0
            for record in SeqIO.parse("Filtered_pipolins/"+pipolin, "genbank"):
                pipolin_length += len(record.seq)

                for pipolin_cds in record.features:
                    if pipolin_cds.type == "CDS":
                        if pipolin_cds.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                            num_pipolb += 1
                            length_pipolb_list.append(len(pipolin_cds.qualifiers["translation"][0]))
                            length_pipolb_sum += len(pipolin_cds.qualifiers["translation"][0])
            
            pipolin_dictionary[pipolin_id].append(str(pipolin_length))
            pipolin_dictionary[pipolin_id].append(str(num_pipolb))
            pipolin_dictionary[pipolin_id].append(str(length_pipolb_list))
            pipolin_dictionary[pipolin_id].append(str(length_pipolb_sum))
            

            #get pipolin contig identifiers.
            for pipolin_old in os.listdir("Results_EP/"+genome_id+"/pipolins/"):
                if pipolin.split(".")[0] in pipolin_old and ".fa" in pipolin_old:
                    pipolin_fa = pipolin_old
            pipolin_contig_ids = [re.id for re in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin_fa, "fasta")]

            #Record genome length [11], number of contigs [12], contig length list [13],  contig length sum [14], contig_name [15], contig_info [16] and accesion_number [17]
            #USES ORIGINAL PIPOLINS
            full_genome_length = 0
            pipolin_contig = []
            pipolin_contig_info = []
            pipolin_contig_accesion_number = []
            pipolin_contig_lenght = []
            pipolin_contig_sum = 0
            for fa in SeqIO.parse("Results_EP/"+genome_id+"/"+genome_id+".fa", "fasta"):
                full_genome_length += len(fa)
                if fa.id in pipolin_contig_ids:
                    pipolin_contig.append(str(fa.id))
                    pipolin_contig_info.append(" ".join(str(fa.description).split(" ")[2:]))
                    pipolin_contig_accesion_number.append((fa.description).split(" ")[1])
                    pipolin_contig_lenght.append(str(len(fa)))
                    pipolin_contig_sum += len(fa)
            pipolin_dictionary[pipolin_id].append(str(full_genome_length))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_lenght))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_sum))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig))

            if "_p_" in str(pipolin_contig):
                pipolin_dictionary[pipolin_id].append("plasmid")
            else:
                pipolin_dictionary[pipolin_id].append("record")

            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_info))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_accesion_number))

            #Record repeats presence [18], number [19], coordinates [20], length [21] and overlaps [22]
            att_presence = False
            att_number = 0
            att_coordinates = []
            att_length = []
            att_mean_length = 0
            att_type = []
            att_overlaps = []
            att_overlaps_clean = []


            for record in SeqIO.parse("Filtered_pipolins/"+pipolin, "genbank"):
                if "trimmed" in pipolin: #we do not include de novo in this case
                    repeats = [feat for feat in record.features if (feat.type == "repeat_region" and feat.qualifiers["rpt_family"][0]=="Att" and feat.qualifiers["note"][0]=="pipolin conserved")]
                else:
                    repeats = [feat for feat in record.features if (feat.type == "repeat_region" and feat.qualifiers["rpt_family"][0]=="Att")]
                if repeats != []:
                    att_presence = True
                    att_number += len(repeats)

                    for repeat in repeats:
                        att_coordinates.append(str(record.id)+"::"+str(repeat.location))
                        att_length.append(int(abs(repeat.location.end-repeat.location.start)))

                        if repeat.qualifiers["note"][0] not in att_type:
                            att_type.append(repeat.qualifiers["note"][0])

                        att_range = range(int(repeat.location.start), int(repeat.location.end))
                        for feat in record.features:
                            if feat.type != "repeat_region" and feat.type != "source":
                                if feat.location.start in att_range or feat.location.end in att_range:
                                    if "product" in list(feat.qualifiers.keys()):
                                        att_overlaps.append(str(record.id)+"::"+str(feat.qualifiers["product"]).replace("'",""))
                                        if str(feat.qualifiers["product"]).replace("'","").replace("[","").replace("]","") not in att_overlaps_clean:
                                            att_overlaps_clean.append(str(feat.qualifiers["product"]).replace("'","").replace("[","").replace("]",""))
                                    else: #In case it has no product
                                        att_overlaps.append(str(record.id)+"::"+str(feat.type))
                                        if str(feat.type) not in att_overlaps_clean:
                                            print(str(feat.type), att_overlaps_clean)
                                            att_overlaps_clean.append(str(feat.type).replace("'","").replace("[","").replace("]",""))
            
            pipolin_dictionary[pipolin_id].append(str(att_presence))
            pipolin_dictionary[pipolin_id].append(str(att_number))
            pipolin_dictionary[pipolin_id].append(str(att_type))
            pipolin_dictionary[pipolin_id].append(str(att_coordinates))
            pipolin_dictionary[pipolin_id].append(str(att_length))
            if len(att_length) > 0:
                pipolin_dictionary[pipolin_id].append(str(sum(att_length)/len(att_length)))
            else:
                pipolin_dictionary[pipolin_id].append("0")
            pipolin_dictionary[pipolin_id].append(str(att_overlaps))
            if att_overlaps_clean == []:
                pipolin_dictionary[pipolin_id].append("NA")
            else:
                pipolin_dictionary[pipolin_id].append(", ".join(att_overlaps_clean))


            #Record [assembly gaps], [assembly gaps_paired-ends], [assembly gaps_pipolin_Structure], [Pipolin_fragments]
            reconstruction_gaps = 0
            gap_paired_ends = 0
            gap_pipolin_structure = 0

            for record in SeqIO.parse("Filtered_pipolins/"+pipolin, "genbank"):
                for pipolin_feat in record.features:
                    if "assembly_gap" in str(pipolin_feat.type):
                        reconstruction_gaps += 1
                        
                        if pipolin_feat.qualifiers["linkage_evidence"][0] == "paired-ends":
                            gap_paired_ends += 1

                        if pipolin_feat.qualifiers["linkage_evidence"][0] == "pipolin_structure":
                            gap_pipolin_structure += 1


            pipolin_dictionary[pipolin_id].append(str(reconstruction_gaps))
            pipolin_dictionary[pipolin_id].append(str(gap_paired_ends))
            pipolin_dictionary[pipolin_id].append(str(gap_pipolin_structure))
            pipolin_dictionary[pipolin_id].append(str(gap_pipolin_structure+1))



conversion_dict = {}
with open("conversion_table.txt", "r") as f:
    for line in f:
        conversion_dict[line.split("\t")[1].split(".")[0]] = line.strip("\n").split("\t")

#Record Genome AN and Genome File 
for key in pipolin_dictionary:
    pipolin_dictionary[key].append(conversion_dict[pipolin_dictionary[key][0]][2])
    pipolin_dictionary[key].append(conversion_dict[pipolin_dictionary[key][0]][3])


rev_genus = ["Aeromonas", "Citrobacter", "Corynebacterium", "Enterobacter", "Escherichia", "Limosilactobacillus", "Pseudosulfitobacter", "Salmonella", "Staphylococcus", "Vibrio"]

""" DEPRECARTED
### 1) Read genome information
genome_info = {}
with open("genome_metadata_pipolins_updatedCont.txt", "r") as f_genome:
    for line in f_genome:
        genome_fields = line.replace("\n","").split("\t")
        G = genome_fields[0]
        organism = genome_fields[2]
        species = genome_fields[16]
        genus = genome_fields[17]
        genome_info[G] = [organism, species, genus]
#Add organism, species name, assembly status from assembly stats file
for key in pipolin_dictionary:
    pipolin_dictionary[key].append(genome_info[pipolin_dictionary[key][0]][0])
    pipolin_dictionary[key].append(genome_info[pipolin_dictionary[key][0]][1])
    pipolin_dictionary[key].append(genome_info[pipolin_dictionary[key][0]][2])
    if genome_info[pipolin_dictionary[key][0]][2] in rev_genus:
        pipolin_dictionary[key].append(genome_info[pipolin_dictionary[key][0]][2])
    else:pipolin_dictionary[key].append("Other")
"""




#print or write dict #add haeder 
table_output =  "Pipolin_ID\tGenome_ID\tGenone_num\tPipolin_file\tFiltering_type\tpipolin_length\tpiPolB_num\tpiPolB_fragments_length\tpiPolB_fragments_length_sum\t"
table_output += "Genome_length\tContigs_length\tContigs_length_sum\tContigs_ids\tContig_record_type\tContig_description\tContig_AN\tAtts\tNumber_atts\t"
table_output += "Att_type\tCoordinates_atts\tLength_atts\tLength_atts_mean\tIntegration_site\tIntegration_site_clean\t"
table_output += "Asembly_gaps\tAssembly_gaps_paired_ends\tAssembly_gaps_pipolin_reconstruction\tPipolin_fragments\t"
table_output += "Genome_AN\tGenome_file\n"
for key in pipolin_dictionary:
    table_output += str(key)+"\t"+"\t".join(pipolin_dictionary[key])+"\n"



with open("pipolin_summary_postprocessed.txt", "w") as f:
    f.write(table_output)




"""
### OLD ACCESION TO BIOSAMPLE METADATA ETC
for key in pipolin_dictionary:
    n_count += 1
    Entrez.email = "vmc11298@gmail.com"
    AN = pipolin_dictionary[key][-2]
    
    check = False
    for i in range(5):
        search_id = Entrez.esearch(db="assembly", term=AN+"[Assembly Accession]")
        try:
            content = Entrez.read(search_id)
            id_num = content['IdList'][0]
            check = True
            break
        except:
            print("*** ERROR ***")
            print(AN, "There was an ERROR accesing Assembly databse. Retrying")
            print(pipolin_dictionary[key])
    
    if not check:
        print("Couldn't access, skipping.")
        continue

    
    fetch_id = Entrez.efetch(db="assembly", id=id_num, rettype="docsum", retmode="xml")
    content2 = Entrez.read(fetch_id)
    #Record species name [Field 20], organism [Fiel 21], genus [Field 22], AssemblyStatus [Field 23]
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["Organism"])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["Organism"].split()[0])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["SpeciesName"])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["AssemblyStatus"])
    print(key, AN)


    #BiosampleAN [Field 24], #BiosampleID [Field 25], #Biosource [Field 26], #Coverage [Field 27]
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["BioSampleAccn"])
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["BioSampleId"])
    pipolin_dictionary[key].append(str(content2['DocumentSummarySet']['DocumentSummary'][0]["Biosource"]))
    pipolin_dictionary[key].append(content2['DocumentSummarySet']['DocumentSummary'][0]["Coverage"])
    

    #Host and isolation Source, SRA and geo location
    biosample_id_num = content2['DocumentSummarySet']['DocumentSummary'][0]["BioSampleId"]
    fetch_id_biosample = Entrez.efetch(db="biosample", id=biosample_id_num, rettype="full", retmode="text")
    tmp_dict = {"host":"NA", "iso_source":"NA", "SRA":"NA", "geo_loc":"NA"}
    for line in fetch_id_biosample:
        if "/host=" in line:
            tmp_dict["host"] = line.split('"')[1]
        if "/isolation source=" in line:
            tmp_dict["iso_source"] = line.split('"')[1]
        if "SRA:" in line:
            tmp_dict["SRA"] = line.split("SRA:")[1].replace("\n", "")
        if "/geographic location=" in line:
            tmp_dict["geo_loc"] = line.split('"')[1]
        if "Error" in line:
            print(line)
            print("Error in", biosample_id_num)
    
    for tmp_key in tmp_dict:
        pipolin_dictionary[key].append(tmp_dict[tmp_key])
    
    print(n_count)
    print(pipolin_dictionary[key])

"""