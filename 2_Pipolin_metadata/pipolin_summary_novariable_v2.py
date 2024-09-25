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
            piPolB_start = "NA"
            for record in SeqIO.parse("Filtered_pipolins/"+pipolin, "genbank"):
                pipolin_length += len(record.seq)

                for pipolin_cds in record.features:
                    if pipolin_cds.type == "CDS":
                        if pipolin_cds.qualifiers["product"][0] == "Primer-independent DNA polymerase PolB":
                            piPolB_start = pipolin_cds.location.start
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

            pipolin_contig_desc = [re.description for re in SeqIO.parse("Results_EP/"+genome_id+"/pipolins/"+pipolin_fa, "fasta")]
            pipolin_contig_ids = {}
            for desc in pipolin_contig_desc:
                g_key, coord_value = desc.split(' ',1) 
                pipolin_contig_ids[g_key] = coord_value

            #Record genome length [11], number of contigs [12], contig length list [13],  contig length sum [14], contig_name [15], contig_info [16], accesion_number [17], and pipolin_coordinates predicted by EP in contig AN [18] and contig custom code [19]
            #USES ORIGINAL PIPOLINS
            full_genome_length = 0
            pipolin_contig = []
            pipolin_contig_info = []
            pipolin_contig_accesion_number = []
            pipolin_contig_coordinates = []
            pipolin_contig_lenght = []
            pipolin_contig_sum = 0
            for fa in SeqIO.parse("Results_EP/"+genome_id+"/"+genome_id+".fa", "fasta"):
                full_genome_length += len(fa)
                if fa.id in list(pipolin_contig_ids.keys()):
                    pipolin_contig.append(str(fa.id))
                    pipolin_contig_info.append(" ".join(str(fa.description).split(" ")[2:]))
                    pipolin_contig_accesion_number.append(str(fa.description).split(" ")[1])
                    pipolin_contig_coordinates.append(str(fa.description).split(" ")[1]+" "+pipolin_contig_ids[fa.id])
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
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_coordinates))
            pipolin_dictionary[pipolin_id].append(str(pipolin_contig_desc))


            #Record repeats presence, number, coordinates, length, overlaps, 5' or 3'
            att_presence = False
            att_number = 0
            att_coordinates = []
            att_length = []
            att_mean_length = 0
            att_type = []
            att_overlaps = []
            att_overlaps_clean = []
            att_overlaps_feat = []
            att_target_end = "NA"


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
                                        att_overlaps_feat.append(feat)
                                        if str(feat.qualifiers["product"]).replace("'","").replace("[","").replace("]","") not in att_overlaps_clean:
                                            att_overlaps_clean.append(str(feat.qualifiers["product"]).replace("'","").replace("[","").replace("]",""))
                                    else: #In case it has no product
                                        att_overlaps.append(str(record.id)+"::"+str(feat.type))
                                        att_overlaps_feat.append(feat)
                                        if str(feat.type) not in att_overlaps_clean:
                                            print(str(feat.type), att_overlaps_clean)
                                            att_overlaps_clean.append(str(feat.type).replace("'","").replace("[","").replace("]",""))
            
            if len(att_overlaps_feat) == 1 and "product" in list(att_overlaps_feat[0].qualifiers.keys()):
                if 'tRNA-' in att_overlaps_feat[0].qualifiers["product"][0] or 'transfer-messenger RNA' in att_overlaps_feat[0].qualifiers["product"][0] :
                    if att_overlaps_feat[0].location.strand == 1:
                        if att_overlaps_feat[0].location.start < piPolB_start:
                            att_target_end = "3'"
                        elif piPolB_start < att_overlaps_feat[0].location.start:
                            att_target_end = "5'"
                    elif att_overlaps_feat[0].location.strand == -1:
                        if att_overlaps_feat[0].location.start < piPolB_start:
                            att_target_end = "5'"
                        elif piPolB_start < att_overlaps_feat[0].location.start:
                            att_target_end = "3'"
                    

            
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
            pipolin_dictionary[pipolin_id].append(att_target_end)


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




#print or write dict #add haeder 
table_output =  "Pipolin_ID\tGenome_ID\tGenone_num\tPipolin_file\tPipolin_scaffold\tPipolin_length\tpiPolB_CDS_num\tpiPolB_CDS_length\tpiPolB_CDS_length_sum\t"
table_output += "Genome_length\tContigs_length\tContigs_length_sum\tContigs_ids\tContig_record_type\tContig_description\tContig_AN\tContig_pipolin_coordinates_ExplorePipolin\tCode_pipolin_coordinates_ExplorePipolin\tAtts\tNumber_atts\t"
table_output += "Att_type\tCoordinates_atts\tLength_atts\tLength_atts_mean\tIntegration_site\tIntegration_site_clean\tTarget_end\t"
table_output += "Pipolin_gaps\tAssembly_gaps_paired_ends\tAssembly_gaps_pipolin_reconstruction\tPipolin_fragments\t"
table_output += "Assembly_AN\tAssembly_sequence_file\n"
for key in pipolin_dictionary:
    table_output += str(key)+"\t"+"\t".join(pipolin_dictionary[key])+"\n"


with open("pipolin_summary_postprocessed.txt", "w") as f:
    f.write(table_output)




