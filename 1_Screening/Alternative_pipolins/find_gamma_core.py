import os
import subprocess
import itertools

'''
Part 1. Find gammaproteobacteria pipolin core genes
'''

output_table = "Genome\thmmsearch_hits_num\thmmsearch_hits_info\tgene_num\tgene_coords\tmin_coord\tmax_coord\tdistance\tpiPolB_present\t"+"\t".join(["piPolB_hmm_12-02-2021","Clu2-IntSXT","Clu3-DUF2787","Clu4-Hyp1","Clu11-IntP2","Clu14-Hyp2","Clu29-Hyp2","Clu32-IntSXT","Clu53-IntP2"])+"\n"

genomes_folder = "results_Enterobacterales_2024-07-02_16.24.48/91347"
for genome in os.listdir(genomes_folder+"_dataset/ncbi_dataset/data"):
    genome_path = genomes_folder+"_dataset/ncbi_dataset/data/"+genome
    if "GCA_" in genome and "protein.faa" in os.listdir(genome_path):
        proteins_path = genome_path+"/protein.faa"
        
        #hmmsearch
        subprocess.run("hmmsearch --tblout "+genome_path+"/"+genome+"_hmmsearch_gamma_core.tbl --domtblout "+genome_path+"/"+genome+"_hmmsearch_gamma_core.domtbl -E 1e-30 --cpu 30  hmms_cluster/gamma_core.hmm "+proteins_path, shell=True, stdout = subprocess.DEVNULL)

        #process output, i.e. select best match for each profile
        match_list = {}
        piPolB = "No"
        with open(genome_path+"/"+genome+"_hmmsearch_gamma_core.tbl", "r") as f1:
            for line in f1:
                if line[0] != "#":
                    hmmsearch_fields = line.split()
                    line_hmm = hmmsearch_fields[2]
                    line_prot_id = hmmsearch_fields[0]
                    line_eval = hmmsearch_fields[4]
                    if line_hmm not in list(match_list.keys()): #add match if profile not present
                        match_list[line_hmm] = [line_prot_id, line_eval] #hmm_profile: [protein_id, e-value]
                    else: #replace by better match if lower evalue
                        if eval(line_eval) < eval(match_list[line_hmm][1]):
                            match_list[line_hmm] = [line_prot_id, line_eval]
                    
                    if line_hmm == "piPolB_hmm_12-02-2021":
                        piPolB = "Yes"

        profile_best_evalue = {"piPolB_hmm_12-02-2021":"NA", "Clu2-IntSXT":"NA", "Clu3-DUF2787":"NA", "Clu4-Hyp1":"NA", "Clu11-IntP2":"NA", "Clu14-Hyp2":"NA", "Clu29-Hyp2":"NA", "Clu32-IntSXT":"NA", "Clu53-IntP2":"NA"}
        for hmm_profile_name in match_list:
            print(match_list[hmm_profile_name][1])
            profile_best_evalue[hmm_profile_name] = match_list[hmm_profile_name][1]
        
        
        
        if match_list != {}:
            #Add coordinates
            prot_id_list = [value[0] for value in match_list.values()]
            prot_id_coords = {}
            with open(genome_path+"/genomic.gff", "r") as f2:
                for line in f2:
                    if "\tCDS\t" in line and line[0]!="#":
                        prot_id_gff = line.split("ID=cds-")[1].split(";")[0]
                        if prot_id_gff in prot_id_list:
                            prot_id_coords[prot_id_gff] = [line.split("\t")[3],line.split("\t")[4]]
            
            #Calculate stats
            min_coord = min([int(value) for sublist in prot_id_coords.values() for value in sublist])
            max_coord = max([int(value) for sublist in prot_id_coords.values() for value in sublist])
            print(min_coord)
            distance = max_coord - min_coord
        
            output_table += "\t".join([genome,str(len(match_list)),str(match_list),str(len(prot_id_coords)), str(prot_id_coords), str(min_coord), str(max_coord), str(distance), piPolB,
                                       "\t".join(list(profile_best_evalue.values()))])+"\n"
            print(genome, len(match_list),match_list)
            print(len(prot_id_coords), prot_id_coords)
            print(min_coord, max_coord, distance)
            print("\n")

with open(genomes_folder+"_Gamma_core_search_results.tsv", "w") as fout:
    fout.write(output_table)


