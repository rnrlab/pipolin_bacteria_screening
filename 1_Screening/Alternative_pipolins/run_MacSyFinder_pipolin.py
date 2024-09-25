import os
from Bio import SeqIO
import subprocess
import itertools

'''
Part 1. Run MacSyFinder
'''

genomes_folder = "results_Enterobacterales_2024-07-02_16.24.48/91347"
output_table = "genome\treplicon	hit_id	gene_name	hit_pos	model_fqn	sys_id	sys_loci	locus_num	sys_wholeness	sys_score	sys_occ	hit_gene_ref	hit_status	hit_seq_len	hit_i_eval	hit_score	hit_profile_cov	hit_seq_cov	hit_begin_match	hit_end_match	counterpart	used_in\n"

for genome in os.listdir(genomes_folder+"_dataset/ncbi_dataset/data"):
    genome_path = genomes_folder+"_dataset/ncbi_dataset/data/"+genome
    if "GCA_" in genome and "protein.faa" in os.listdir(genome_path):
        
        #Order proteins according to gff order
        """
        n_prot = 1
        id_order = {}
        with open(genome_path+"/genomic.gff", "r") as fg:
            for line in fg:
                if "\tCDS\t" in line and "ID=cds-" in line and "pseudo=true" not in line:
                    cds_id = line.split("ID=cds-")[1].split(";")[0]
                    id_order[cds_id] = n_prot
                    n_prot += 1

        with open(genome_path+"/Order_list.tsv", "w") as fi:
            for id in id_order:
                fi.write(id+"\t"+str(id_order[id])+"\n")


        record_dict = {}
        for record in SeqIO.parse(genome_path+"/protein.faa", "fasta"):
            record_dict[str(record.id)] = record 

        multi_fasta_ordered = []
        for id in id_order:
            multi_fasta_ordered.append(record_dict[id])

        SeqIO.write(multi_fasta_ordered, genome_path+"/protein_ordered.faa", "fasta")
        """
        #Run MacSyFidner
        subprocess.run("macsyfinder --db-type ordered_replicon --sequence-db "+genome_path+"/protein_ordered.faa --models model_pipolin_ent_hhe -o "+genome_path+"/MSF_results_hhe -w 30 --force", shell=True)

        #Parse results
        with open(genome_path+"/MSF_results_hhe/best_solution.tsv", "r") as outMSF:
            for line in outMSF:
                if line[0] != "#" and "protein_ordered" in line:
                    output_table += genome+"\t"+line


with open(genomes_folder+"_pipolinMSF_hhe_search_results.tsv", "w") as fout:
    fout.write(output_table)