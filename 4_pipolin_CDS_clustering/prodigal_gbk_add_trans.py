import os	
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

subprocess.run("mkdir -p pipolin_prodigal_cds_gbk", shell=True)
subprocess.run("rm -r pipolin_prodigal_cds_gbk/*", shell=True)


for faa_file in os.listdir("pipolin_prodigal_cds_faa"):
    
    #A) get dna seq record
    dna_seq_file = "../pipolin_mges_graph/Filtered_pipolins_5k/"+faa_file.replace(".prodigal.faa",".fa")
    dna_record = [record for record in SeqIO.parse(dna_seq_file, "fasta")][0]#This is the dna seq record
    dna_record.annotations["molecule_type"] = "DNA"

    #B) add prodigal cds record as seqfeature to dna record
    for protein_record in SeqIO.parse("pipolin_prodigal_cds_faa/"+faa_file, "fasta"):
        protein_record_description = str(protein_record.description).split(" # ")
        start_n = int(protein_record_description[1])
        end_n = int(protein_record_description[2])
        protein_feature = SeqFeature(FeatureLocation(start=start_n, end=end_n), type='CDS',qualifiers={"translation":str(protein_record.seq).replace("*","")})
        dna_record.features.append(protein_feature)

    #C) Write gbks
    SeqIO.write(dna_record, "pipolin_prodigal_cds_gbk/"+faa_file.replace(".faa",".gbk"), 'genbank')


       