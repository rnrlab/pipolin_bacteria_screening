import os
import subprocess
"""
subprocess.run("mkdir MSA_Phrogs_M50_FASTA_func", shell=True)
subprocess.run("mkdir PHROGS_hmm_func", shell=True)



n = 0
with open("phrog_annot_v4.tsv", "r") as f:
    for line in f:
        if "category" not in line and "unknown function" not in line and line.split("\t")[0]:
            n += 1
            subprocess.run("cp MSA_Phrogs_M50_FASTA/phrog_"+line.split("\t")[0]+".fma MSA_Phrogs_M50_FASTA_func/", shell=True)


subprocess.run('ls MSA_Phrogs_M50_FASTA_func/* | parallel --verbose "hmmbuild PHROGS_hmm_func/{/.}.hmm {}"',  shell=True)
"""
subprocess.run("cat PHROGS_hmm_func/* > PHROGS_func.hmm", shell=True)
subprocess.run("hmmpress PHROGS_func.hmm", shell=True)
subprocess.run("hmmsearch -E 1E-3 --cpu 30 --domtblout All_MGE_5-500kbp_RefSeq_PHROGS.domtblout PHROGS_func.hmm ../../MGE_pipolins.faa", shell=True)
