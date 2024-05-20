import os
import subprocess

"""
for hmm in os.listdir("hmms"):
    if hmm[-4:] == ".hmm":
        subprocess.run("hmmsearch -E 1E-3 --cpu 30 --domtblout results_"+hmm.replace(hmm[-4:],".domtbl")+" hmms/"+hmm+" ../../MGE_pipolins.faa",shell=True)
"""

subprocess.run("cat *.domtbl > MGE_pipolin_YR-SR.domtbl",shell=True)


domtbl_out = "Gene_ID\tRecombinase\tScore\n"

with open("MGE_pipolin_YR-SR.domtbl", "r") as f:
    for line in f:
        if "#" != line[0]:
            cds_id = line.split()[0]
            tlen = line.split()[2]
            rec_name = line.split()[3]
            E_val = eval(line.split()[6])
            score = line.split()[7]
            qcov = (eval(line.split()[18])-eval(line.split()[17]))/eval(tlen)
            if qcov >= 0.5 and E_val<=0.001:
                domtbl_out += cds_id+"\t"+rec_name+"\t"+score+"\n"



with open("MGE_pipolin_YR-SR_results.txt", "w") as f:
    f.write(domtbl_out)
