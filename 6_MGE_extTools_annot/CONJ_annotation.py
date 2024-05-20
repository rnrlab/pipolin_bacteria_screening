import os
import subprocess


domtbl_out = "Gene_ID\tCONJ_hmm\tScore\n"

with open("MGE_pipolin_CONJScan.domtblout", "r") as f:
    for line in f:
        if "#" != line[0]:
            cds_id = line.split()[0]
            tlen = line.split()[2]
            q_name = line.split()[3]
            E_val = eval(line.split()[6])
            score = line.split()[7]
            qcov = (eval(line.split()[18])-eval(line.split()[17]))/eval(tlen)
            if qcov >= 0.5 and E_val<=0.001:
                domtbl_out += cds_id+"\t"+q_name+"\t"+score+"\n"



with open("MGE_pipolin_CONJ_results.txt", "w") as f:
    f.write(domtbl_out)
