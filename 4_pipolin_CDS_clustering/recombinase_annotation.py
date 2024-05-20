import os
import subprocess

for hmm in os.listdir("hmms"):
    if hmm[-4:] == ".hmm":
        subprocess.run("hmmsearch -E 1E-5 --tblout results_"+hmm.replace(hmm[-4:],".tbl")+" hmms/"+hmm+" ../bacteria_dec2022_pipolin_cds_partial_filtered.faa",shell=True)

rec_list = {}
for file in os.listdir():
    if ".tbl" in file:
        rec = file.split(".")[0].replace("results_", "")
        with open(file, "r") as f:
            for line in f:
                if "G_" in line:
                    cds_id = line.split()[0]
                    rec_name = line.split()[2]
                    E_val = eval(line.split()[4])
                    if cds_id not in list(rec_list.keys()):
                        rec_list[cds_id] = {"Name":rec_name,"E-value":E_val}
                    else:
                        if E_val < rec_list[cds_id]["E-value"]:
                            rec_list[cds_id] = {"Name":rec_name,"E-value":E_val}

rec_list_output = ""
for cds_id in rec_list:
    rec_list_output += cds_id + "\t" + rec_list[cds_id]["Name"] + "\t" + str(rec_list[cds_id]["E-value"])  + "\n"

with open("hmmsearch_extra_rec_results.txt", "w") as f:
    f.write(rec_list_output)