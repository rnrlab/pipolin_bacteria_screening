import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

filtered_fa_total = ""
fitlered_fa_subset = ""
n_subset = 0 #Subset for EggNog mapper
n_CDS_subset = 0
len_list = []

for record in SeqIO.parse("bacteria_dec2022_pipolin_cds_with_partial.faa", "fasta"):
    len_list.append(len(record.seq))
    if len(record.seq) >= 30 and len(record.seq) <= 2000:
        n_CDS_subset += 1
        filtered_fa_total += ">"+record.id+"\t"+record.description+"\t"+str(len(record.seq)-1)+"\n"+str(record.seq).replace("*","")+"\n"
        fitlered_fa_subset += ">"+record.id+"\t"+record.description+"\t"+str(len(record.seq)-1)+"\n"+str(record.seq).replace("*","")+"\n"

        if n_CDS_subset == 99000:
            n_subset += 1
            n_CDS_subset = 0
            with open("bacteria_dec2022_pipolin_cds_partial_filtered_subset_"+str(n_subset)+".faa", "w") as f:
                f.write(fitlered_fa_subset)
            fitlered_fa_subset = ""


#last subset
n_subset += 1
n_CDS_subset = 0
with open("bacteria_dec2022_pipolin_cds_partial_filtered_subset_"+str(n_subset)+".faa", "w") as f:
    f.write(fitlered_fa_subset)
fitlered_fa_subset = ""


with open("bacteria_dec2022_pipolin_cds_partial_filtered.faa", "w") as f:
    f.write(filtered_fa_total)


len_list_arr = np.array(len_list)
print(len(len_list_arr[len_list_arr>1500]))
"""
plt.hist(len_list_arr, bins=5000)
plt.show()
"""