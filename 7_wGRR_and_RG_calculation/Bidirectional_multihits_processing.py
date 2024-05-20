import sys
import pandas as pd

direction = sys.argv[2]

hit_check = 0
BH_solved = ""
tmp_id = "NA"
tmp_alis = {"lines":[],"qseqid_list":[],"sseqid_list":[]}
with open(sys.argv[1],"r") as f:
    for line in f: 
        if hit_check%100000 == 0 and hit_check!=0:
            print(hit_check)

        if direction == "query":
            hit_id = line.split("\t")[-2].replace("\n","")
        elif direction == "subject":
            hit_id = line.split("\t")[-1].replace("\n","")
        hit_to_mge = "_".join(hit_id.split("_")[:3])


        if hit_to_mge != tmp_id:
            #1) Close previous case
            hit_check += 1
            if tmp_alis["lines"] != []:
                #1A: create dataframe for the multihit case
                df = pd.DataFrame([ali.split("\t")[:-4] for ali in tmp_alis["lines"]], columns=["qseqid", "sseqid", "pident", "qcov", "scov", "bitscore"])
                max_row = df.loc[df['bitscore'].idxmax()]
                #Write results
                max_row_as_string = "\t".join(max_row.astype(str))
                BH_solved += max_row_as_string+"\n"

            #1B) Open new case
            tmp_id = hit_to_mge
            tmp_alis["lines"]=[line.replace("\n","")]

        else:
            tmp_alis["lines"].append(line.replace("\n",""))


#Last case
df = pd.DataFrame([ali.split("\t")[:-4] for ali in tmp_alis["lines"]], columns=["qseqid", "sseqid", "pident", "qcov", "scov", "bitscore"])
max_row = df.loc[df['bitscore'].idxmax()]
max_row_as_string = "\t".join(max_row.astype(str))
BH_solved += max_row_as_string+"\n"

with open(sys.argv[1].replace(".tsv",direction+"_solved.tsv"),"w") as f:
    f.write(BH_solved)
