import sys
import pandas as pd
 
hit_check = 0
bidirectional_qn_sn = ""
tmp_id = "NA"
tmp_alis = {"lines":[],"qseqid_list":[],"sseqid_list":[]}
with open(sys.argv[1],"r") as f:
    for line in f: 
        if hit_check%100000 == 0 and hit_check!=0:
            print(hit_check)
        hit_id = line.split("\t")[-1].replace("\n","")
        if hit_id != tmp_id:
            #1) Close previous case
            hit_check += 1
            #1A: check if nq > 1 or ns > 1
            if len(set(tmp_alis["qseqid_list"]))!=len(tmp_alis["qseqid_list"]) or len(set(tmp_alis["sseqid_list"]))!=len(tmp_alis["sseqid_list"]):
                hit_check += 1
                #If True, compute n_q and n_s
                df = pd.DataFrame([ali.split("\t")[:-1] for ali in tmp_alis["lines"]], columns=["qseqid", "sseqid", "pident", "qcov", "scov", "bitscore"])
                q_n_df = df["qseqid"].value_counts()
                s_n_df = df["sseqid"].value_counts()
                #Now join dfs
                df_qn = pd.merge(df, q_n_df, on='qseqid', how='left')
                df_qn_sb = pd.merge(df_qn, s_n_df, on='sseqid', how='left')
                #Write results
                df_as_string = df_qn_sb.apply(lambda row: "\t".join(map(str, row)), axis=1)
                bidirectional_qn_sn += '\n'.join(df_as_string)+"\n"
            else:
                for line_nq1ns1 in tmp_alis["lines"]:
                    bidirectional_qn_sn += "\t".join(line_nq1ns1.split('\t')[:-1])+"\t1\t1\n"

            #2º) Open new case
            tmp_id = hit_id
            tmp_alis["lines"]=[line.replace("\n","")]
            tmp_alis["qseqid_list"]=[line.split("\t")[0]]
            tmp_alis["sseqid_list"]=[line.split("\t")[1]]

        else:
            tmp_alis["lines"].append(line.replace("\n",""))
            tmp_alis["qseqid_list"].append(line.split("\t")[0])
            tmp_alis["sseqid_list"].append(line.split("\t")[1])


#Add last case
if len(set(tmp_alis["qseqid_list"]))!=len(tmp_alis["qseqid_list"]) or len(set(tmp_alis["sseqid_list"]))!=len(tmp_alis["sseqid_list"]):
    hit_check += 1
    #If True, compute n_q and n_s
    df = pd.DataFrame([ali.split("\t")[:-1] for ali in tmp_alis["lines"]], columns=["qseqid", "sseqid", "pident", "qcov", "scov", "bitscore"])
    q_n_df = df["qseqid"].value_counts()
    s_n_df = df["sseqid"].value_counts()
    #Now join dfs
    df_qn = pd.merge(df, q_n_df, on='qseqid', how='left')
    df_qn_sb = pd.merge(df_qn, s_n_df, on='sseqid', how='left')
    #Write results
    df_as_string = df_qn_sb.apply(lambda row: "\t".join(map(str, row)), axis=1)
    bidirectional_qn_sn += '\n'.join(df_as_string)+"\n"
else:
    for line_nq1ns1 in tmp_alis["lines"]:
        bidirectional_qn_sn += "\t".join(line_nq1ns1.split('\t')[:-1])+"\t1\t1\n"
print(hit_check)
with open(sys.argv[1]+"_nq_ns", "w") as f:
    f.write(bidirectional_qn_sn)
