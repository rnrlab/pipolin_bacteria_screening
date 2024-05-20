#!/bin/bash

#./mmseqs_clustering.sh bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70.faa bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db 30

INPUT="${1}" #bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70.faa
NAME=$(basename "${INPUT}" .faa)
RELEASE=$(date +%Y%m%d)
OUTPUT="${2}" #bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db/bacteria_dec2022_pipolin_cds_partial_filtered_mobileOG70_MMseqsClustered_db
NSLOTS="${3}" #30
MMSEQS=/usr/bin/mmseqs

rm -r  mmseqs_clustering
mkdir -p mmseqs_clustering


# create the seqDB for mmeseqs (for the cascade clustering)
"${MMSEQS}" createdb "${INPUT}" mmseqs_clustering/"${NAME}"_db

mkdir -p mmseqs_clustering/tmp1 
mkdir -p mmseqs_clustering/"${NAME}"_MMseqsClustered_db

# run the cascade clustering
"${MMSEQS}" cluster mmseqs_clustering/"${NAME}"_db mmseqs_clustering/"${OUTPUT}"_clu mmseqs_clustering/tmp1 --threads "${NSLOTS}" --min-seq-id 0.35 -c 0.6 --cluster-mode 0  -s 7.5  --cluster-steps 9  

# # compute output files
"${MMSEQS}" createtsv mmseqs_clustering/"${NAME}"_db mmseqs_clustering/"${NAME}"_db mmseqs_clustering/"${OUTPUT}"_clu mmseqs_clustering/"${OUTPUT}"_clu.tsv 
"${MMSEQS}" createseqfiledb mmseqs_clustering/"${NAME}"_db mmseqs_clustering/"${OUTPUT}"_clu mmseqs_clustering/"${OUTPUT}"_clu_seq
"${MMSEQS}" result2flat mmseqs_clustering/"${NAME}"_db mmseqs_clustering/"${NAME}"_db mmseqs_clustering/"${OUTPUT}"_clu_seq mmseqs_clustering/"${OUTPUT}"_clu.fasta


# #rearrange and compute some results
awk -f ${PWD}/convert_long_wide.awk mmseqs_clustering/"${OUTPUT}"_clu.tsv > mmseqs_clustering/"${OUTPUT}"_clu_wide.tsv
awk -F'\t' '{print NR"\t"$1"\t"NF-1}' mmseqs_clustering/"${OUTPUT}"_clu_wide.tsv > mmseqs_clustering/"${OUTPUT}"_clu_wide_cl_name_rep_size.tsv

# # extract cluster sequences
mkdir -p mmseqs_clustering/clusters
python3 DeepClustering_mmseqs2clu_fasta.py mmseqs_clustering/"${OUTPUT}"_clu.fasta

# compute stats and annotations
cd cluster_annotation_hhblits_pfam
python3 hhblits_pfam_annotation.py
cd ..
python3 Cluster_annotation.py
python3 Bipartite_network.py
