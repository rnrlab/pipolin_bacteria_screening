#!/bin/bash

#0) Set variables and create folders
INPUT_DB="${1}" #pipolin_test
OUTPUT=$(basename "${INPUT_DB}")


rm -r  linclust_all_mge
mkdir -p linclust_all_mge

#1) Filter dataset
seqkit seq -m 50 "${INPUT_DB}".faa > "${INPUT_DB}"_50.faa 

#1.1 Create database
mmseqs createdb "${INPUT_DB}"_50.faa linclust_all_mge/"${OUTPUT}"_db 

#1.2) Cluster database
mmseqs linclust linclust_all_mge/"${OUTPUT}"_db  linclust_all_mge/"${OUTPUT}"_db_clu linclust_all_mge/tmp -c 0.75 --min-seq-id 0.75 -e 0.001 --cluster-mode 1  
mmseqs createtsv linclust_all_mge/"${OUTPUT}"_db linclust_all_mge/"${OUTPUT}"_db linclust_all_mge/"${OUTPUT}"_db_clu linclust_all_mge/"${OUTPUT}"_db_clu.tsv 


#2) Get interacting MGEs
python3 filter_clusters.py linclust_all_mge/"${OUTPUT}"_db_clu.tsv "${OUTPUT}".faa

#3) MMSeqs all vs all 
rm -r  mmseqs_all_vs_all
mkdir -p mmseqs_all_vs_all

#3.1 Encode proteins
python3 CDS_id_encoding.py "${INPUT_DB}"_filtered.faa

#3.2 Create protein database
mmseqs createdb "${INPUT_DB}"_filtered_encoded.faa mmseqs_all_vs_all/"${OUTPUT}"_db #--dbtype 1 #_encoded

rm "${INPUT_DB}"_50.faa 

#3.3 Cluster full database
mmseqs cluster mmseqs_all_vs_all/"${OUTPUT}"_db mmseqs_all_vs_all/"${OUTPUT}"_db_clu mmseqs_all_vs_all/tmp --min-seq-id 0.35 -c 0.5  -s 7.5   #alternative: --cluster-mode 1 --single-step-clustering
mmseqs createtsv mmseqs_all_vs_all/"${OUTPUT}"_db mmseqs_all_vs_all/"${OUTPUT}"_db mmseqs_all_vs_all/"${OUTPUT}"_db_clu mmseqs_all_vs_all/"${OUTPUT}"_clu.tsv
mmseqs createseqfiledb mmseqs_all_vs_all/"${OUTPUT}"_db mmseqs_all_vs_all/"${OUTPUT}"_db_clu  mmseqs_all_vs_all/"${OUTPUT}"_db_clu_seq
mmseqs result2flat mmseqs_all_vs_all/"${OUTPUT}"_db mmseqs_all_vs_all/"${OUTPUT}"_db mmseqs_all_vs_all/"${OUTPUT}"_db_clu_seq mmseqs_all_vs_all/"${OUTPUT}"_clu.fasta

awk -f convert_long_wide.awk mmseqs_all_vs_all/"${OUTPUT}"_clu.tsv > mmseqs_all_vs_all/"${OUTPUT}"_clu_wide.tsv
awk -F'\t' '{print NR"\t"$1"\t"NF-1}' mmseqs_all_vs_all/"${OUTPUT}"_clu_wide.tsv > mmseqs_all_vs_all/"${OUTPUT}"_clu_wide_cl_name_rep_size.tsv

mkdir -p mmseqs_all_vs_all/clusters
python3 mmseqs2clu_fasta.py mmseqs_all_vs_all/"${OUTPUT}"_clu.fasta

#3.4 carry out all vs all comparisons
python3 all_vs_all_cluster.py

echo "mmseqs finished"

#4) Process output
echo "Adding hit id (alpha ord) and sort by it"
for file in mmseqs_all_vs_all/clusters_all_vs_all/Clu_*/*_ali.m8; do
    cat "$file" | awk -F'\t' '$1!=$2' | awk 'BEGIN{FS=OFS="\t"} {print $0, ($1 <= $2 ? $1 "_" $2 : $2 "_" $1)}' > "$file"_hitID
    sort -V -t $'\t' -k7 --parallel=30 "$file"_hitID > "$file"_hitID_sorted
done

echo "Computing best hits"
for ((i=1; i<=9; i++)); do
    ls mmseqs_all_vs_all/clusters_all_vs_all/Clu_${i}*/*_ali.m8_hitID_sorted | parallel --verbose 'python3 Bidirectional_hits.py {}'
done

echo "Sorting by pair id (MGE-MGE) "
for file in mmseqs_all_vs_all/clusters_all_vs_all/Clu_*/*_ali.m8; do
    sort -V -t $'\t' -k7 --parallel=30 "$file"_hitID_sorted_BH > "$file"_hitID_sorted_BH_sorted
done

echo "Counting multihits"
for ((i=1; i<=9; i++)); do
    ls mmseqs_all_vs_all/clusters_all_vs_all/Clu_${i}*/*_ali.m8_hitID_sorted_BH_sorted | parallel --verbose 'python3 Bidirectional_hits_count_nq_ns.py {}'
done

echo "Splitting BH in single and multi"
for file in mmseqs_all_vs_all/clusters_all_vs_all/Clu_*/*_ali.m8; do
    echo "$file"
    cat "$file"_hitID_sorted_BH_sorted_nq_ns | awk -F'\t' '$7==1 && $8==1' | cut -f 1,2,3,4,5 > "$file"_BBH_single_hit
    cat "$file"_hitID_sorted_BH_sorted_nq_ns | awk -F'\t' '$7+$8>2' > "$file"_BBH_multi_hit
    rm "$file"_hitID "$file"_hitID_sorted "$file"_hitID_sorted_BH "$file"_hitID_sorted_BH_sorted
done

#5) Solve multihits 
echo "Solving multihits"
rm -f BBH_multi_hit.tsv
for mh_file in mmseqs_all_vs_all/clusters_all_vs_all/Clu_*/*multi_hit; do
    cat "$mh_file" >> BBH_multi_hit.tsv
done

#R --vanilla < Bidirectional_multihits_processing.R 
cat BBH_multi_hit.tsv | awk '{print $0"\t"$1"_"$2"\t"$2"_"$1}' > BBH_multi_hit_ids.tsv
sort -V -t $'\t' -k9 --parallel=30 BBH_multi_hit_ids.tsv > BBH_multi_hit_ids_sorted_q.tsv
sort -V -t $'\t' -k10 --parallel=30 BBH_multi_hit_ids.tsv > BBH_multi_hit_ids_sorted_s.tsv
python3 Bidirectional_multihits_processing.py BBH_multi_hit_ids_sorted_q.tsv query
python3 Bidirectional_multihits_processing.py BBH_multi_hit_ids_sorted_s.tsv subject 
cat BBH_multi_hit_ids_sorted_qquery_solved.tsv BBH_multi_hit_ids_sorted_ssubject_solved.tsv > BBH_multi_hit_ids_solved_all.tsv
R --vanilla < Bidirectional_multihits_processing_part2.R
rm BBH_multi_hit_ids.tsv BBH_multi_hit_ids_sorted_qquery_solved.tsv BBH_multi_hit_ids_sorted_ssubject_solved.tsv



#6) Join with single hits
rm -f BBHsg.tsv
for gh_file in mmseqs_all_vs_all/clusters_all_vs_all/Clu_*/*_single_hit; do
    cat "$gh_file" >> BBHsg.tsv
done

echo "Merging with single-hits"
cat BBH_multi_hit_solved.tsv BBHsg.tsv > BBH.tsv
rm BBHsg.tsv

#6) Compute wGRR
echo "Computing wGRR"
R --vanilla < wGRR_short.R

