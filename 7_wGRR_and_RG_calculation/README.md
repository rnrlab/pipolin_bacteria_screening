## 7 - wGRR computation and RG/NRG classification

This folder includes all scripts necessary for wGRR calculation between MGEs and classification of their CDSs into recombining (RGs) and non-recombining genes (NRGs).

Computation of wGRR values is carried out by "Full_wGRR_script.sh", which will call the remaining scripts in the folder (which in turn do specific tasks). To name a few (the most important):
- CDS_id_encoding.py: encodes MGE ids in a three letter code to reduce memory usage.
- filter_clusters.py: keeps clusters of MGE genes where at least 3 pipolins (and 3% of cluster) are involved. 
- all_vs_all_cluster.py: performs the all-vs-all aligments within clusters necessary to calculate wGRR. 

RG/NRG classification and plotting scripts are in a separate folder (7b_RG_assignation)
