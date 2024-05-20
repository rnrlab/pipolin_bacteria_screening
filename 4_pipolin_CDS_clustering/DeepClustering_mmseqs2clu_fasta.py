from collections import OrderedDict
import sys

#python3 DeepClustering_mmseqs2clu_fasta.py ${DBNAME}.fasta
#e.g. python3 DeepClustering_mmseqs2clu_fasta.py ${DBNAME}_deep_clustered_db.fasta

a=open(sys.argv[1],"r")
dicc={}
prev=""
for n in a:
	if n.startswith(">"):
		if prev.startswith(">"):
			prot = n
			dicc[n]=[]
	else:
		dicc[prot].append(prev)
		dicc[prot].append(n)
	prev=n

clu=0
for k in sorted(dicc, key=lambda k: len(dicc[k]), reverse=True):
	clu=clu+1 #old name
	clu_name = k.split()[0].split("|")[0].replace(">","")
	if len([n for n in dicc[k] if n.startswith(">")])>= 5:
		w = open("mmseqs_clustering/clusters/Cluster_"+str(clu)+"_"+clu_name+".fa","w+")
		for n in dicc[k]:
			w.write(n)
		w.close()
