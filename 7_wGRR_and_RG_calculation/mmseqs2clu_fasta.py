from collections import OrderedDict
import sys

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
	if len([n for n in dicc[k] if n.startswith(">")])>= 2:
		w = open("mmseqs_all_vs_all/clusters/Clu_"+str(clu)+".fa","w+")
		for n in dicc[k]:
			w.write(n)
		w.close()
