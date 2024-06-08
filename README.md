# Identification and characterization of pipolins in Bacteria

**Pipolins** constitute a new group of mobile genetic elements (MGEs) that encode a primer-independent B-family DNA polymerase (piPolB). Previously, pipolins had been identified as integrative MGEs in diverse bacteria and as circular plasmids in mitochondria and a few gram positive species. Detailed analysis of *E. coli* pipolins revealed pipolins are present in diverse phylogroups, encode a diverse repertoire of DNA metabolism genes, and show evidence of recent horizontal transfer between different *E. coli* strains. 
> [**Redrejo-Rodríguez, M., *et al.*** Primer-independent DNA synthesis 
>by a family B DNA polymerase from self-replicating Mobile genetic elements. 
>*Cell reports*, 2017](https://doi.org/10.1016/j.celrep.2017.10.039)
> [**Chuprikova, L., Mateo-Cáceres, V., de Toro, M. & Redrejo-Rodríguez, M.** ExplorePipolin: reconstruction and annotation of piPolB-encoding bacterial mobile elements from draft genomes. *Bioinformatics Advances* 2022.](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbac056/6659502)

In this repository, we have made available the code developed during our latest work on pipolins. We carried out a pipolin screening of the Assembly (NCBI) database. The analysis of the structure of pipolins revealed that they are commonly integrative elements, usually flanked by direct repeats, sharing known mobile elements integration hotspots (*e.g.* tRNA genes). Remarkably, integrase dynamics correlates with alternative integration spots and enables diverse lifestyles, ranging from integrative to mobilizable plasmid pipolins. Pipolins harbor a minimal core and a large cargo module enriched for defense factors, which are actively exchanged with other mobile elements. These findings indicate pipolins act as orthogonal reservoirs of defense genes that play a key role in the exchange mechanisms for defense genes in bacterial populations.

Pre-print of this work is available in biorxiv:
> [**Mateo-Cáceres, V., Redrejo-Rodríguez, M.** Pipolins are bimodular platforms that maintain a reservoir of defense systems exchangeable with various bacterial genetic mobile elements. 
>*bioxriv*, 2024](https://doi.org/10.1101/2024.05.22.595293)

Each folder in this repository contains the code used in each project sub-task:
1. **Screening**: Download genomes from Assembly using Datasets (NCBI) and pipolin detection with ExplorePipolin
1. **Pipolin metadata**: Parse genome metadata and ExplorePipolin results to calculate and plot statistics.
1. **piPolB phylog**: Infer and plot piPolB phylogeny
1. **pipolin CDS clustering**: Cluster proteins encoded in pipolins and reannotation. 
1. **candidate recombinases**: Create cluster presence-abscence matrix and find candidate pipolin recombinases.
1. **MGE extTools annot**: Parse information of different annotation tools used on pipolins and other MGEs (plasmids, phages, ciMGEs). Plot annotation resuts.
1. **wGRR and RG calculation**: Calculate wGRR (wighted gene repertoire relatedness) among pipolins and other annotated MGEs to detect recent events of gene exchange (aka recombining genes or RGs).

*Discalimer:
- Work in progress.
- Some scripts may require additional software/packages, files or databases not included in this repository. Since some input and output files surpass Github size, the correct execution of several scripts uploaded here is greatly impaired. Feel free to request the authors any missing file that might interest the reader and do not hesitate to ask for help if you encounter any difficulies while using the code provided in this repository*
