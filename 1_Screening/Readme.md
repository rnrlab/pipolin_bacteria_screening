## 1 - Screening

The two scripts "screening_v13.py" and "genome_downloading.py" are responsible for the cyclic concurrent download and screening of pipolins in genome assemblies (in batches of 50). Screening starts by running only "screening_v13.py", which will use "genome_downloading.py" later. May need code adjustments if NCBI Datasets version is distinct than 13.43.2.

After the revisions, we also included the results of detecting alternative versions of pipolins (i.e. piPolB-less pipolins) with MacSyFinder and the results of the screening in Archea.
