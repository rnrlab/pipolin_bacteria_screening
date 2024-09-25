#Note: Run in the same folder than  genome_downloading.py
'''
Part 1. Preparations must be made
'''
import multiprocessing
import os
import sys
import subprocess
import time
from datetime import datetime
from Bio import SeqIO

start_time = time.time()
current_day_h = str(datetime.now()).replace(" ","_").replace(":",".")[:-7]
term = input("Enter taxon name or id: ").lower()
filename = term.replace(" ","_")+".zip"
out_prefix = input("Input prefix of results directory: ")
os.mkdir("results_"+out_prefix+"_"+current_day_h)
os.chdir("results_"+out_prefix+"_"+current_day_h)


'''
Part 2. NCBI Datasets download proteins
'''
start_time_datasets = time.time()

#Database selection
database_question = input("Do you want to download genomes from GenBank, RefSeq or Both (GenBank default)? ")
if database_question == "Both":
    print("Downloading genomes from both GenBank and RefSeq.")
    database_filter = "all"
elif database_question == "RefSeq":
    print("Downloading genomes from RefSeq.")
    database_filter = "RefSeq"
else:
    print("Downloading genomes from GenBank.")
    database_filter = "GenBank"

assembly_level_question = input("Limit to genomes at one or more assembly levels (comma-separated: chromosome,complete,contig,scaffold) -default all- :")
assembly_level_param = "--assembly-level "+assembly_level_question
if assembly_level_question == "":
    assembly_level_param = ""

assembly_date_question = input("Limit to genomes before date (YYYY-MM-DD) -default today- :")
assembly_date_param = "--released-before "+assembly_date_question
if assembly_date_question == "":
    assembly_date_param = ""

final_datasets_command = " ".join(["datasets_v16.22.1", "download", "genome", 
                                   "taxon", term, 
                                   "--filename", filename, 
                                   "--exclude-atypical", 
                                   "--assembly-source",  database_filter, 
                                   assembly_level_param,
                                   assembly_date_param,
                                   "--include protein,gff3", 
                                   "--assembly-version latest", 
                                   "--dehydrated"])

print("Final command: ", final_datasets_command)

subprocess.run(final_datasets_command, shell=True) 

#pre: datasets_v16.22.1 download genome taxon 192989 --filename test_nano.zip --exclude-atypical --assembly-source RefSeq --assembly-level chromosome,complete --assembly-version latest --preview
#--include genome, , excluded for now
#alter --released-before 2022-11-11

#Unzip dehydrated zip
fetch_folder_name = filename[:-4]+"_dataset"
subprocess.run(["unzip", filename, "-d", fetch_folder_name])
subprocess.run("datasets rehydrate --directory "+fetch_folder_name+" --max-workers 30", shell=True)


