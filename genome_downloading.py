### Imports
import os
import sys
import subprocess
import time
from Bio import SeqIO



### Parameters and variables
split_n = sys.argv[1]
fetch_path = sys.argv[2]
fetch_folder_name = sys.argv[3]
file = sys.argv[4]
n_tries = int(sys.argv[5])

#for the AN:ID dictionary aka conversion table
if "conversion_table.txt" in os.listdir():
    with open("conversion_table.txt", "r") as f:
        content = f.readlines()
        n_identif = len(content)
        conversion_table = "".join(content)
else:
    n_identif = 0
    conversion_table = ""

#Stores genomes (fetch.txt format) that couldn't be downloaded after 5 tries
if "not_downloaded.txt" in os.listdir():
    with open("not_downloaded.txt", "r") as f:
        content = f.readlines()
        download_error = "".join(content)
else:
    download_error = ""



### Functions
def rehydrate_genomes(split_n, n_identif, fetch_path, fetch_folder_name, file, n_tries, conversion_table, download_error):
    ### PART 0: CREATE RECOVERY FILE ###
    old_n_recovery = int(file.replace(".txt","").split("_")[-1])-1
    old_recovery_file = "_".join(file.replace(".txt","").split("_")[:-1])+"_"+str(old_n_recovery)+".txt.recovery.txt"
    if old_recovery_file in os.listdir(fetch_path):
        subprocess.run(["rm", fetch_path+old_recovery_file])
        subprocess.run(["cp", fetch_path+file, fetch_path+file+".recovery.txt"])
    else:
        subprocess.run(["cp", fetch_path+file, fetch_path+file+".recovery.txt"]) #if screening ends

    ### PART A: REHYDRATE ###
    subprocess.run(["mv", fetch_path+file, fetch_path+"fetch.txt"])
    subprocess.run("datasets rehydrate --directory "+fetch_folder_name+" --max-workers 30", shell=True, stdout=subprocess.DEVNULL)
    print("Split download completed. Verifying downloaded files. ")

    ### PART B: VERIFY AND RETRY DOWNLOADS ###
    vh_return = verify_hydration(n_tries, fetch_path, n_identif, download_error)

    ### PART C: GENOME CONTIG RENAMING ###
    gp_return = genome_processing(fetch_path, conversion_table, n_identif)

    ### PART D: UPDATE FILES ###
    with open("conversion_table.txt", "w") as f_ct:
        f_ct.write(gp_return)

    with open("not_downloaded.txt", "w") as f_nd:
        f_nd.write(vh_return)

    ### PART E: FINAL ###
    #Delete fetch.txt for the new splits to come
    subprocess.run("rm "+fetch_path+"fetch.txt", shell = True)
    print("Rehydration and renaming complete.")

def verify_hydration(n_tries, fetch_path, n_identif, download_error):
#Note: would it be better to retry downloads in each split or each X splits? (just a few files are expected to fail)
    print("Files that couldn't be downloaded will be added to a retry list.")
    #Number of retries
    for i in range(n_tries): 
        #Check if downloaded genomes exist
        retry_list = ""
        incomplete = False
        with open(fetch_path+"fetch.txt", "r") as file:
            for line in file:
                if not os.path.exists(fetch_path+line.strip("\n").split("\t")[-1]):
                    retry_list += line
                    incomplete = True
                    print("File "+line.strip("\n").split("\t")[-1].split("/")[-1]+" couldn't be downloaded and has been added to the retry list.")
        #Now downloaded genomes are retried
        if incomplete:
            with open(fetch_path+"fetch.txt", "w") as file:
                file.write(retry_list)
            print("Retrying download (try "+str(i+1)+" out of 5):")
            subprocess.run("datasets rehydrate --directory "+fetch_folder_name+" --max-workers 30", shell=True, stdout=subprocess.DEVNULL)
        #If there is still genomes that couldn't be downloaded the user is warned
    if incomplete:
        download_error += retry_list
    return download_error

def genome_processing(fetch_path, conversion_table, n_identif):
    #Concatenate and modify file name and contig name
    for genome_folder in os.listdir(fetch_path+"data/"):
        if os.path.isdir(fetch_path+"data/"+genome_folder):
            n_identif += 1
            fasta_files = [fa for fa in os.listdir(fetch_path+"data/"+genome_folder) if ".fna" in fa]
            corrected_file = "tmp_genomes/G_"+str(n_identif)+".fna"
            subprocess.run("cat "+fetch_path+"data/"+genome_folder+"/*.fna > "+fetch_path+"data/"+genome_folder+"/G_"+str(n_identif)+".fna", shell=True)
            with open(fetch_path+"data/"+genome_folder+"/G_"+str(n_identif)+".fna", "r") as original, open(corrected_file, 'w') as corrected:
                record_num = 0  
                plasmid_num = 0
                records = SeqIO.parse(fetch_path+"data/"+genome_folder+"/G_"+str(n_identif)+".fna", 'fasta')
                for record in records:
                    if "plasmid" in str(record.description):
                        plasmid_num += 1
                        record.id = "G_"+str(n_identif)+'_p_'+str(plasmid_num)
                    else:
                        record_num += 1 
                        record.id = "G_"+str(n_identif)+'_r_'+str(record_num)
                    SeqIO.write(record, corrected, 'fasta')

            #Record in conversion_table.txt 
            conversion_table += str(n_identif)+"\t"+"G_"+str(n_identif)+".fna"+"\t"+genome_folder+"\t"+";".join(fasta_files)+"\n"

            #Remove original genomes
            subprocess.run(["rm", "-r", fetch_path+"data/"+genome_folder])

    #Stats
    print("Number of sequences: ", n_identif)
    print("Number of not downloaded: ")
    subprocess.run("wc -l not_downloaded.txt", shell=True)
    
    #Get back conversion table
    return conversion_table 

### Run
start_genome_set=time.time()
rehydrate_genomes(split_n, n_identif, fetch_path, fetch_folder_name, file, n_tries, conversion_table, download_error)
print("Rehydration process finished in:", str(time.time() - start_genome_set), "seconds.")
