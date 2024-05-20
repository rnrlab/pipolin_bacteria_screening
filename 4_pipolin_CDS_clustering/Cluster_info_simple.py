new_cluster_information = ""

#transposases
Yrec_terms = ["PHAGE_INT", "PHAGE_INTEGRASE"]
Srec_terms = ["RESOLVASE,", "RESOLVASE'", ",RESOLVASE", "'RESOLVASE"]
transposase_terms = ["TRANSPOSASE", "TNP_", "'RVE'", "RVE'", ",RVE", "RVE,", "TNPB_"]

#AN metab
NA_mod = ["DNA METHYLASE", "METHYLASE_S", "_MTASE", "DNA METHYLTRANSFERASE", "ECORI_METHYLASE", "UDG", "DNA-SULFUR MODIFICATION",
          "RNA METHYLTRANSFERASE"]
NA_cleav = ["RESTRICTION ENZYME","RESTRICTION-MODIFICATION", "RESTRICTION MODIFICATION", "HSDR", "REASE_AHJR", "MRR_CAT", "RESTRICTION SYSTEM COMPONENT", 
            "RESTRICTION ENDONUCLEASE", "HEPN", "HNH_", "RAMA", "PDDEX", "PD-(D/E)XK", "NERD", "'PIN_", "PIN DOMAIN", "ENDORIBONUCLEASE",
            "MISMATCH ENDONUCLEASE", ]

nuc_metab = ["DEOXYCYTIDYLATE DEAMINASE", "NUCLEOTIDE PYROPHOSPHOHYDROLASE", "EAL", "PRA-CH", "GGDEF", "DEOXYNUCLEOSIDE KINASE",
             "NTP_TRANSFERASE", "GUANYLATE_CYC"]

other_hydrolases = ["ABHYDROLASE", "PHOSPHOFYDROLASE", "CN_HYDROLASE", "PHOSPHOESTERASE", "PHOSPHODIESTERASE", "HYDROLASE", "PHOSPHATASE"]


other_dna_bind = ["HELIX-TURN-HELIX", "BETR", "WYL", "ARM DNA-BINDING DOMAIN", "NUCLEIC-ACID-BINDING", "DNA-BINDING", "DNA BINDING",
                  "CSD", "H-NS", "NA37", "TRANSCRIPTION REGULATOR", "TRANSCRIPTIONAL REGULATOR", "SOPB_HTH", "DNA RECOMBINATION",
                  "MOBC", "MERR", "SIR2", "PROPHAGE PROTEIN COM", "NUCLEIC-ACID-BINDING"]

#membrane and atpases
AAA_Helicase = ["AAA", "HELICASE"]
membrane = ["ABC_TRAN", "PERMEASE", "CHANNEL", "MEMBRANE", "PUMP", "BAND_7", "ABC-3C", "EXPORTER", "ANTIPORTER", "SECRETION SYSTEM",
                "STAS", "TRANSPORTER", "FIMBRIAL", "PILUS", "USHER", "PAPD", "PAPC"]

#Pass info
pass_clusters = []
with open("../Cluster_association/all_context_genes/list_MMSeqs_cluster_pass.txt", "r") as f:
    for line in f:
        if "Cluster_" in line:
            pass_clusters.append(line.split("\t")[0])

n = 0

with open("Cluster_information.txt", "r") as f:
    for line in f:
        if "G_" in line:
            n += 1
            fields = line.split("\t")
            #print(line, n, fields[1])
            cluster = fields[1]

            pfam = fields[7]
            piPolB = fields[9]
            hhblits = fields[14]
            if piPolB != "0":
                new_cluster_information += line.replace("\n","\tpiPolB\tff2600\n")
            else:

                Yrec_check = False
                for term in Yrec_terms:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tYrec\t6b3700\n")
                        Yrec_check = True
                        break    
                if Yrec_check:
                    continue
                Srec_check = False
                for term in Srec_terms:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tSrec\t967d1a\n")
                        Srec_check = True
                        break    
                if Srec_check:
                    continue
                transposase_check = False
                for term in transposase_terms:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tTransposase\te3ad74\n")
                        transposase_check = True
                        break    
                if transposase_check:
                    continue


                if "RELAXASE" in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tRelaxase\t74e3b3\n")
                        continue

                if "PEPTIDASE" in hhblits.upper() or "PROTEASE" in hhblits.upper() or "PROTEINASE" in hhblits.upper():
                        new_cluster_information += line.replace("\n","\tPeptidase\tfa78d1\n")
                        continue


                if "EXCISIONASE" in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tExcisionase\td578fa\n")
                        continue

                #DNA mod and cut
                NA_mod_check = False
                for term in NA_mod:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tNA Modification\tfffd99\n")
                        NA_mod_check = True
                        break    
                if NA_mod_check:
                    continue    
                
                NA_cleav_check = False
                for term in NA_cleav:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tNA Cleavage\tf5f118\n")
                        NA_cleav_check = True
                        break    
                if NA_cleav_check:
                    continue    
                
               
                nuc_metab_check = False
                for term in nuc_metab:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tNucleotide metabolism\t1780d1\n")
                        nuc_metab_check = True
                        break    
                if nuc_metab_check:
                    continue

                other_hydrolases_check = False
                for term in other_hydrolases:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tOther hydrolases\tf5e9dc\n")
                        other_hydrolases_check = True
                        break    
                if other_hydrolases_check:
                    continue

                #Membrane
                membrane_check = False
                for term in membrane:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tMembrane\t3b8c3e\n")
                        membrane_check = True
                        break    
                if membrane_check:
                    continue 

                AAA_Helicase_check = False
                for term in AAA_Helicase:
                    if term in pfam.upper()+hhblits.upper():
                        new_cluster_information += line.replace("\n","\tHelicase-AAA\tffb042\n")
                        AAA_Helicase_check = True
                        break    
                if AAA_Helicase_check:
                    continue   
                

                

                #Other DNA bind 
                other_dna_bind_check = False
                for term in other_dna_bind:
                    if term in hhblits.upper():
                        new_cluster_information += line.replace("\n","\tOther DNA binding\ta4ccde\n")
                        other_dna_bind_check = True
                        break    
                if other_dna_bind_check:
                    continue  



                #This is for testing if only DUF/UPF found
                all_DUF_UPF = True
                for hhblits_match in hhblits.split(","):
                     if "Uncharacterised protein family" not in hhblits_match and "Domain of unknown function" not in hhblits_match and\
                        "Protein of unknown function" not in hhblits_match and "Family of unknown function" not in hhblits_match and\
                        "Bacterial protein of unknown funciton" not in hhblits_match and "Uncharacterized protein conserved in bacteria" not in hhblits_match:
                          all_DUF_UPF = False
                if all_DUF_UPF:
                        new_cluster_information += line.replace("\n","\tDUF and UPF\tc4c4c4\n")
                        continue                


                if hhblits == "[]\n" and pfam[0:5] == "{'-':" and ", '" not in pfam:
                     new_cluster_information += line.replace("\n","\tNo PFAM\t4d4d4d\n")
                     continue
            
                new_cluster_information += line.replace("\n","\tOther\tffffff\n")



        else:
            new_cluster_information += line.replace("\n","\tCustom_category\tColor_code\n") 

with open("Cluster_info_simplified.txt", "w") as f:
     f.write(new_cluster_information)



