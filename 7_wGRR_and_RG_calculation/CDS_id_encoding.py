import random
import string
import sys

def generate_unique_identifiers(num_names):
    identifiers = set()  # Using a set to ensure uniqueness
    while len(identifiers) < num_names:
        identifier = ''.join(random.choices(string.ascii_letters + string.digits, k=3))
        if identifier[0:2] != "uc": #uc identifiers give errors in mmseqs
            identifiers.add(identifier)
    return list(identifiers)

unique_identifiers = generate_unique_identifiers(100000)

mge_to_code = {}
coded_mge = ""
n = 0
with open(sys.argv[1], "r") as f:
    for line in f:
        if ">" in line:
            cds_id = line.split()[0].replace(">","")
            mge_id = "_".join(cds_id.split("_")[:-1])
            if mge_id not in list(mge_to_code.keys()):
                n += 1
                mge_to_code[mge_id] = unique_identifiers[n]
            coded_mge += line.replace(mge_id, mge_to_code[mge_id])
        else:
            coded_mge += line

with open(sys.argv[1].replace(".faa","_encoded.faa"), "w") as f:
    f.write(coded_mge)

with open("mge_codes.tsv", "w") as f:
    f.write("MGE_id\tcode\n")
    for mge in mge_to_code:
        f.write(mge +"\t" + mge_to_code[mge] + "\n")

print("Encoding finished")