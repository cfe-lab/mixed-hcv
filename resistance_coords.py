"""
Massages the resistance AA positions so that it is compatible with HCVDeli.py

We want to go from cols:  PROT	GENO	WT	POS	AA	DRUG
to cols:  Gene,StartNuc_0based,AfterEndNuc_0based.

Only care about Genotype 1A - NS5a right now.
"""
import csv
import os

RESISTANCE_CSV = os.path.dirname(os.path.realpath(__file__)) + "/data/HCV_ResMutList.csv"
RESISTANCE_TARGET_CSV = os.path.dirname(os.path.realpath(__file__)) + "/data/HCV_ResMutList_TargetNucCoord.csv"
DESIRED_SUBTYPES = ["1A"]
DESIRED_PROT = ["NS5a"]

seen_row = set()
with open(RESISTANCE_CSV, 'rU') as fh_in, open(RESISTANCE_TARGET_CSV, 'w') as fh_out:
    reader = csv.DictReader(fh_in)
    writer = csv.DictWriter(fh_out, fieldnames=["Gene", "StartNuc_0based", "AfterEndNuc_0based"])
    writer.writeheader()
    for line in reader:
        subtype = line["GENO"]
        protein = line["PROT"]
        if subtype not in DESIRED_SUBTYPES or protein not in DESIRED_PROT:
            continue

        aa_pos_1based = int(line["POS"])
        nuc_pos_0based =  (aa_pos_1based - 1) * 3
        outrow = dict()
        outrow["Gene"] = protein
        outrow["StartNuc_0based"] = nuc_pos_0based
        outrow["AfterEndNuc_0based"] = nuc_pos_0based + 3

        if (protein, nuc_pos_0based, nuc_pos_0based + 3) not in seen_row:
            seen_row.add((protein, nuc_pos_0based, nuc_pos_0based + 3))
            writer.writerow(outrow)



