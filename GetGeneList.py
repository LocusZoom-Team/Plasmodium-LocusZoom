import pandas as pd
import re

records = []

geneTypes = {"protein_coding_gene" : 1, "ncRNA_gene" : 0, "pseudogene":0}

with open("Pfalciparum3D7.gff") as f:
    for line in f:
        if line.startswith("#"):
            continue
        section = line.strip().split("\t")
        if section[2] not in geneTypes:
            continue

        chromosome = section[0]
        startPosition = int(section[3])
        endPosition = int(section[4])
        attributes = section[8]

        chromosomeMatch = re.search(r'Pf3D7_(\d+)_v3', chromosome)

        if not chromosomeMatch:
            continue

        chromosomeNumber = int(chromosomeMatch.group(1))

        attributesDict = {}

        for attribute in attributes.split(";"):
            key, value = attribute.split("=",1)
            attributesDict[key.strip()] = value.strip()
        
        if attributesDict.get("Name"):
            gene_id = attributesDict.get("Name")
        else:
            gene_id = attributesDict.get("ID")

        records.append({"Gene": gene_id, "Chrom": chromosomeNumber, "Start": startPosition, "End": endPosition, "Coding": geneTypes[section[2]]})

df = pd.DataFrame(records)
df = df.sort_values("End").drop_duplicates(subset="Gene", keep="last")
df = df.sort_values(["Chrom", "Start"]).reset_index(drop=True)

df.to_csv("Pfalciparum_gene_list.tsv", sep="\t", index=False)