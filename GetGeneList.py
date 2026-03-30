import pandas as pd #Import this for the final file generation
import re #Import package to search for chromosome name

records = [] #Created an empty list to hold things in

geneTypes = {"protein_coding_gene" : "Coding", "ncRNA_gene" : "Non-Coding", "pseudogene": "Non-coding"} #Gene types to grab our total genes, number means coding or not

with open("Pfalciparum3D7.gff") as f: #Open file
    for line in f: #go line by line
        if line.startswith("#"): #Anything with a # is a commment, so skip
            continue #go to next line
        section = line.strip().split("\t") #split section based on tab since tab delimitedp
        if section[2] not in geneTypes:
            continue

        chromosome = section[0]
        startPosition = int(section[3])
        endPosition = int(section[4])
        attributes = section[8]

        chromosomeMatch = re.search(r'Pf3D7_\d+_v3', chromosome)

        if not chromosomeMatch:
            continue

        

        attributesDict = {}

        for attribute in attributes.split(";"):
            key, value = attribute.split("=",1)
            attributesDict[key.strip()] = value.strip()
        
        if attributesDict.get("Name"):
            gene_id = attributesDict.get("Name")
        else:
            gene_id = attributesDict.get("ID")

        records.append({"Gene": gene_id, "Chrom": chromosome, "Start": startPosition, "End": endPosition, "Coding": geneTypes[section[2]]})

df = pd.DataFrame(records)
df = df.sort_values("End").drop_duplicates(subset="Gene", keep="last")
df = df.sort_values(["Chrom", "Start"]).reset_index(drop=True)

df.to_csv("Pfalciparum_gene_list.tsv", sep="\t", index=False)