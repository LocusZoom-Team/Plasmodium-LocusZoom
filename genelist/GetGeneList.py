import pandas as pd #Import this for the final file generation
import re #Import package to search for chromosome name

records = [] #Created an empty list to hold things in

geneTypes = {"protein_coding_gene" : "Coding", "ncRNA_gene" : "Non-Coding", "pseudogene": "Non-Coding"} #Gene types to grab our total genes, number means coding or not

with open("Pfalciparum3D7.gff") as f: #Open file
    for line in f: #go line by line
        if line.startswith("#"): #Anything with a # is a commment, so skip
            continue #go to next line
        section = line.strip().split("\t") #split section based on tab since tab delimitedp
        if section[2] not in geneTypes: #Filter out anything that's not a gene
            continue #skip ahead to the next line

        chromosome = section[0] #Save chromosome 
        startPosition = int(section[3]) #Save start position, and convert from int to string
        endPosition = int(section[4]) #Save end position, and convert from int to string
        attributes = section[8] #Save attributes, this will have our gene name (if present)

        chromosomeMatch = re.search(r'Pf3D7_\d+_v3', chromosome) #Use re package to filter our apicoplast and mitochondrial chromosomes out

        if not chromosomeMatch: #If apicoplast or mitochrondrial, skip
            continue #continue to next line

        

        attributesDict = {} #Create attribute dictionary to store attributes for gene name

        for attribute in attributes.split(";"): #Split the attribute based on semicolon, that's how attribute is split
            key, value = attribute.split("=",1) #Then split it again based on equal sign, use tuple unpacking
            attributesDict[key.strip()] = value.strip() #add name or ID for key in dictionary, and value, strip whitespace
        
        if attributesDict.get("Name"): #See if name is in attributes dict
            gene_id = attributesDict.get("Name") #Get value for name
        else:
            gene_id = attributesDict.get("ID") #Other get value for ID

        records.append({"Gene": gene_id, "Chrom": chromosome, "Start": startPosition, "End": endPosition, "Coding": geneTypes[section[2]]}) #Add to our file

df = pd.DataFrame(records) #Change file to dataframe to sort lines
df = df.sort_values("End").drop_duplicates(subset="Gene", keep="last") #get rid of duplicates, keep the last one
df = df.sort_values(["Chrom", "Start"]).reset_index(drop=True) #Sort df based on chrom and start

df.to_csv("Pfalciparum_gene_list.txt", sep="\t", index=False) #write to file and we are done
