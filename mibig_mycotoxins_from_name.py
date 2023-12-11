import argparse
import logging
import glob
import os
import json
import csv
import re

from mibig_classes import Compound

mycotoxins_tier1 = ["aflatoxin","citreoviridin","citrinin","fumonisin","fusarin","nitropropionic acid","ochratoxin","patulin",
                   "sterigmatocystin","t-2","neosolaniol","diacetoxyscirpenol","deoxynivalenol","nivalenol",
                   "fusarenon","zearalenone"]
mycotoxins_tier2 = ["aflatrem","alternariol","altertoxin","beauvericin","butenolide","cladofulvin","culmorin","cyclopiazonic acid",
                    "cyclosporin","cytochalasin","chaetoglobosin","phomacin","flavichalasin","diplodiatoxin","diplonine","dipmatol","emodin",
                    "enniatin","ergotamine","gliotoxin","fumagillin","fumitremorgin","fusaproliferin","fusaric acid","janthitrem","lolitrem",
                    "luteoskyrin","moniliformin","mycophenolic acid","paspalitrem","paxilline","penicillic acid","penitrem","pr-toxin",
                    "roquefortine","secalonic acid","tenuazonic acid","trypacidin","verrucosidin","verruculogen"]
mycotoxins_all = mycotoxins_tier1 + mycotoxins_tier2

def get_mycotoxin_mibig_data(mibig_path,mycotoxin_names):
    mycotoxin_mibig_data=[]

    mibig_entries=set(glob.glob(mibig_path+"/*.json"))
    for entry in mibig_entries:
        #Get name of mibig accession from file path
        filename=os.path.basename(entry)
        mibig_accession=os.path.splitext(filename)[0]

        #Load json data
        with open(entry, 'r') as json_file:
            bgc_data = json.load(json_file)

        organism_name=bgc_data["cluster"]["organism_name"]
        #Find all associated compounds
        compound_count=len(bgc_data["cluster"]["compounds"])
        for i in range(0,compound_count):
            compound_name=bgc_data["cluster"]["compounds"][i]["compound"]
            #Store compound info in a new Compound object
            compound=Compound(mibig_accession,i,compound_name,organism_name)
            #If the compound is a mycotoxin, store compound properties (MIBiG ID, organism, index, name and inchikey)
            compound.mycotoxin_from_name(mycotoxin_names, exact_match = False)
            if compound.mycotoxin == True:
                summary=[compound.query,compound.name,compound.mibig_accession,compound.organism,compound.index]
                mycotoxin_mibig_data.append(summary)
    return mycotoxin_mibig_data

def add_missing_compounds(mycotoxin_mibig_data):
    for compound_query in mycotoxins_all:
        found = False
        for row in mycotoxin_mibig_data:
            compound_in_mibig=row[0]
            if re.search(compound_query,compound_in_mibig,re.IGNORECASE):
                found = True
                continue
        if found == False:
            summary=[compound_query,"NA","NA","NA","NA"]
            mycotoxin_mibig_data.append(summary)
    return mycotoxin_mibig_data

def main(mibig_path,output):
    mycotoxin_mibig_data = get_mycotoxin_mibig_data(mibig_path,mycotoxins_all)
    mycotoxin_mibig_data = add_missing_compounds(mycotoxin_mibig_data)
    with open(output,"w") as bgc_file:
        writer = csv.writer(bgc_file)
        #Write a header
        header = ["query","match","mibig_accession","organism","compound_index"]
        writer.writerow(header)
        writer.writerows(mycotoxin_mibig_data)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("mibig_path", type=str, help="Path to a MIBiG database of annotations in .json format")
    parser.add_argument("output", type=str, help="Path to write bgc output to in .csv format")
    parser.add_argument("-v","--verbose",help="Enable verbose mode",action="store_const",const=logging.INFO,dest="loglevel")
    args = parser.parse_args()
    #Set logging level
    logging.basicConfig(level=args.loglevel)
    #Run the main script
    main(args.mibig_path,args.output)