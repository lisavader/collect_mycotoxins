import argparse
import logging
import glob
import os
import json
import csv

from mibig_classes import Compound

mycotoxin_names = ["aflatoxin","citreoviridin","citrinin","fumonisin","fusarin","nitropropionic acid","ochratoxin","patulin",
                   "sterigmatocystin","trichothecene","t-2","ht-2","neosolaniol","diacetoxyscirpenol","deoxynivalenol","nivalenol",
                   "fusarenon","zearalenone"]

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
                summary=[compound.mibig_accession,compound.organism,compound.index,compound.name]
                mycotoxin_mibig_data.append(summary)
    return mycotoxin_mibig_data

def main(mibig_path,output_path):
    mycotoxin_mibig_data = get_mycotoxin_mibig_data(mibig_path,mycotoxin_names)
    with open(output_path,"w") as file:
        writer = csv.writer(file)
        #Write a header
        header = ["MIBiG Accession","Organism","Compound Index","Compound Name"]
        writer.writerow(header)
        writer.writerows(mycotoxin_mibig_data)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("mibig_path", type=str, help="Path to a MIBiG database of annotations in .json format")
    parser.add_argument("output_path", type=str, help="Path to write output to in .csv format")
    parser.add_argument("-v","--verbose",help="Enable verbose mode",action="store_const",const=logging.INFO,dest="loglevel")
    args = parser.parse_args()
    #Set logging level
    logging.basicConfig(level=args.loglevel)
    #Run the main script
    main(args.mibig_path,args.output_path)