import argparse
import logging
import pandas as pd
import json
import glob
import os
import csv

from mibig_classes import DatabaseID, Compound

def prefilter_mibig_entries(mibig_path,mycotoxin_formulas):
    """ Find entries that contain a compound with a molecular formula matching a mycotoxin, or for which the molecular formula is unknown."""
    #List all .json files in MIBiG database
    print("Finding entries in database...")
    mibig_entries=set(glob.glob(mibig_path+"/*.json"))
    print("Total nr. of MIBiG entries found: "+str(len(mibig_entries)))

    print("Prefiltering entries based on molecular formula...")
    filtered_entries=set()
    missing_formula_count=0

    #Load json data for each entry
    for entry in mibig_entries:
        with open(entry, 'r') as json_file:
            bgc_data = json.load(json_file)
        #Find number of associated compounds
        compound_count=len(bgc_data["cluster"]["compounds"])
        #If the compounds molecular formula matches a mycotoxin molecular formula, add entry to filtered list
        for i in range(0,compound_count):
            try:
                molecular_formula=bgc_data["cluster"]["compounds"][i]["molecular_formula"]
                if molecular_formula in mycotoxin_formulas:
                    filtered_entries.add(entry)
                    break
            except KeyError:
                #Add to filtered list if molecular formula is unknown
                filtered_entries.add(entry)
                missing_formula_count+=1
                break

    print("Nr. of entries kept after filtering: "+str(len(filtered_entries)))
    print("Nr. of entries kept due to missing molecular formula: "+str(missing_formula_count))
    return(filtered_entries)

def get_mycotoxin_mibig_data(filtered_entries,mycotoxin_inchikeys):
    """ Find compounds in the MIBiG database that are considered mycotoxins, according to an inchikey match.
        Return list with compound properties (MIBiG accession, index, name and inchikey)."""
    mycotoxin_mibig_data=[]
    missing_id_count=0

    print("Finding mycotoxin entries based on inchikey...")
    print("Searching databases for inchikeys...")
    for entry in filtered_entries:
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
            try:
                #Retrieve a list of all database ids
                id_list=bgc_data["cluster"]["compounds"][i]["database_id"]
                #Store ids as DatabaseID objects, and add to database_ids attribute of the Compound object
                for database_id in id_list:
                    database,id=database_id.split(":")
                    compound.database_ids.append(DatabaseID(database,id))
                #If the compound is a mycotoxin, store compound properties (MIBiG ID, organism, index, name and inchikey)
                compound.mycotoxin_from_inchikey(mycotoxin_inchikeys)
                if compound.mycotoxin == True:
                    mycotoxin_mibig_data.append(compound.summary())

            #Keep track of compounds without any database id
            except KeyError:
                missing_id_count+=1
                pass

    if missing_id_count != 0:
        logging.warning(str(missing_id_count)+" compounds with missing database id!")

    #Return data for all mycotoxin compounds
    return(mycotoxin_mibig_data)


def main(comptox_path,mibig_path):
    #Load CompTox table with mycotoxin data
    print("Loading mycotoxin data from CompTox table...")
    comptox_mycotoxins=pd.read_csv(comptox_path)
    #Extract molecular formulas and inchikeys of mycotoxins
    mycotoxin_formulas=set(comptox_mycotoxins["MOLECULAR FORMULA"])
    mycotoxin_inchikeys=set(comptox_mycotoxins["INCHIKEY"])
    #Prefilter MiBIG entries based on molecular formula
    filtered_entries=prefilter_mibig_entries(mibig_path,mycotoxin_formulas)
    #Find which MiBIG entries are linked to mycotoxin production
    mycotoxin_mibig_data=get_mycotoxin_mibig_data(filtered_entries,mycotoxin_inchikeys)
    #Save as csv
    with open("mycotoxin_mibig_data.csv","w") as file:
        writer = csv.writer(file)
        #Write a header
        header = ["MIBiG Accession","Organism","Compound Index","Compound Name","Inchikey"]
        writer.writerow(header)
        writer.writerows(mycotoxin_mibig_data)


if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("comptox_path", type=str, help=("Path to a .csv file containing a"
                                                        " CompTox list of mycotoxins"
                                                        " (url:https://comptox.epa.gov/dashboard/chemical-lists/MYCOTOX2)"))
    parser.add_argument("mibig_path", type=str, help="Path to a MIBiG database of annotations in .json format")
    parser.add_argument("-v","--verbose",help="Enable verbose mode",action="store_const",const=logging.INFO,dest="loglevel")
    args = parser.parse_args()
    #Set logging level
    logging.basicConfig(level=args.loglevel)
    #Run the main script
    main(args.comptox_path,args.mibig_path)
