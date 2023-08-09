import argparse
import logging
import pandas as pd
import json
import glob
import os
import requests
import time
import chemspipy
import csv
from bs4 import BeautifulSoup

#For access to ChemSpider an API key is needed (see: https://developer.rsc.org/get-started)
chemspider_api_key = "" #<-- input your API key here

class DatabaseID:
    """ An identifier within a database of chemical structures """
    def __init__(self,database,id):
        self.database: str = database
        self.id: str = id
        self.url: str = ""

    def lookup_url(self) -> requests.Response:
        """ Send HTTP request and handle HTTP errors"""
        for retries in range(3):
            try:
                response=requests.get(self.url)
                response.raise_for_status()
            except requests.HTTPError:
                #Npatlas throws a 429 error once we've made too many requests, it seems waiting for a minute is the only solution
                if response.status_code == 429:
                    print(self.database+":"+self.id+": HTTPError "+str(response.status_code)+": Sleeping for 60 seconds...")
                    time.sleep(60)
                    continue
                #Other http errors for which we can retry the connection
                elif response.status_code in (408,503):
                    print(self.database+":"+self.id+": HTTPError "+str(response.status_code)+": Sleeping for 1 second...")
                    time.sleep(1)
                    continue
            break
        return(response)

    def get_inchikey(self) -> str:
        """ Search database to retrieve inchikey """
        inchikey = None

        try:
            if self.database == "chembl":
                self.url="https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_chembl_id__exact="+self.id+"&format=json"
                response=self.lookup_url()
                inchikey=response.json()["molecules"][0]["molecule_structures"]["standard_inchi_key"]

            elif self.database == "chebi":
                self.url="https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId="+self.id
                response=self.lookup_url()
                soup = BeautifulSoup(response.text, 'xml')
                inchikey=soup.find("inchiKey").text

            elif self.database == "npatlas":
                self.url="https://www.npatlas.org/api/v1/compound/"+self.id
                response=self.lookup_url()
                inchikey=response.json().get("inchikey")

            elif self.database == "chemspider":
                global chemspider_api_key
                if not chemspider_api_key:
                    logging.warning("No ChemSpider API key found")
                    chemspider_api_key = input("Enter ChemSpider API key:")
                session=chemspipy.ChemSpider(chemspider_api_key)
                compound=session.get_compound(self.id)
                inchikey=compound.inchikey

            elif self.database == "pubchem":
                self.url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+self.id+"/property/InChIKey/TXT"
                response=self.lookup_url()
                inchikey=response.text.rstrip()

        except Exception as e:
            logging.warning(self.database+":"+self.id+": "+repr(e))

        return inchikey

class Compound:
    """ A compound within the MIBiG database """
    def __init__(self,mibig_accession,index,name,organism):
        self.mibig_accession: str = mibig_accession
        self.index: int = index
        self.name: str = name
        self.organism: str = organism
        self.database_ids: List[DatabaseID] = []
        self.inchikey: str = None
        self.mycotoxin: bool = False

    def is_mycotoxin(self,mycotoxin_inchikeys):
        """ Use associated database ids to find the inchikey and check if the compound is a mycotoxin """
        #Reorder database ids (to search in order)
        database_order = {"chembl":0,"chebi":1,"npatlas":2,"chemspider":3,"pubchem":4}
        self.database_ids.sort(key=lambda val: database_order[val.database])

        #Search databases for inchikey
        for database_id in self.database_ids:
            logging.info("Searching "+database_id.database+": "+database_id.id)
            self.inchikey=database_id.get_inchikey()
            #Stop searching once a valid inchikey is found
            if self.inchikey is not None:
                break

        #If the inchikey occurs in a list of mycotoxin inchikeys, add mycotoxin label
        if self.inchikey in mycotoxin_inchikeys:
            self.mycotoxin = True

    def summary(self) -> list:
        """ Create summary of compound properties """
        summary=[self.mibig_accession,self.organism,self.index,self.name,self.inchikey]
        return summary

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
            try:
                compound_name=bgc_data["cluster"]["compounds"][i]["compound"]
                #Store compound info in a new Compound object
                compound=Compound(mibig_accession,i,compound_name,organism_name)
                #Retrieve a list of all database ids
                id_list=bgc_data["cluster"]["compounds"][i]["database_id"]
                #Store ids as DatabaseID objects, and add to database_ids attribute of the Compound object
                for database_id in id_list:
                    database,id=database_id.split(":")
                    compound.database_ids.append(DatabaseID(database,id))
                #If the compound is a mycotoxin, store compound properties (MIBiG ID, index, name and inchikey)
                compound.is_mycotoxin(mycotoxin_inchikeys)
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
