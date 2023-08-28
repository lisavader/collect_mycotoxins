import requests
import time
import chemspipy
import re
import logging
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
        self.mycotoxin: bool = None

    def inchikey_from_database(self):
        """ Use associated database ids to find the inchikey"""
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
        return self.inchikey

    def mycotoxin_from_inchikey(self,mycotoxin_inchikeys):
        if self.inchikey is None:
            self.inchikey_from_database()
        if self.inchikey in mycotoxin_inchikeys:
            self.mycotoxin = True
        else:
            self.mycotoxin = False
        return self.mycotoxin

    def mycotoxin_from_name(self,mycotoxin_names,exact_match = False):
        for mycotoxin_name in mycotoxin_names:
            if exact_match == True:
                if mycotoxin_name.casefold() == self.name.casefold():
                    self.mycotoxin = True
            elif exact_match == False:   #partial match (e.g. matches aflatoxin in compound name 'Aflatoxin A')
                if re.search(mycotoxin_name,self.name,re.IGNORECASE):
                    self.mycotoxin = True
        if self.mycotoxin != True:
            self.mycotoxin = False
        return self.mycotoxin