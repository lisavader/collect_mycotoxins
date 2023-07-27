import argparse
import glob
import os
import json
import csv

def list_missing_database_ids(mibig_path):
    missing_database_ids = []
    mibig_entries=set(glob.glob(mibig_path+"/*.json"))
    for entry in mibig_entries:
        with open(entry, 'r') as json_file:
            bgc_data = json.load(json_file)
        #Get name of mibig accession from file path
        filename=os.path.basename(entry)
        mibig_accession=os.path.splitext(filename)[0]
        #Find all associated compounds
        compound_count=len(bgc_data["cluster"]["compounds"])
        for i in range(0,compound_count):
            compound_name=bgc_data["cluster"]["compounds"][i]["compound"]
            try:
                id_list=bgc_data["cluster"]["compounds"][i]["database_id"]
            except KeyError:
                missing_database_ids.append([mibig_accession,i,compound_name])
    return missing_database_ids

def main(mibig_path):
    missing_database_ids=list_missing_database_ids(mibig_path)
    with open("results/missing_database_ids.csv","w") as file:
        writer = csv.writer(file)
        header = ["accession","compound_index","compound_name"]
        writer.writerow(header)
        writer.writerows(missing_database_ids)

if __name__ == "__main__":
    #Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("mibig_path", type=str, help="Path to a MIBiG database of annotations in .json format")
    args = parser.parse_args()
    #Run the main script
    main(args.mibig_path)
