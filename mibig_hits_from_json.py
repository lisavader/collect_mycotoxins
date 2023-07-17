import json
import argparse


def cc_mibig_topresult(antismash_results: dict, record: int, region: int, cutoff: float=0.9):
    #For the specified record and region, extract MIBiG clustercompare RegionToRegion results from json
    cc_mibig_RTR=antismash_results["records"][record]["modules"]["antismash.modules.cluster_compare"]["db_results"]["MIBiG"]["by_region"][str(region)]["RegionToRegion_RiQ"]
    #The score and reference BGC for the top hit
    topscore=list(cc_mibig_RTR["scores_by_region"].values())[0]
    reference=list(cc_mibig_RTR["scores_by_region"].keys())[0]
    #More information related to the reference BGC
    reference_info=cc_mibig_RTR["reference_regions"][reference]
    compound=reference_info["description"]
    organism=reference_info["organism"]

    if topscore >= cutoff:
        topresult=(topscore,reference,compound,organism)
        return(topresult)


def main(json_path):
    #Import json file with antismash results
    with open(json_path, "r") as json_file:
        antismash_results = json.load(json_file)
    #Loop over records (contigs)
    last_record=len(antismash_results["records"])
    for record in range(0,last_record):
        #Skip if the record contains no BGCs:
        if not antismash_results["records"][record]["areas"]:
            pass
        else:
            #Loop over regions (BGCs)
            last_region=len(antismash_results["records"][record]["areas"])
            for region in range(1,last_region+1):
                #calculate top mibig result
                result=cc_mibig_topresult(antismash_results,record,region)
                if result:
                    print(result)
        
  
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser("mycotoxins_from_json")
    parser.add_argument("json_path", type=str, help="Path to a json file containing antiSMASH results. AntiSMASH options required: cc-mibig.")
    args = parser.parse_args()
    main(args.json_path)
