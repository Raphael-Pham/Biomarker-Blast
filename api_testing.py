import requests
import json

# MarkerDB
def markerdb_get(req_type,
                 NAME=None,
                 PAGE=None,
                 MARKERDB_ID=None,
                 CATEGORY=None,
                 BIOMARKER_TYPE=None):
    URL = "http://markerdb.ca/api/v1/"
    PARAMS = {}
    API_KEY = "cc69fb7c1a4809bbae7c0701d5188c17"

    if req_type == "condition":
        URL += "geneapi/generequest"
        PARAMS = { 'api_key': API_KEY, 'name': NAME, 'page': PAGE }
    elif req_type == "chemical":
        URL += "chemicalapi/chemicalrequest"
        PARAMS = { 'api_key': API_KEY, 'name': NAME, 'markerdb_id': MARKERDB_ID, 'page': PAGE }
    elif req_type == "gene":
        URL += "geneapi/generequest"
        PARAMS = { 'api_key': API_KEY, 'name': NAME, 'markerdb_id': MARKERDB_ID, 'page': PAGE }
    elif req_type == "protein":
        URL += "proteinapi/proteinrequest"
        PARAMS = { 'api_key': API_KEY, 'name': NAME, 'markerdb_id': MARKERDB_ID, 'page': PAGE }
    elif req_type == "karyotype":
        URL += "karyotypeapi/karyotyperequest"
        PARAMS = { 'api_key': API_KEY, 'markerdb_id': MARKERDB_ID}
    elif req_type == "all":
        URL += "generalapi/generalrequest"
        PARAMS = { 'api_key': API_KEY, 'category': CATEGORY, 'biomarker_type': BIOMARKER_TYPE, 'page': PAGE }
    else:
        print("Invalid request type!")

    r = requests.get(url=URL, params=PARAMS)
    return json.dumps(r.json(), sort_keys=True, indent=2)


def fetch_all_info(condition_name):
    URL = "http://markerdb.ca/api/v1/generalapi/generalrequest"
    PARAMS = { 'api_key': "cc69fb7c1a4809bbae7c0701d5188c17",
               'category': 'Predictive',
               'biomarker_type': 'SequenceVariant',
               'page': 1 }
    r = requests.get(url = URL, params = PARAMS)
    return json.dumps(r.json(), sort_keys=True, indent=3)

def get_gene_ids(biomarkers):                           # ! ISSUE: how to check all available pages??
    gene_ids = []
    for biomarker_id, biomarker_name, biomarker_type in biomarkers:
        if biomarker_type == "Chemical":                # ! no gene_id property
            #gene_id = markerdb_get("chemical", NAME=biomarker_name, MARKERDB_ID=biomarker_id, PAGE=1)
            #gene_ids += gene_id
            pass
        elif biomarker_type == "Protein":               # ! no gene_id property
            #gene_id = markerdb_get("protein", NAME=biomarker_name, MARKERDB_ID=biomarker_id, PAGE=1)
            #gene_ids += gene_id
            pass
        elif biomarker_type == "Genetic":
            gene_ids += markerdb_get("gene", NAME=biomarker_name, MARKERDB_ID=biomarker_id, PAGE=1)["gene"]["id"]
        elif biomarker_type == "Karyotype":
            gene_ids += markerdb_get("karyotype", MARKERDB_ID=biomarker_id)["gene"]["id"]
        else:
            print("Invalid biomarker type!")
            break
    return gene_ids

def main():
    #all_data = fetch_all_info("Alzheimer's Disease")
    #print(all_data)

    print(markerdb_get("chemical", NAME="1-Methylhistidine", MARKERDB_ID="MDB00000001", PAGE=1))

'''
    1. Query all biomarkers to match a single condition name
    2. Iterate through list of biomarkers, querying corresponding (biomarker type)
        API call to fetch gene_id
    3. Return list of gene_id's to be looked up on GenBank database
'''

'''
    1. Cystic fibrosis (chemical)
    2. Sickle cell (chemical + protein) 
    3. Working: down syndrome
'''

if __name__ == '__main__':
    main()