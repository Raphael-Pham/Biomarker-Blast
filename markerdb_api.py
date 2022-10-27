import requests
import json
import constants

def markerdb_get(req_type,
                 QUERY=None,
                 NAME=None,
                 PAGE=None,
                 MARKERDB_ID=None,
                 CATEGORY=None,
                 BIOMARKER_TYPE=None):
    URL = constants.MARKERDB_API
    PARAMS = {}
    API_KEY = "cc69fb7c1a4809bbae7c0701d5188c17"

    if req_type == "condition":
        URL += "conditionapi/conditionrequest"
        PARAMS = { 'api_key': API_KEY, 'query': QUERY, 'page': PAGE }
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
    return r.json()

def main():
    print(markerdb_get("chemical", NAME="1-Methylhistidine", MARKERDB_ID="MDB00000001", PAGE=1))

if __name__ == '__main__':
    main()