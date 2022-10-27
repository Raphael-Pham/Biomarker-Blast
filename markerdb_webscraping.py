import requests
from bs4 import BeautifulSoup, SoupStrainer
import lxml
import cchardet
import re
import time

from api_testing import markerdb_get

def build_soup(url, requests_session, type):
    page = requests_session.get(url)

    if type == "first_biomarker":
        return BeautifulSoup(page.text, "lxml")
    elif type == "biomarker":
        return BeautifulSoup(page.text, "lxml", parse_only=SoupStrainer("div", {"id": "Sequence_Variants"}))
    elif type == "sequence_variant":
            return BeautifulSoup(page.text, "lxml", parse_only=SoupStrainer("a"))

def get_biomarker_ids(condition_id, requests_session):
    gene_ids = set()
    page_num = 1

    soup = build_soup("https://markerdb.ca/conditions/" + condition_id + "?page=1", requests_session, "first_biomarker")
    if "We're sorry, but something went wrong." in str(soup.find_all("div", class_="dialog")):
        print("Condition not available in MarkerDB.")
        exit(1)

    seq_variants = soup.find_all("div", {"id": "Sequence_Variants"}, limit=1)
    while len(seq_variants) > 0:
        seq_variants = seq_variants[0].find_all('a', href=re.compile(r'(.*)sequence_variants(.*)'))

        for variant in seq_variants:
            biomarker_id = str(variant).split('"')[1].split("/")[-1]
            soup = build_soup("https://markerdb.ca/sequence_variants/" + biomarker_id, requests_session, "sequence_variant")
            gene_id = soup.contents[-1].text.split(' ')[-1]
            gene_ids.add(gene_id)
        
        if page_num % 5 == 0:
            print(page_num)
            
        page_num += 1
        soup = build_soup("https://markerdb.ca/conditions/" + condition_id + "?page=" + str(page_num), requests_session,
            "biomarker")
        seq_variants = soup.find_all("div", {"id": "Sequence_Variants"}, limit=1)

    return gene_ids
    
def scrape_marker_db():
    requests_session = requests.Session()

    condition_id = str(markerdb_get("condition", QUERY="Homocystinuria", PAGE=1)['conditions'][0]['id'])
    biomarker_ids = get_biomarker_ids(condition_id, requests_session)
    print(biomarker_ids)


if __name__ == "__main__":
    scrape_marker_db()