import os
import sys
import requests
from bs4 import BeautifulSoup, SoupStrainer
import lxml
import cchardet
import re
import time
from Bio import Entrez, SeqIO
from Bio.SeqIO import FastaIO

from markerdb_api import markerdb_get
import constants

def build_soup(url, requests_session, type):
    page = requests_session.get(url)

    if type == "first_biomarker":
        return BeautifulSoup(page.text, "lxml")
    elif type == "biomarker":
        return BeautifulSoup(page.text, "lxml", parse_only=SoupStrainer("div", {"id": "Sequence_Variants"}))
    elif type == "sequence_variant" or type == "genbank":
        return BeautifulSoup(page.text, parse_only=SoupStrainer("a"), features="xml")
    elif type == "fasta":
        return BeautifulSoup(page.text, features="xml")

def get_biomarker_ids(condition_id, requests_session):
    gene_ids = set()
    page_num = 1

    soup = build_soup(constants.MARKERDB_CONDITIONS + condition_id + "?page=1", requests_session, "first_biomarker")
    if "We're sorry, but something went wrong." in str(soup.find_all("div", class_="dialog")):
        print("Condition not available in MarkerDB.")
        exit(1)

    print()
    print("Biomarker fetching progress:")
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
        soup = build_soup(constants.MARKERDB_CONDITIONS + condition_id + "?page=" + str(page_num), requests_session,
            "biomarker")
        seq_variants = soup.find_all("div", {"id": "Sequence_Variants"}, limit=1)
    print()

    return gene_ids

def save_gene_sequences(biomarker_ids, requests_session):
    print()
    print("Saving fasta files...")
    saved_file_names = []
    saved_fastas = open('saved_fasta_names.txt', 'w')

    for biomarker_id in biomarker_ids:
        soup = build_soup(constants.GENBANK_GENE + biomarker_id, requests_session, "genbank")
        ncbi_string = r'(.*)ncbi_uid=' + re.escape(str(biomarker_id)) + r'(.*)report=fasta'

        file_extension = soup.find_all('a', ref=re.compile(ncbi_string), limit=1)[0]['href']
        accession_num = file_extension[9:].split('?')[0]
        from_to = file_extension.split('&')
        from_idx = from_to[0].split('=')[-1]
        to_idx = from_to[1].split('=')[-1]

        handle = Entrez.efetch(db='nuccore', id=accession_num, seq_start=from_idx, seq_stop=to_idx, rettype='fasta', retmode='text')
        fasta_filename = "/genbank_fasta/" + accession_num + '_' + from_idx + '_' + "to" + '_' + to_idx + ".fasta"
        SeqIO.write(SeqIO.read(handle, 'fasta'), os.path.dirname(os.path.realpath(__file__)) + fasta_filename, 'fasta')
        time.sleep(1)                                   # Prevent NCBI from closing connection due to frequency of efetch requests
        print(fasta_filename + " saved!")
        saved_file_names.append(fasta_filename[1:] + "\n")
    
    saved_fastas.writelines(saved_file_names)
    saved_fastas.close()

def scrape_marker_db(user_input):
    requests_session = requests.Session()
    Entrez.email = "raphaelpham14@gmail.com"

    condition_html = markerdb_get("condition", QUERY=user_input, PAGE=1)['conditions'][0]
    print("\nRetrieved " + condition_html['name'] + " data...")
    print("Fetching all associated biomarker GenBank IDs...")
    biomarker_ids = get_biomarker_ids(str(condition_html['id']), requests_session)
    print("Biomarker IDs: " + str(biomarker_ids))
    save_gene_sequences(biomarker_ids, requests_session)

if __name__ == "__main__":
    scrape_marker_db(' '.join(sys.argv[1:]))