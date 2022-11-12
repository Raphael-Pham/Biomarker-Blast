import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

import constants

FILE_2 = "known_genomes/known_seq.fasta"

def main():
    blastn_path = os.path.dirname(os.path.realpath(__file__)) + "/packages/blast-BLAST_VERSION+/bin/blastn.exe"
    saved_fastas = open('saved_fasta_names.txt')
    fasta_names = [name.strip() for name in saved_fastas.readlines()]
    hsp_results = []

    for fasta_file in fasta_names:
        output = NcbiblastnCommandline(cmd=blastn_path, query=fasta_file, subject=FILE_2, task=constants.BLAST_TYPE, evalue=0.000002, outfmt=5)()[0]
        blast_output_record = NCBIXML.read(StringIO(output))

        for alignment in blast_output_record.alignments:
            sorted_hsps = sorted(alignment.hsps, key=lambda x: x.identities / x.align_length, reverse=True)
            hsp = sorted_hsps[0]
            hsp_results.append((fasta_file.split('/')[1], hsp))

    print(hsp_results)

if __name__ == "__main__":
    main()