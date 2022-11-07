import os
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

BLAST_TYPE = "blastn"
FILE_1 = "genbank_fasta/NG_007882.2_5009_to_12181.fasta"
FILE_2 = "known_genomes/known_seq.fasta"

# input: string, string
# output: iterator
def convert_file_to_record_iterator(in_file, file_type):
    return SeqIO.read(in_file, format=file_type)

# input: string, string, string 
# output: None (BLAST output is printed to screen)
def blast_sequences(input_seq, query_seq, db):
    result = NCBIWWW.qblast(program=BLAST_TYPE, database=db, sequence=input_seq, query_file=query_seq)
    print(result.read())

# def main():
#     print("converting files...")
#     condition = convert_file_to_record_iterator(FILE_1, "fasta")
#     known_genome = convert_file_to_record_iterator(FILE_2, "fasta")

#     print("blast-ing sequences...")
#     blast_sequences(condition, known_genome, "nt")

def main():
    blastn_path = os.path.dirname(os.path.realpath(__file__)) + "/packages/blast-BLAST_VERSION+/bin/blastn.exe"
    output = NcbiblastnCommandline(cmd=blastn_path, query=FILE_1, subject=FILE_2, outfmt=5)()[0]
    blast_output_record = NCBIXML.read(StringIO(output))
    for alignment in blast_output_record.alignments:
        for hsp in alignment.hsps:
            print()
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

if __name__ == "__main__":
    main()