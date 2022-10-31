from Bio import SeqIO
from Bio.Blast import NCBIWWW

BLAST_TYPE = "blastp"
FILE_1 = "condition.gbk"
FILE_2 = "known_genome.faa"

# input: string, string
# output: iterator
def convert_file_to_record_iterator(in_file, file_type):
    return SeqIO.read(in_file, format=file_type)

# input: string, string, string 
# output: None (BLAST output is printed to screen)
def blast_sequences(input_seq, query_seq, db):
    result = NCBIWWW.qblast(BLAST_TYPE, db, input_seq, query_file=query_seq)
    print(result.read())

def main():
    print("converting files...")
    condition = convert_file_to_record_iterator(FILE_1, "genbank")
    known_genome = convert_file_to_record_iterator(FILE_2, "fasta")

    print("blast-ing sequences...")
    blast_sequences(condition, known_genome, "nt")

if __name__ == "__main__":
    main()