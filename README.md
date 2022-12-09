# Biomarker Blast
Biomarker Blast is a useful tool for clinicians to predict susceptibility for diseases with an alignment comparison to a reference gene, finding genetic biomarkers such as single nucleotide polymorphisms and frameshift mutations.

### **Prerequisites**
Python 3.x+ must be installed.

### **Setup**
Install the necessary packages by running (if not already installed locally):
```
pip install bs4 lxml cchardet biopython
```
Next, run:

*Linux*:
```
sudo apt install ncbi-blast+
```
*Mac*:
```
brew install blast
```

Then, add a sequence of interest's FASTA file to the `known_genomes` directory.

### **Execution**
To run the tool, simply type into the command line:
```
python3 biomarker_blast.py -d [DISEASE_NAME] -f [FILE_NAME]
```
where DISEASE_NAME is the disease/condition of interest, and
FILE_NAME is the name of the file containing the sequence of interest (FASTA format).

For help with using this tool:
```
python3 biomarker_blast.py -h
```


### **Sample Run**
Here is an example of how to use this tool:

```
python3 biomarker_blast.py -d Usher Syndrome Type I -f known_seq.fasta
```
Sample output:

```
Successfully retrieved Usher Syndrome Type I data from MarkerDB...
Fetching all associated biomarker GenBank IDs. Please be patient...

Unique biomarker IDs: {'64072', '65217', '83715', '124590', '4647', '10083'}

Saving fasta files...
/genbank_fasta/NG_008835.1_4974_to_424001.fasta saved!
...

Printing alignment hit results:
1. File name:  NG_008835.1_4974_to_424001.fasta
Score:  59.8958 bits (65.0)
Expect:  5.91473e-08
Identities:  34/35 (97.1%)
Positives:  34/35 (97.1%)
Gaps:  0/35 (0.0%)
Query  62121  CCTCCCAAAGTGCTGGGATTACAGGTGTGAGCCAC  62155
Sbjct  37359  CCTCCCAAAGTGCTGGGATTACAGGTGTAAGCCAC  37325
...
```