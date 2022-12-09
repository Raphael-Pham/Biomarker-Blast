from markerdb_webscraping import scrape_marker_db
from blast_api import get_blast_results
import sys
from os.path import exists

def main():
    args = sys.argv[1:]
    condition_name, file_name = parse_user_input(args)

    scrape_marker_db(condition_name)
    hsps = get_blast_results(file_name)
    display_hsp_results(hsps)

def parse_user_input(args):
    start_disease_idx, end_disease_idx, start_file_idx = -1, -1, -1
    
    # check if there is a help flag
    if args[0] == "-h" and len(args) == 1:
        help()
    elif "-d" not in args:
        raise TypeError("Missing -d flag.")
    elif "-f" not in args:
        raise TypeError("Missing -f flag.")
    
    # split args by flags
    for i in range(len(args)):
        if args[i] == "-d":
            start_disease_idx = i + 1
        elif args[i] == "-f":
            end_disease_idx = i
            start_file_idx = i + 1

    if args[start_disease_idx:end_disease_idx+1] == ["-f"]:
        raise TypeError("Missing disease name.")
    elif start_file_idx == len(args):
        raise TypeError("Missing file name.")
    
    condition_name = ' '.join(args[start_disease_idx:end_disease_idx])
    file_name = args[start_file_idx:][0]

    if not exists('known_genomes/' + file_name):
        raise TypeError("File " + file_name + " does not exist within the known_genomes directory.")

    return condition_name, file_name

def help():
    print("Usage: python3 biomarker_blast.py -d [DISEASE_NAME] -f [FILE_NAME]")
    sys.exit(0)

def display_hsp_results(hsps):
    print("Printing alignment hit results:")
    i = 1
    for file_name, hsp in hsps:
        print(str(i) + '. ' + "File name: ", file_name)
        print("Score: ", str(hsp.bits) + " bits " + '(' + str(hsp.score) + ')') #
        print("Expect: ", hsp.expect) 
        print("Identities: ", str(hsp.identities) + '/' + str(hsp.align_length) + ' (' + str(round(hsp.identities/hsp.align_length*100, 1)) + '%' + ')') 
        print("Positives: ", str(hsp.positives) + '/' + str(hsp.align_length) + ' (' + str(round(hsp.positives/hsp.align_length*100, 1)) + '%' + ')') 
        print("Gaps: ", str(hsp.gaps) + '/' + str(hsp.align_length) + ' (' + str(round(hsp.gaps/hsp.align_length*100, 1)) + '%' + ')') 

        print("Query  " + str(hsp.query_start) + "  " + hsp.query + "  " + str(hsp.query_end))
        print("Sbjct  " + str(hsp.sbjct_start) + "  " + hsp.sbjct + "  " + str(hsp.sbjct_end))
        print()

        i += 1


if __name__ == '__main__':
    main()