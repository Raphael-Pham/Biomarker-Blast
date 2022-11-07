#!/bin/bash

export PYTHONPATH="$PWD/packages"
#export PATH=$PATH:/$PWD/packages/blast-BLAST_VERSION+/bin

python markerdb_webscraping.py Usher Syndrome Type I
python blast_api.py