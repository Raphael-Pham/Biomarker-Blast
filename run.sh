#!/bin/bash

export PYTHONPATH="$PWD/packages"
clear
python3 markerdb_webscraping.py
python3 blast_api.py