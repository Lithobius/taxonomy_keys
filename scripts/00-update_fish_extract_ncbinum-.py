#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generating an NCBI search to get taxonomy from list of accession numbers.
This version is to update an existing search from NCBi with new terms only so that the entire search doesn't need rerun every time.

Input: the CRABS database, previous search
output: List of accession numbers that were not previously searched

Author: Kate Sheridan
2024 version 0.1.0

2024-march-11; work paused to try another method
"""

# load-in
from pyprojroot import here # for 'here' like R

# log-setup
# check logging package for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'update-ncbi-search_log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up classes and functions

# global variable for terms
already_found = []
to_search = []

# read in fasta
def get_terms_fasta(raw_fasta):
    # open old fasta to read and new fasta to write
    with open(raw_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                asvnum = line[1:-1]
                logging.info("found " + asvnum)
                if asvnum not in to_search:
                    to_search.append(asvnum)
                else:
                    pass
            else:
                pass
    #return acc_list

def new_terms_only(old_search):
    for line in old_search:
        logging.info("found: "  +line)

# write accessions to txt
def save_accns(outfile):
    with open(outfile, "w") as out:
        for acc in to_search:
            out.write(acc+"\n")
            logging.info('wrote: '+acc)


# script
if __name__ == "__main__":

    # input FASTA and txt
    # extract accession numbers
    crabs_fasta = here("./database12s/crabs/blastdb/20240310_12s_crabs.fasta")
    get_terms_fasta(crabs_fasta)

    # input existing list
    existing_terms = here(
        "./rawdata/12s_references/crabs/20230917_12s_combined-ncbi-accn.txt")
    new_terms_only(existing_terms)

    new_search = here("./rawdata/12s_references/crabs/20240311_12s_ncbi-accn_new.txt")
    save_accns(new_search)

# end
logging.info('Fin!')
