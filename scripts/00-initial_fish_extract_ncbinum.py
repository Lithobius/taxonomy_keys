#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generating an NCBI search to get taxonomy from list of accession numbers

Author: Kate Sheridan
2022 version 0.1.0
"""

# load-in
from pyprojroot import here # for 'here' like R

# log-setup
# check logging package for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'fastanames2csv-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up classes and functions

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


# read in txt
def get_terms_txt(txtfile):
    with open(txtfile) as f:
        for line in f:
            asvnum = f.readline()
            logging.info("found " + asvnum)
            if asvnum not in to_search:
                to_search.append(asvnum)
            else:
                pass

# write accessions to txt
def save_accns(outfile):
    with open(outfile, "w") as out:
        for acc in to_search:
            out.write(acc+"\n")
            logging.info('wrote: '+acc)


# script
if __name__ == "__main__":

    # list of search terms
    to_search = []
    # input FASTA and txt
    mitofish_fasta = here("./rawdata/12s_references/crabs/mitofishdb.fasta")
    outgroups_fasta = here("./rawdata/12s_references/crabs/20230905_12s_outgroups-ncbi-filter2.fasta")
    ncbifish_fasta = here("./rawdata/12s_references/crabs/20230910_12s_fish-ncbi-filter.fasta")
    ncbifish_toolong = here("./rawdata/12s_references/crabs/20230910_12s_ncbi-fish_toolong.txt")
    # extract accession numbers
    get_terms_fasta(mitofish_fasta)
    get_terms_fasta(outgroups_fasta)
    get_terms_fasta(ncbifish_fasta)
    get_terms_txt(ncbifish_toolong)
    # outfile
    search_terms = here(
        "./rawdata/12s_references/crabs/20230917_12s_combined-ncbi-accn.txt")
    save_accns(search_terms)

# end
logging.info('Fin!')
