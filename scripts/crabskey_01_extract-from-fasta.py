#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Make a key based only on the CRABS database.
Extract fasta headers and format into CSV for taxize/worrms.

Input: FASTA from CRABS database
Output: CSV with headers accession, taxonomy

Author: Kate Sheridan
2024 version 0.1.0
"""

# load-in
import os
import csv
from pyprojroot import here # for 'here' like R

# log-setup
# check logging package for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'crabsdb-extractfasta_log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')

# set up classes and functions

# read in fasta, write out csv
def name_extract(raw_fasta, csv2write):
    # open old fasta to read and new fasta to write
    with open(raw_fasta, 'r') as f, \
        open(csv2write, 'w') as out:
        accwriter = csv.writer(out, delimiter = ",")
        accwriter.writerow(["accession", "taxonomy"])
        for line in f:
            if line.startswith('>'):
                seqname = line[1:-1]
                logging.info("found " + seqname)
                seq_acc = seqname.split("\t")
                logging.info(seq_acc)
                accwriter.writerow(seq_acc)
            else:
                pass


# script
if __name__ == "__main__":

    # input FASTA
    fasta = here("./database12s/crabs/blastdb_make/20241015_crabs_mifish.fasta")
    # csv to write out
    new_csv = here("./database12s/crabs/20241015_12s_crabs_headers.csv")
    # run fasta function
    name_extract(fasta, new_csv)


# end
logging.info('Fin!')
