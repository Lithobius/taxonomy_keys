#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script extracts mitofish information from CRABS output

This assumes you've run CRABS with mitofish and assigned taxonomy to the original database.

By extracting all taxonomy from this output, it is possible to eliminate species that CRABS will extract from mitofish, since it is written in such a way that the whole mitofish database is always pulled without filter. Mitofish guarantees a 12S sequence of a length not too long for cutadapt, so these sequences should be preferred.

NOTE this does not yet filter for unique species, it is unique by accession number.

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
insert copyright, license, additional credits etc if necessary
"""

# module import
# built-in
import os

# common libraries
import csv
from pyprojroot import here  # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'mitofish-extract-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions


def extractMitofish(crabs_output, outfile_path):
    with open(crabs_output, 'r') as f:
        with open(outfile_path+'20230511_mitofish-species.csv', 'w') as out:
            mitofish_tsv = csv.reader(f, delimiter="\t")
            mitofish_csv = csv.writer(out)
            mitofish_csv.writerow(["accn", "uid", "family", "species"])
            for line in mitofish_tsv:
                accn = line[0]
                uid = line[1]
                family = line[6]
                species = line[8]
                mitofish_csv.writerow([accn, uid, family, species])
                #logging.info(line[0])
            # no headers
            # tab separated
            # extract column 0, 1, 6, 8
            # accession, ?uid?, family, species

# basic-body

if __name__ == "__main__":

    # read in mitofish output from crabs, extract species
    crabs_mitofish = ("processeddata/species_lists/ncbi/20230511_mitofishdb_tax.tsv")
    mitofish_out_path = ("./processeddata/species_lists/ncbi/")
    extractMitofish(crabs_mitofish, mitofish_out_path)



# end
logging.info('Fin!')
