#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is for extracting a csv of names out of an NCBI download.
Specifically: send to > complete record > file > tinyseq

The output of that extraction is a xml file with numbers, name of record, accession numbers, base pairs etc.

Here we want to get accession numbers, taxa ID, and taxa names only. They can be messy since they'll be run through a cleaning and higher taxonomy script in R later.
But we'll be able to run the result through 'classify' from taxize to give us the full taxonomy from NCBI for each taxonomy id, then fix it with worrms, then match fixed taxonomy to accession number.

Right now this solution has an extra column added as some entries shift the taxonomy id and organism name over by one for how we're calling the fields. It will be fixed during loadin in R. In the future I'll update this to call the fields directly.


Author: Kate Sheridan
2023 August version 0.1.0
"""

# load-in
import os
import csv
import xml.etree.ElementTree as ET
from pyprojroot import here # for 'here' like R

# log-setup
# check logging package for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'ncbi-extract-tinyseq-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# read in xml, sort, write out
def taxonomymatch(ncbi_xml, to_export):
    # make xml etree
    tree = ET.parse(ncbi_xml)
    root = tree.getroot()
    logging.info(root[0][1].text)
    with open(to_export, 'w') as out:
        ncbiwriter = csv.writer(out, delimiter = ",")
        # write header
        ncbiwriter.writerow(["accession", "tax_id", "organism", "defline"])
        for i in root:
            # write: accession number, taxonomy id, organism, bonus field for later adjustment
            ncbiwriter.writerow([i[1].text, i[2].text, i[3].text, i[4].text])
            logging.info(i[2].text)


# script

if __name__ == "__main__":

    # input FASTA
    ncbi_file = here("./species_lists/ncbi/2023_ncbi/20230819_ncbi_fish_12s_tinyseq.xml")
    export = here("./species_lists/ncbi/2023_ncbi/20230828_ncbi_fish_12s_tinyseq-extract.csv")
    # run fasta function
    taxonomymatch(ncbi_file, export)





# end
logging.info('Fin!')
