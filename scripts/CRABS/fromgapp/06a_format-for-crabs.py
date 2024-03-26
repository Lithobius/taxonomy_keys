#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This updates our output from 05 to read into crabs

note: this is only necessary if the output from ncbi has decimals or any other info

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
insert copyright, license, additional credits etc if necessary
"""

# module import
# built-in
import os
import re

# common libraries
from pyprojroot import here  # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'crabsready-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')

# set up functions
def extractMitofish(crabs_output, ready_out):
    with open(crabs_output, 'r') as f:
        with open(ready_out, 'w') as out:
            for line in f:
                if line.startswith('>'):
                    # extract ID
                    accn = line[1:-1]
                    # modify ID
                    accn_mod = re.sub('\.[1-9]', '', line)
                    #logging.info(accn_mod)
                    # writes accnmod newline then sequence then newline
                    out.write(accn_mod)
                else:
                    out.write(line)

# basic-body

if __name__ == "__main__":

    # read in fasta from crabs: pga results
    # adjust fasta headers
    fish_ncbi = ("database12s/crabs/20240212_12s_fish-ncbi-filter.fasta")
    crabsready_fish = ("./database12s/crabs/20240223_fish-short_crabs-ready.fasta")
    extractMitofish(fish_ncbi, crabsready_fish)

    outgroups_ncbi = ("database12s/crabs/20240212_12s_outgroups-ncbi-filter2.fasta")
    crabsready_outgroups = ("./database12s/crabs/20240223_outgroups-short_crabs-ready.fasta")
    extractMitofish(outgroups_ncbi, crabsready_outgroups)




# end
logging.info('Fin!')
