#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This updates our output from 06b to match taxonomy
Currently input from custom databases gets 'CRABS:' appended.
This is potentially a bug in CRABS to report in github issues.

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
from pyprojroot import here  # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'format-taxonomy-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# set up functions
def removeCrabs(crabs_output, ready_out_path):
    with open(crabs_output, 'r') as f:
        with open(ready_out_path + '20230511_pgaout.fasta', 'w') as out:
            for line in f:
                # save sequence
                seq = next(f)
                if line.startswith('>'):
                    # extract ID
                    accn = line[1:-1]
                    # modify ID
                    accn_mod = line.replace('CRABS:', '')
                    #logging.info(accn_mod)
                    # writes accnmod newline then sequence then newline
                    out.write(accn_mod + seq)
                else:
                    pass

# basic-body

if __name__ == "__main__":

    # read in fasta from crabs: pga results
    # adjust fasta headers
    crabs_mitofish = ("processeddata/crabs/pgaout.fasta")
    crabsready_out_path = ("./processeddata/crabs/")
    removeCrabs(crabs_mitofish, crabsready_out_path)

# end
logging.info('Fin!')
