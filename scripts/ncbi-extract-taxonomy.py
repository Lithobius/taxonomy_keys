#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is for extracting a csv of names out of an NCBI download.
Specifically: send to > complete record > file > summary

The output of that extraction is a text file with numbers, name of record, accession numbers, base pairs etc.

Here we want to get accession numbers and taxa names only. They can be messy since they'll be run through a cleaning and higher taxonomy script in R later.

Any license info also goes here
Author: "Kate Sheridan"
Version: "0.1.0"
Status: dev
"""

# module import
# built-in
import os

# common libraries
import csv
import re
from pyprojroot import here # like here in R

# less-known libraries
# own libraries/functions

# log-setup
# check logging module for more settings
# I'm mostly using the debug options for development
import logging

LOG_FILENAME = 'ncbi-extract-log.txt'

logging.basicConfig(filename=LOG_FILENAME,
                    level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s - %(message)s')
logging.info('Debut!')


# Storing the name and accession fields
class ncbiRecord(object):

    def __init__(self, name, acc):
        self.name = name
        self.acc = acc

# set up functions

# first function reads in the NCBI records and makes an object
def readNCBI(ncbi_file):
    with open(ncbi_file, 'r') as f:
        for line in f:
            if line.strip() == "":
                logging.info('blank')
            #if line.startswith('@'):
                pass
            else:
                fields = line.strip().split(',')
                name = fields[0]
                # advance two lines to get accession
                bp = next(f, "none")
                acc = next(f, "none")
                logging.debug(name + " " + acc)
                yield ncbiRecord(name, acc)


def writeNCBI(records):
    with open(here("./species_lists/ncbi/2023_ncbi/20230828_ncbi_fish_12s_extract.csv"), 'w') as out:
        ncbiwriter = csv.writer(out)
        for r in records:
            # get a regex for the most obvious offenders
            # we'll strip out the random capital letters, etc. later
            namewrite = re.sub('[0-9]|sp\\.|\\.|S ribosomal RNA gene|isolate|voucher|mitochondrion|mitochondrial|:|-|gene|for|rRNA|_|ribosomal RNA|and|tRNA-Val genes|small subunit|FACEN|DNA|tRNA|Kuroshio Biological Research Foundation|partial|MParis|cytochrome b| \\(CytB\\)|\\(cytb\\)|CYTB|cytochrome c|\\(eral\\)|bwicor|haplotype|sspo|Phe |genome| genomic| complete | chaperone| unplaced| scaffold| \\(|\\)| chromosome| part |Hap |mitochondiral|PREDICTED| seed storage| protein|methyltransferase|methylcytidine| adult heart library| prime similar to|arachidonate| lipoxygenase|strain|Dloop|LodgeLab|Contig|linkage group|UNVERIFIED| sequence|\\<|//>|Chiba museum| Zap Express library| similar to| transcribed RNA|haplogroup| personal|mitochrodrial| clone|control region|collectiondate|popvariant| contains |aquarium| Fish ', "", r.name)
            acc2 = r.acc.split()
            accwrite = acc2[0]
            ncbiwriter.writerow([namewrite,accwrite])
            logging.debug(namewrite + " " + accwrite)



if __name__ == "__main__":

    print("hello world")
    ncbi_in = here("./species_lists/ncbi/2023_ncbi/20230828_ncbi_fish_12s.txt")
    #ncbi_in = here("./rawdata/12s_references/ncbi_test.txt")
    ncbi_data = readNCBI(ncbi_in)
    writeNCBI(ncbi_data)


# end
logging.info('Fin!')
