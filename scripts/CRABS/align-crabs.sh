#!/bin/bash

#==============================================================================
#title        :  align-crabs.sh
#description  :  sequence alignment for the database produced by CRABS, to investigate the missing sequences
#author       :  Kate Sheridan
#date         :  20230912
#version      :  0.1.0
#notes        :  requires CRABS.
#==============================================================================

# starts in project root
# move to folder for data download/processing
cd rawdata/12s_references/crabs

# sequence alignment
clustalo --in mitofishdb.fasta --out 20230917_mitofish-align.fasta -v
