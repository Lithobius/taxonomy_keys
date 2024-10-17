#!/bin/bash

#==============================================================================
#title        :  01d_crabs_mitofish-and-taxonomy.sh
#description  :  This script uses crabs to download mitofish and the NCBI taxonomy file. Note this taxonomy file is required for future steps of taxonomy assignment, and sometimes updates from NCBI.
#author       :  Kate Sheridan
#date         :  20230912
#version      :  0.1.0
#notes        :  requires CRABS.
#==============================================================================

# starts in project root
# move to folder for data download/processing
cd database12s/crabs

# 1: db_download
# crabs can download from several databases

# mitofish db is just all of mitofish
crabs db_download --source mitofish --output mitofishdb.fasta --keep_original yes

# taxonomy
# This file needs to be downloaded if we use
# NCBI etc data, even a custom database
crabs db_download --source taxonomy
