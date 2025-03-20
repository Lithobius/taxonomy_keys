#!/bin/bash
#==============================================================================
#title        :  04_crabs_visualizations.sh
#description  :  This script uses the visualization parameters in CRABS to view parameters of the custom database. Input - cleaned database from step 02
#author       :  Kate Sheridan
#date         :  20230912
#version      :  0.1.0
#notes        :  requires CRABS.
#==============================================================================

# starts in project root
# move to folder for data download/processing
cd rawdata/12s_references/crabs

# 8 visualization
# basic initial visuals
# diversity by taxonomic level
# barplot with bar for # species and # sequences
crabs visualization --method diversity --input 20230905_dbcleaned.tsv --level order

# length of sequences, frequency plot
crabs visualization --method amplicon_length --input 20230905_dbcleaned.tsv --level order

# db_completeness
# you have to provide the species names in a txt
# gives you summary stats for sp of interest; assesses reference database but also all of NCBI
crabs visualization --method db_completeness --input 20230905_dbcleaned.tsv --output 20230905_mifishcomplete.tsv --species species.txt --taxid nodes.dmp --name names.dmp

# primer_efficiency
# makes a bar graph with proportion base pair occurrences
# generates a fasta with sequences that contributed to graph
crabs visualization --method primer_efficiency --input 20230905_dbcleaned.tsv --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --fwd_name mifish_miya_f --rev_name mifish_miya_r --raw_file 20230905_12s_crabs-merge.fasta --tax_group Gobiidae --output gobiidae_bind.fasta

# phylo
# not working for me right now but should build a tree
# I think I would rather do this myself tbh
## yeah muscle is still giving me an error
#crabs visualization --method phylo --input 20230905_dbcleaned.tsv --level family --species species.txt --taxid nodes.dmp --name names.dmp
