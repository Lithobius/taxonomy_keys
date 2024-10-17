#!/bin/bash

# starts in project root
# move to folder for data download/processing
cd processeddata/crabs

# 1: db_download
# crabs can download from several databases

# mitofish db is just all of mitofish
crabs db_download --source mitofish --output mitofishdb.fasta --keep_original yes

# taxonomy
# This file needs to be downloaded if we use
# NCBI etc data, even a custom database
crabs db_download --source taxonomy

# 2: db_import
# we can import our own custom db
#crabs db_import --input ./processeddata/species_lists/ncbi/20230504_test.fasta --output processeddata/crabs/output.fasta --seq_header species --delim ' '

# 3: db_merge
# when we want to join multiple databases
# tbd

# 4: in silico pcr
# 12S	MiFish-U, FW:GTCGGTAAAACTCGTGCCAGC, R:CATAGTGGGGTATCTAATCCCAGTTTG
crabs insilico_pcr --input mitofishdb.fasta --output pcrout.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --error 4

# 4.2 pga
# uses pairwise global alignment to test if
# any of the references just had the primers trimmed out
# this appends to pcr results
crabs pga --input mitofishdb.fasta --output pgaout.fasta --database pcrout.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --speed medium --percid 0.95 --coverage 0.95 --filter_method strict

# 5 assign_tax
# taxonomy assignment by accession number
# --missing for custom names, same format as lineage input
crabs assign_tax --input pgaout.fasta --output taxassigned.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp

# taxonomy to whole mitofish database
#crabs assign_tax --input mitofishdb.fasta --output mitofishdb_tax.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp

# 6 dereplicate
# options:
# strict: only unique sequences
# single_species: one sequence per species
# uniq_species: all unique sequences for each species
crabs dereplicate --input taxassigned.tsv --output dereplicated.tsv --method uniq_species

# 7 clean database with filters
# min and max length
# --maxns # ambiguous basepairs allowed
# enviro/species: discard environmental/no sp name, yes/no
crabs seq_cleanup --input dereplicated.tsv --output dbcleaned.tsv --discard discardedseq.tsv --minlen 75 --maxlen 1000 --maxns 2 --enviro yes --species yes --nans 0

# crabs db_subset ; can use text file to restrict output

# 8 visualization
# basic initial visuals
# diversity by taxonomic level
# barplot with bar for # species and # sequences
crabs visualization --method diversity --input dbcleaned.tsv --level order

# length of sequences, frequency plot
crabs visualization --method amplicon_length --input dbcleaned.tsv --level order

# db_completeness
# you have to provide the species names in a txt
# gives you summary stats for sp of interest; assesses reference database but also all of NCBI
crabs visualization --method db_completeness --input dbcleaned.tsv --output mifishcomplete.tsv --species species.txt --taxid nodes.dmp --name names.dmp

# primer_efficiency
# makes a bar graph with proportion base pair occurrences
# generates a fasta with sequences that contributed to graph
crabs visualization --method primer_efficiency --input dbcleaned.tsv --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --fwd_name mifish_miya_f --rev_name mifish_miya_r --raw_file mitofishdb.fasta --tax_group Gobiidae --output gobiidae_bind.fasta

# phylo
# not working for me right now but should build a tree
# I think I would rather do this myself tbh
#crabs visualization --method phylo --input dbcleaned.tsv --level family --species species.txt --taxid nodes.dmp --name names.dmp


# 9 export
# exports databse into a format for taxonomy assignment
# in format chose the one you want based on the program you will use
## not needed for this project
#crabs tax_format --input dbcleaned.tsv --output dbvsearch.fasta --format sintax
