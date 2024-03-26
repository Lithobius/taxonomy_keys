#!/bin/bash

# starts in project root
# move to folder for data download/processing
cd processeddata/crabs/mitofish_nep

# taxonomy
# This file needs to be downloaded if we use
# NCBI etc data, even a custom database
crabs db_download --source taxonomy

# mitofish db is just all of mitofish
crabs db_download --source mitofish --output mitofishdb.fasta --keep_original yes

# 2: db_import
# we can import our own custom db
#crabs db_import --input 20230504_test.fasta --output dbimport.fasta --seq_header species --delim ' '

# 3 db db_merge
# add db merge step for mitofish + custom

# 4: in silico pcr
# 12S	MiFish-U, FW:GTCGGTAAAACTCGTGCCAGC, R:CATAGTGGGGTATCTAATCCCAGTTTG
crabs insilico_pcr --input mitofishdb.fasta --output mitofish_pcrout.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --error 3

# 4.2 pga
# uses pairwise global alignment to test if
# any of the references just had the primers trimmed out
# this appends to pcr results
# note: --input is original db, --database is PCR results
crabs pga --input mitofishdb.fasta --output mitofish_pgaout.fasta --database mitofish_pcrout.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --speed medium --percid 0.95 --coverage 0.95 --filter_method strict

# 5 assign_tax
# taxonomy assignment by accession number
# --missing for custom names, same format as lineage input
## note content from db_import may need to be updated first
## open a github issue later
crabs assign_tax --input mitofish_pgaout.fasta --output mitofish_taxassigned.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp

# 6 dereplicate
# options:
# strict: only unique sequences
# single_species: one sequence per species
# uniq_species: all unique sequences for each species
crabs dereplicate --input mitofish_taxassigned.tsv --output mitofish_dereplicated.tsv --method uniq_species


# 7 clean database with filters
# min and max length
# --maxns # ambiguous basepairs allowed
# enviro/species: discard environmental/no sp name, yes/no
crabs seq_cleanup --input mitofish_dereplicated.tsv --output mitofish_dbcleaned.tsv --discard mitofish_discardedseq.tsv --minlen 25 --maxlen 1000 --maxns 2 --enviro yes --species yes --nans 0

# 7.2 crabs db_subset ; can use text file to restrict output
# put input of species in region here

# 8 visualization
# basic initial visuals
# diversity by taxonomic level
# barplot with bar for # species and # sequences
crabs visualization --method diversity --input mitofish_dbcleaned.tsv --level order

# length of sequences, frequency plot
crabs visualization --method amplicon_length --input dbcleaned.tsv --level order

# db_completeness
# you have to provide the species names in a txt
# gives you summary stats for sp of interest; assesses reference database but also all of NCBI
crabs visualization --method db_completeness --input mitofish_dbcleaned.tsv --output mifishcompleteness.tsv --species 20230511_fishsp.txt --taxid nodes.dmp --name names.dmp

# primer_efficiency
# makes a bar graph with proportion base pair occurrences
# generates a tsv with sequences that contributed to graph
crabs visualization --method primer_efficiency --input mitofish_dbcleaned.tsv --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --fwd_name mifish_miya_f --rev_name mifish_miya_r --raw_file mitofishdb.fasta --tax_group Gobiidae --output Gobiidae_bind.tsv
