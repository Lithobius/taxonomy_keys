#!/bin/bash

#==============================================================================
#title        :  02_crabs_combine-mifish-ncbifish-outgroups.sh
#description  :  This script combines the output from querying NCBI and extracting mitofish using CRABS.
#author       :  Kate Sheridan
#date         :  20230912
#version      :  0.1.0
#notes        :  requires CRABS, taxonomy file from 01a.
#==============================================================================

# starts in project root
# move to folder for data download/processing
cd database12s/crabs

# 3: db_merge
# when we want to join multiple databases
crabs db_merge --output 20240223_12s_crabs-merge.fasta --uniq yes --input 20240226_mitofish_12s-filter.fasta 20240212_12s_fish-ncbi-filter.fasta 20240212_12s_outgroups-ncbi-filter2.fasta #20230905_12s_ncbi-outgroups_crabs-ready.fasta

# 4: in silico pcr
# 12S	MiFish-U, FW:GTCGGTAAAACTCGTGCCAGC, R:CATAGTGGGGTATCTAATCCCAGTTTG
crabs insilico_pcr --input 20240223_12s_crabs-merge.fasta --output 20240223_pcrout.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --error 6

# 4.2 pga
# uses pairwise global alignment to test if
# any of the references just had the primers trimmed out
# this appends to pcr results
crabs pga --input 20240223_12s_crabs-merge.fasta --output 20240223_pgaout.fasta --database 20240223_pcrout.fasta --fwd GTCGGTAAAACTCGTGCCAGC --rev CATAGTGGGGTATCTAATCCCAGTTTG --speed medium --percid 0.9 --coverage 0.9 --filter_method relaxed

# 5 assign_tax
# taxonomy assignment by accession number
# --missing for custom names, same format as lineage input
crabs assign_tax --input 20240223_pgaout.fasta --output 20240310_taxassigned.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp

# taxonomy to whole mitofish database
#crabs assign_tax --input mitofishdb.fasta --output mitofishdb_tax.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp

# 6 dereplicate
# options:
# strict: only unique sequences
# single_species: one sequence per species
# uniq_species: all unique sequences for each species
crabs dereplicate --input 20240310_taxassigned.tsv --output 20240310_dereplicated.tsv --method uniq_species

# 7 clean database with filters
# min and max length
# Min, there's some small segments on NCBI that should match if they align
# Max, I have set to around a mitochondrial genome, but I may need to expand
# --maxns # ambiguous basepairs allowed
# enviro: discard environmental? yes
# species: discard 'sp.'? no, we want to permit LCA etc.
# nans: number of missing taxonomic levels, we're going to allow for every level to be nan so that we can get records that have higher taxonomy but not lower taxonomy, or are just missing some of the taxonomy in NCBI but we can get it back with WoRMS. It seems like some of it is just an issue with ncbi taxonomy but isn't necessarily the sequence itself missing information. So keep it all, repopulate in r, bring it back to crabs?
crabs seq_cleanup --input 20240310_dereplicated.tsv --output 20240310_dbcleaned.tsv --discard 20240310_discardedseq.tsv --minlen 20 --maxlen 20000 --maxns 5 --enviro yes --species no --nans 7

# crabs db_subset ; can use text file to restrict output
