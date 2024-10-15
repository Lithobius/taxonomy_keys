#!/bin/bash
#==============================================================================
#title        :  03_crabs_export-blastdb.sh
#description  :  This script outputs the cleaned database into a format that can be used for taxonomic assignment, then processes that format into a local BLAST database. Input - cleaned database from step 02
#author       :  Kate Sheridan
#date         :  20230912
#version      :  0.1.0
#notes        :  requires CRABS, makeblastdb from NCBI.
#==============================================================================

# starts in project root
# move to folder for data download/processing
cd database12s/crabs

pwd

# 9 export
# exports databse into a format for taxonomy assignment
# in format chose the one you want based on the program you will use
## options:
## sintax: >KY305681;tax=d:Eukaryota,p:Chordata,c:Actinopteri,o:Cypriniformes,f:Cyprinidae,g:Barbodes,s:Barbodes_binotatus
## rdp: >KY213961	root;Eukaryota;Chordata;Actinopteri;;Centropomidae;Lates;Lates_japonicus
## qiif: fasta with only accessions, .txt with key to higher taxonomy: KY305681	k__Eukaryota;p__Chordata;c__Actinopteri;o__Cypriniformes;f__Cyprinidae;g__Barbodes;s__Barbodes_binotatus
## dad: >Eukaryota;Chordata;Actinopteri;Cypriniformes;Botiidae;Botia
### Note I don't see accession or species in DAD format
## dads: >KY250421 Oncorhynchus Oncorhynchus_masou
## idt: >Eukaryota;Chordata;Actinopteri;Gadiformes;Gadidae;Gadus;Gadus_macrocephalus
crabs tax_format --input 20241010_dbcleaned.tsv --output blastdb/20241010_12s_crabs.fasta --format rdp

# makeblastdb from ncbi command line tools
# manual and instructions: https://www.ncbi.nlm.nih.gov/books/NBK569861/
makeblastdb -in blastdb/20241010_12s_crabs.fasta -out blastdb/crabs-mifish -title 20241010_crabs-mifish -dbtype nucl -parse_seqids
