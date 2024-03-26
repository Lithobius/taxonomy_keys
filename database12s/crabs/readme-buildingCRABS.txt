# README

Here's the instructions for building or updating a 12S CRABS database.

While CRABS itself is fairly simple, I have taken extra steps in building the database to ensure it contains sufficient outgroups and contains as many fish as possible.

The pipeline is as follows:
	- 1 extract all relevant files and sequences to put into the database
		- extract fish from ncbi
		- extract outgroups from ncbi
		- extract mitofish
			- remove non 12s sequences from mitofish
		- extract ncbi taxonomy file
	- 2 combine all sequences
		- merge into one fasta
		- in-silico pcr
		- alignment
		- dereplicate
		- cleanup
	- 3 export blast database
	- 4 generate visualizations


## 1: extract all relevant files and sequences to put into the database
- 1a crabs_ncbi-12s-fish-shortonly.py
Input: Query
Output: FASTA of fish

First, we are going to run an entrez search using biopython to find 12S fish.

We keep anything under 50,000 basepairs as this is the maximum possible for mtDNA for fish. This means anything in our search that is 'too long' is most likely an annotated homology in a complete chromosome from the nuclear genome. Previous searches have confirmed this, but accession numbers that are found but determined to be too long have been saved into a 'too long' file. It is still possible to get homologous nuclear genes here but there should be fewer.

This automatically formats the fasta for crabs by removing the decimals at the end of the accession numbers.

The query is: "((((((Myxini[Organism]) OR Actinopterygii[Organism]) OR Dipnomorpha[Organism]) OR Chondrichthyes[Organism]) OR Hyperoartia[Organism]) OR Coelacanthimorpha[Organism]) AND 12S"

While many of these will undoubtably be Mitofish duplicates, Mitofish by itself is not reliable, as the way it queries NCBI can result in some versions having different assortments of species. This query was written August 2023 and may require updating to include any fish groups if higher taxonomy is updated at some point.

- 1b crabs_ncbi_12s-outgroups-query
Input: results from previous blast searches
Output: .txt of outgroup search terms

This generates the outgroup query based on previous BLAST searches to NCBI. We search the description of each to make sure it does not contain a fish, then collate the accession numbers into a string matching Entrez query format.


- 1c crabs_ncbi-outgroups-short-only.py
Input: .txt of outgroup search terms
Output: FASTA of outgroups

Similar to 1a, this script uses entrez to query NCBI for outgroups. The output from 1b is the query.

This automatically formats the fasta for crabs by removing the decimals at the end of the accession numbers.

- 1d crabs_mitofish-and-taxonomy.sh
Input: NA
Output: mitofish.fasta, ncbi accession2taxid, names.dmp, nodes.dmp, etc.

This shell script uses CRABS to extract mitofish and the NCBI taxonomy.
Note that this will extract the ENTIRE of both of those databases.

- 1e mitofish-extract_ncbinum.py
Input: mitofishdb.fasta
Output: .txt of mitofish accessions

The mitofish database contains a lot of non-12S sequences. This creates problems during the stage of database creation that attempts to align sequences, and takes up unnecessary space. To create a filter to get 12s, we will query the NCBI accessions from this fasta file to extract their descriptions and filter them. 

This script reads the fasta, saves the fasta names, and then writes them to a .txt file as a list separated by linebreaks. Note that if for whatever reason the fasta does not have plain accession numbers, the get_terms_fasta function must be updated to extract only the accession.

Last run of this took > 1 hour.

- 1f mitofish-efetch-taxonomy.py
Input: .txt of mitofish accessions
Output: csv of NCBI information for mitofish

Here we will query NCBI for additional information about our accessions so that we can reduce the size and scope of the fasta file. This will only have the accession number and description of the sequence. The description has information about which genes are included so we can filter it in the next step

Note; update this to fix header names and remove desc_verbatim

Last run of this took > 3 hours.

- 1g mitofish-create-filter.Rmd
Input: .csv with accession and description from NCBI
Output: .txt with list of accessions to keep

- 1h mitofish

## 2: Combine all sequences


