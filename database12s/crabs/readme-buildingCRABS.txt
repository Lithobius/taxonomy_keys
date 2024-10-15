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

	NOTE use _if updating if updating to load in the past query and filter by it. Otherwise you will query blast for 10s of thousands of outgroups when you may only need a few thousand.


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

Last run of this took > 1 hour. The FASTA is big. And it looks like it isn't running, but you have to close and open the log to see the updates.

- 1f if-update mitofish pre-ncbi
Input: past-dated taxonomy output
Output: .txt of mitofish accessions trimmed to only new records

The last download of mitofish had 800k records, and these take a long time to search so there's no real need to search all 800k when the information from the last run hasn't changed. This just extracts the NCBI accessions from last time and filters this time by that.

- 1f mitofish-efetch-taxonomy.py
Input: .txt of mitofish accessions
Output: csv of NCBI information for mitofish

Here we will query NCBI for additional information about our accessions so that we can reduce the size and scope of the fasta file. This will only have the accession number and description of the sequence. The description has information about which genes are included so we can filter it in the next step

Note; update this to fix header names and remove desc_verbatim

Last run of this took > 3 hours.

- 1g mitofish-create-filter.Rmd
Input: .csv with accession and description from NCBI
Output: .txt with list of accessions to keep

Mitofish has more than just 12S. If other mitochondrial genes are left in then future steps will be greatly impacted. First here  After we get the "taxonomy" from NCBI, we will classify every sequence by its contents and keep ones that likely contain 12S. These are: containing the words "12S", the whole mitochondrial genome - it will be here no matter what, ribosomal RNA - 12S is ribosomal, and small subunit - this is the portion of the ribosome 12S is on.

- 1h mitofish filter fasta
Input: FASTA of mitofish database, List from 1g
output: FASTA of only Accessions identified by 1g



## 2: Combine all sequences


