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
		- combine all sequences into one FASTA
	- 2 in-silico pcr
	- 3 alignment
	- 4 taxonomic lineage
	- 5 curating with filtering
		- dereplicate
		- cleanup
	- 7 export blast database
	- 6 generate visualizations


Note that the new version of CRABS does not use FASTA as intermediate files. Each later step converts to FASTA from txt before running.

## 1: extract all relevant files and sequences and combine into one FASTA
- 1a crabs_ncbi-12s-fish-shortonly.py
Input: Query
Output: FASTA of fish

NOTE I may be able to do this directly in CRABS now; test this.

First, we are going to run an entrez search using biopython to find 12S fish.

We keep anything under 50,000 basepairs as this is the maximum possible for mtDNA for fish. This means anything in our search that is 'too long' is most likely an annotated homology in a complete chromosome from the nuclear genome. Previous searches have confirmed this, but accession numbers that are found but determined to be too long have been saved into a 'too long' file. It is still possible to get homologous nuclear genes here but there should be fewer.

This automatically formats the fasta for crabs by removing the decimals at the end of the accession numbers.

The query is: "((((((Myxini[Organism]) OR Actinopterygii[Organism]) OR Dipnomorpha[Organism]) OR Chondrichthyes[Organism]) OR Hyperoartia[Organism]) OR Coelacanthimorpha[Organism]) AND 12S"

While many of these will undoubtably be Mitofish duplicates, Mitofish by itself is not reliable, as the way it queries NCBI can result in some versions having different assortments of species. This query was written August 2023 and may require updating to include any fish groups if higher taxonomy is updated at some point.


- 1b crabs_ncbi-12s-outgroups-update-existing-list
Input: Previous FASTAS
Output: CSVs of existing accessions

Since we have run this before, when updating we want to load in the past queries and filter by it. Otherwise you will query blast for 10s of thousands of outgroups when you may only need a few thousand. The output here will filter the results from step 1c

- 1c crabs_ncbi_12s-outgroups-query
Input: results from previous blast searches, CSV of existing accessions
Output: .txt of outgroup search terms

This generates the outgroup query based on previous BLAST searches to NCBI. We search the description of each, where possible, to make sure it does not contain a fish, then collate the accession numbers into a string matching Entrez query format. Then it is filtered by outgroups we already have in our FASTA.

- 1d crabs_ncbi-outgroups-short-only.py
Input: .txt of outgroup search terms
Output: FASTA of outgroups

Similar to 1a, this script uses entrez to query NCBI for outgroups. The output from 1b is the query.

This automatically formats the fasta for crabs by removing the decimals at the end of the accession numbers.

- 1d crabs_mitofish-and-taxonomy.sh
Input: NA
Output: mitofish.fasta, ncbi accession2taxid, names.dmp, nodes.dmp, etc.

This shell script uses CRABS to extract mitofish and the NCBI taxonomy.
Note that this will extract the ENTIRE of both of those databases.

Now there are updated commands for both of these options, script updated 16 Oct 2024 to reflect this

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

I've additionally removed some records that claim to be ok based on my filters but seem to be only one or two genes we aren't interested in. Hopefully this improves the quality of the database overall.

- 1h mitofish filter fasta
Input: FASTA of mitofish database, List from 1g
output: FASTA of only Accessions identified by 1g

Note that this may run in the background, like 1e, so reopen the log to see where it is at. This simply filters the mitofish database by the accessions we identified in 1g.

## 2: Combine all sequences
- 02 import-and-combine.sh
input: all fastas to import and combine
output: combined txt of everything.


Now that we've assembled all we need, its time to combine the sequences into a database to BLAST against.

First we have to import them into crabs format, using the NCBI taxonomy. Then combine what we imported. 



## 3: in silico pcr
Input: merged sequences from 2
output: trimmed sequences to that marker

NOTE you need python 3.9+ at least to avoid a few specific errors

Here we will be running an in-silico pcr. We provide the file as well as the primers, and it will test each sequence forward and reverse and reverse complement for the primers, then trim out and export the portion of the sequence that would amplify.
It also saves the untrimmed sequences as well, i'll assess these later.

You can specify threads or it will autodetect.

## 4: Pairwise global alignment
input: trimmed pcr output, merged sequences from 2
output: aligned sequences appended to pcr output

This catches sequences that might be missed via in silico pcr by aligning them with the amplicons from the in-silico pcr. Then anything with a successful alignment is part of the database.

NOTE: THIS TAKES A LONG TIME. Last run took 4 days 12 hours, and included 'outgroups'.

## 5: filter
Step 1 Dereplicate
input: pga output
output: dereplicated sequences

options for this include: one per species, all unique sequences ignoring taxonomy, and all unique sequences per species. The last one is the recommended option, and what we use.

Step 2 Subset
input: dereplicated sequences
output: subset of those

I don't do this currently, but it can constrain the database to only specific taxa or sequences before export. It might be a good idea to filter out some known 'bad' sequences at some point though.

## 6: export
input: filtered dataset
output: blast db, FASTA, etc.

It is possible to export directly as a blast database, however I had issues with my BLAST not recognizing the output as having taxonomy.

Previously I export as a FASTA for RDP and use makeblastdb on that. Doing this again seemed to work.

There are other formats you can export to, and I'll check back with this later to see if any of them are "better"

## 7: visualizations

CRABS has several visualizations built in.
--diversity-figure gives you a bar chart of how many sequences you have at a specified taxonomic level
