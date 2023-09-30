# taxonomy_keys
Central location to deposit current work on keys for taxonomy; revisions, NCBI issues, etc. Then I don't have to search every repository I have for each version. 

## latest-keys
Go here for the latest versions



## Scripts 
To generate keys and lists.
Hakai-rosetta folder with workflow for rosettas. 
Older scripts largely outdated but collected here for reference when reviewing old projects.

Currently (2023-09-27) the workflow needs updated here as we should create a crabs database in this directory.
The latest version was generated from a CRABS database from my chapter 1 repository. 
Right now there are a few possible methods in the folder for extracting NCBI info.
	`00` scripts are to process initial downloads/etc. 
		`extract-taxonomy` would be from a NCBI summary file
		`extract-taxonomy_tinyseq` would be from an NCBI tinyseq file
		`extract_ncbinum` is from the crabs database FASTA
	`01_fish_ncbi-efetch-taxonomy.py` takes the output from `extract-ncbinum` or any txt file that's a list of ncbi numbers to query ncbi for taxonomy using efetch. Currently it does not have a method to handle problems from NCBI and will sometimes throw an error and stop the loop. But there is nothing wrong with the record, it just needs to be restarted.
	`01a_trimlist.py` will trim the list of accession numbers for when efetch-taxonomy hangs up and stops. This way the script can be easily restarted and won't have to rewrite numbers
	`02_fish_ncbi2worms.Rmd` Goes through the output from efetch taxonomy to: indicate gene(s) from accession number, generate a verbatim name to match, match with worms

## Species lists
latest lists

ncbi

BC = lists from partners in BC, OBIS/GBIF pulls filtered to latitudes

Hakai = blast database and rosetta process
	2022rosetta > current draft is where output from scripts go, then are moved to latest-keys. Includes intermediate files, best to use latest-keys instead
	> problems: contains output from the first few scripts in the rosetta process. includes outstanding problems and conflicts currently being explored.
	> taxize : initial taxize output before merging with each other and resolved problems. separated by phylum. >> this means if only one phyla needs updated the whole database doesn't need run through taxize again
	
Fish = right now (august 2022) most lists are in scratchwork as I confirm which ones are valid/useful/not outdated. Mostly copypasted from previous projects

Fossil = 2021 extraction of BC fossils from PaleoBioDB

Invertebrates = same as fish

	
## references

pdfs with species lists for regions

