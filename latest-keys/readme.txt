Latest Keys readme

These are the current drafts, organized by date.

Please notify me of errors or missing data asap with a csv of verbatim ids.

How to use:
1) verbatim_identification (COI) or species_verbatim (MiFish) to match blast hit
2) replace with found_taxa
3) filter by freswater_only/terrestrial_only == 0, 
4) select only found_taxa, verbatim_identification, rank and higher taxonomy
5) rename verbatim_identification to match your dataset
6) merge back with eDNA data
7) remove verbatim_id if desired, rename found_taxa 

THERE IS NO SPECIES COLUMN (August 2022)
It's already a big file, R is slow with big files, and species is redundant with found_taxa in over 50% (COI) or 90% (MiFish) of entries.
Either a) filter to species
filter(rank == 'species') or 
b) make your own if desired
mutate(species = ifelse(rank == 'species', found_taxa, NA_character_))



COI:
	[date]_coi_rosetta_draft.csv
		contains the most recent compilation of the taxonomy key for the hakai database
	Contents:
			verbatim_identification should match BLAST hit
			identification is the verbatim_identification without numbers and 'sp.' etc. In some cases this query may modify the verbatim ID so future comparison may be important.
			found_taxa : replacement from either open tree of life or worms
			rank : what taxonomic level it is to. COI starts at Phylum
			worms_aphiaid, bold_id, gbif_key, ott_id: identification number for that taxa in various databases. You only need this if you will be querying those databases again.
			status + unacceptreason: WoRMS output accepted or not and what reason is provided. Informative in most situations where found_taxa and identification don't match.
			is_synonym + flags: Open tree of life output. Not always up to date, but informative if using output for OTL
			valid_authority: WoRMS output, publication of type specimen
			worms_sciname, tol_sciname, bold_taxon, gbif_sciname: useful if querying those databases again, seeing where names differ. Can use as replacements but warning that it not all taxa are present in all databases so many NAs, also advise using id number to requery higher taxonomy. Do not need to retain otherwise.
			freshwater_only, terrestrial_only: 1 = true 0 = false. Discard taxa with 1 from findings. Compiled generally from WoRMS output (worrms, wm_records_names()), but manually for species not found in WoRMS, sometimes at higher taxonomy than species.
			
			
MIFISH: 

(note august 2022: there are a few minor errors still being fixed, an intermediate version with updates for esa exists for a few taxa in the haida gwaii dataset, but the current version is still dated 0729; formatting etc on the intermediate is wrong, its missing columns, etc. Also some fish got true for terrestrial only and I'm going to figure out why...)

The biggest difference between MiFish and COI files is that Mifish uses species_verbatim instead of verbatim_identification. This is because entries were all to species in the beginning; since then I have added some higher taxonomy.

1)	[date]_mifish_all_draft.csv
		contains EVERYTHING in the Midori 12S classifier. Including bacteria, protists, plants...
		Just enough to sort the taxonomy, not a full key.
		
2)	[date]_mifish_fish-only_rosetta_draft.csv
		Contains only fish, defined as class == Actinopteri, Coelocanthi, Dipneusti, Elasmobranchii, Holocephali, Leptocardii, Myxini, and Petromyzontii
		NOTE that in other databases than WoRMS these classes may appear as other ranks
	contents:
		species_verbatim: match to ID from classifier
		found_taxa is the match from either WoRMS or Open Tree of Life
		taxa_query: modified species_verbatim to search databases with
		other than that columns are the same as COI. 
	

FRESHWATER_TERRESTRIAL
	[date]_freshwater_terrestrial.csv
		this is extracted from the COI rosetta. Technically all the information is in the COI rosetta file but it is technically a smaller file. From a previous workflow, may become deprecated in the future.

	
REMOVE_FROM_SPLISTS

These are in process lists of specific problem taxa / frequent substitutions. Right now it is not a substitution list and actively being adjusted

The sequence for the ASV in question should just be re-blasted if its either not_bc or not_east_pacific.
Not_nearshore is for when we're trying to filter out any taxa that may be true hits but aren't going to be found in the nearshore
other_problem is the key field here; species that appear here are in the problems_above_my_paygrade.csv or possible_taxonomic_revision_key in in_process_spreadsheets. 
	

	
	
	