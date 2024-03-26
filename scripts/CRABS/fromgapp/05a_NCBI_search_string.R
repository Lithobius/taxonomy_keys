#format for NCBI nucleotide search

library(here)
library(tidyverse)
library(taxonomizr)
library(Biostrings)

#fish
fish <- read.csv("./processeddata/species_lists/region/20230504_obis-fish_fishbase_NEP.csv") %>%
  dplyr::select(-c("X", "region")) %>%
 # dplyr::select(c("ncbi_id")) %>%
  distinct()

with_NCBI_IDs <- filter(fish, !is.na(ncbi_id))
without_NCBI_IDs <- filter(fish, is.na(ncbi_id)) #many of these do have NCBI IDs

fish_1 <- fish %>%
  dplyr::select(c("ncbi_id")) %>%
  distinct() %>%
  filter(!is.na(ncbi_id))# %>%
#  .[1:30,] %>%
#  as.data.frame() %>%                to make a subset for testing
#  `colnames<-`("ncbi_id")

a1 <- fish_1 %>%
  mutate(chr = paste("txid", ncbi_id, "[ORGN]", sep = ""))
a2 <- toString(a1$chr)

write_file(a2, "./scripts/Python/ncbi_IDs.txt")


######other stuff ###########

b2 <- fish_1 %>%
  mutate(chr = paste("txid", ncbi_id, "[ORGN]", " AND mitochondrial[WORD] NOT UNVERIFIED[WORD] NOT PREDICTED[WORD]", sep = ""))
b3 <- toString(b2$chr)
write_lines(b3, "./scripts/Python/entrez_queries.txt")


a3 <- gsub(',', ' or',a2)
a4 <- paste("(",a3," AND mitochondrial [WORD] NOT UNVERIFIED [WORD] NOT PREDICTED [WORD] AND 100:20000 [SLEN]")
a4 #paste into search here: https://www.ncbi.nlm.nih.gov/nuccore/advanced



################

#messing around
#test datasets on few sspecies
prop <- readDNAStringSet("./rawdata/NCBI_searches/PROP.fasta")
word <- readDNAStringSet("./rawdata/NCBI_searches/WORD.fasta")
mito <- readDNAStringSet("./rawdata/NCBI_searches/mito.fasta")

seq_name_p = names(prop)
seq_name_w = names(word)
seq_name_m = names(mito)

accessions_p <- sapply(strsplit(seq_name_p, " "), "[", 1)
accessions_w <- sapply(strsplit(seq_name_w, " "), "[", 1)
accessions_m <- sapply(strsplit(seq_name_m, " "), "[", 1)

sequence = paste(prop)
dataframe <- data.frame(accessions_p, seq_name_p, sequence)

t1 <- accessionToTaxa(dataframe$accessions_p)
accessionToTaxa(c("Z17430.1","Z17429.1","X62402.1",'NOTREAL'),sqlFile)






#sequence_p = paste(prop)
#dataframe_p <- data.frame(seq_name, sequence)

accessions_p <- sapply(strsplit(seq_name_p, " "), "[", 1)
accessions_w <- sapply(strsplit(seq_name_w, " "), "[", 1)
accessions_m <- sapply(strsplit(seq_name_m, " "), "[", 1)

pw <- accessions_p[!(accessions_p %in% accessions_w)]
pm <- accessions_p[!(accessions_p %in% accessions_m)]
wm <- accessions_w[!(accessions_w %in% accessions_m)]
mp <- accessions_m[!(accessions_m %in% accessions_p)]









#inverts
inverts <- read.csv('./processeddata/species_lists/region/20230412_obis-invertebrates_nep.csv')%>%
  select(c("ncbi_id")) %>%
  distinct()

b1 <- inverts %>%
  mutate(chr = paste("txid", ncbi_id, "[ORGN]", sep = ""))
b2 <- toString(b1$chr)

b3 <- gsub(',', ' or',b2)
b4 <- paste("(",b3,"AND mitochondrial [WORD] NOT UNVERIFIED [WORD] NOT PREDICTED [WORD] AND 100:20000 [SLEN]")
b4 

