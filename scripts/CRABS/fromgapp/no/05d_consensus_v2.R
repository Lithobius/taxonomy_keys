#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DECIPHER", force = TRUE)
#BiocManager::install("msa", force = TRUE)
library(here)
library(tidyverse)
library(taxonomizr)
library(Biostrings)
library(msa)
library(phytools)



#merge top500 (after updating taxonomy) and fasta, group by family, make allignment, make consensus sequence - make consensus for all

top500_tax <- read_csv('rawdata/top500_20230405/20230413_top500-worrms.csv') 

ASVs <- read.delim("rawdata/top500_20230405/blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out", 
                   h = TRUE, fill = TRUE) %>% 
  # standardize names and remove overhanging x
  janitor::clean_names() %>%
  rename_with(., ~str_remove_all(string = .,
                                 pattern = "x_")) %>%
  separate_wider_delim(taxonomy, delim = ' / ', 
                       names = c('kingdom', 'phylum', 'class', 'order', 
                                 'family', 'genus', 'species')) %>%
  # column to match existing scripts
  dplyr::rename(asv = query_id) %>%
  filter(kingdom == "Eukaryota") %>%
  filter(identity_percentage > 98)

fastaFile <- readDNAStringSet("rawdata/top500_20230405/12S_ASV_sequences.length_var.fasta")
asv = names(fastaFile)
sequence = paste(fastaFile)
asv_sequences <- data.frame(asv, sequence)

#use NCBI_species
t1 <- merge(ASVs[c("asv", "species")],top500_tax, by.x = "species", by.y = "ncbi_species")

t2 <- merge(t1, asv_sequences, by = "asv", all.x = T) %>%
  filter(!is.na(order))
unique(t2$order)

unique(fish$family)

#msaConvert()

df <- data.frame(matrix(ncol = 2, nrow = 0)) %>% `colnames<-`(c("order", "sequence"))

for (i in unique(t2$order)) {
  group <- t2 %>%
    filter(order == i) %>%
    dplyr::select(c("sequence")) %>%
    distinct() %>%
    filter(nchar(sequence) < 190) #something is up here - primers not removed? error in sequencing? - filter out oddities anyways
  consensus = ifelse(nrow(group) == 1,
                     group,
                     msaConsensusSequence(
                       msa(
                         DNAStringSet(group$sequence), 
                         method = "ClustalW"), 
                       type = "Biostrings"))
  df_2 <- data.frame(order = i, sequence = as.character(consensus))
  
  df <- rbind(df, df_2)
}
#writing a fasta file
Xfasta <- character(nrow(df) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", df$order)
Xfasta[c(FALSE, TRUE)] <- df$sequence

writeLines(Xfasta, "rawdata/top500_20230405/top500_consensus.fasta")


#confirming other 12s labelled sequences

fastaFile1 <- readDNAStringSet("rawdata/temp/sequence.fasta")
asv1 = names(fastaFile1)
sequence1 = paste(fastaFile1)
asv_sequences1 <- data.frame(asv = asv1, sequence = sequence1)


n1 <- group[1,]
n2 <- asv_sequences1[1,2]
n3 <- t(merge(n1,n2)) %>% as.data.frame()
group_set1 <- DNAStringSet(n3$V1)
group_aliC1 <- msa(group_set1, method = "ClustalW")
group_aliC1


#alignments and consensus for specific groups
group <- t2 %>%
  filter(order == "Clupeiformes") %>%
  dplyr::select(c("sequence")) %>%
  distinct() %>%
  filter(nchar(sequence) < 180)
group_set <- DNAStringSet(group$sequence)
group_aliC <- msa(group_set, method = "ClustalW")
group_aliC
consensus <- msaConsensusSequence(group_aliC, type = "Biostrings")
consensus
df_2 <- data.frame(order = "Clupeiformes", sequence = as.character(consensus))

df <- rbind(df, df_2)





#mitofish <- readDNAStringSet("rawdata/MitoFish/mito-all.fasta")
#mito_align <- msa(mitofish, method = "ClustalW") 

#consensus <- msaConsensusSequence(alignment, type = "Biostrings")
#consensus

