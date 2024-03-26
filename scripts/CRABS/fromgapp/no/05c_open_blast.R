#figuring out how to open blast files

#library(XML)
library(tidyverse)
library(here)

#results <- xmlParse("scripts/python/blast.xml")
#df <- xmlToDataFrame(results)
#list<- xmlToList(results)
#results

#columns are:             from: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
  # Entrez query: taxa ID, and other filters
  # Query ID (ASV###)
  # GI and accession numbers
  # Percent identity
  # alignment length (sequence overlap)
  # mismatch (number of mismatches)
  # gapopen (number of gap openings)
  # qstart (start of allignment in query)
  # qend (end of allignment in query)
  # sstart (start of alignment in subject)
  # send (end of alignment in subject)
  # e-value
  # Bit Score?

#fish in region
fish <- read.csv("./processeddata/species_lists/region/20230504_obis-fish_fishbase_NEP.csv") %>%
  dplyr::select(-c("X", "region")) %>%
  # dplyr::select(c("ncbi_id")) %>%
  distinct() %>%
  filter(!is.na(ncbi_id))
  
fish_1 <- fish %>%
  dplyr::select(c("ncbi_id")) %>%
  distinct() %>%
  filter(!is.na(ncbi_id)) %>%
  .[1:192,] %>%
  as.data.frame() %>%            #    to make a subset for testing
  `colnames<-`("ncbi_id")

t1 <- read_delim("scripts/Python/blastn_all.tab",
                 col_names = FALSE, delim = "\t")  %>%
  `colnames<-`(c("query", "consensus_group", "ID", "percent_identity", "align_length", "mismatch", 
                 "gapopen", "qstart", "qend", "sstart", "send", "e_value", "bit_score")) %>%
  separate(col = ID, into = c("del1", "gi", "del2", "accession", "del3"), sep = "\\|") %>%
  dplyr::select(-c("del1", "del2", "del3")) %>%
  mutate(ncbi_id = parse_number(query), .keep = "unused") %>%
  distinct()

length(unique(t1$ncbi_id))

#there are a few parameters we can use here:
# percent ID, alignment length, and e-value are probably the most useful
t2 <- t1 %>%
  filter(#align_length >= 120, # 12s mifish region is ~170bp
         e_value < 0.01, #
         percent_identity > 60) # can play around with these numbers

length(unique(t2$ncbi_id))
hist(t2$percent_identity)
hist(t2$align_length)


t3 <- t2$ncbi_id

inclusion <- read_file("./scripts/Python/ncbi_IDs.txt") %>%
  strsplit(.,', ') %>%
  as.data.frame() %>%
  .[1:192,] %>% 
  as.data.frame() %>%
  `colnames<-`("query") %>%
  mutate(ncbi_id = parse_number(query), .keep = "unused") %>%
  mutate(inclusion = ifelse(ncbi_id %in% t2$ncbi_id, "y", "n")) %>%
  merge(., fish, by = "ncbi_id") %>%
  distinct()

missing_taxa <- filter(inclusion, inclusion == "n")

trial <- missing_taxa[2:11,] %>%
  rbind(filter(inclusion[1:11,], inclusion == "y"))
write.csv(trial, "./scripts/Python/trial.csv")

t4 <- t1 %>%
  group_by(ncbi_id) %>%
  summarise(max_align = max(align_length),
            min_align = min(align_length),
            max_e_val = max(e_value),
            min_e_val = min(e_value),
            max_ID = max(percent_identity))

trial <- missing_taxa[1:10,] %>%
  rbind(filter(inclusion[1:11,], inclusion == "y")) %>%
  dplyr::select(c(1:3)) %>%
  merge(t4, by = "ncbi_id", all.x = T)
write.csv(trial, "./scripts/Python/trial.csv")

t5 <- filter(t1, ncbi_id %in% trial$ncbi_id) %>%
  merge(., trial, by = "ncbi_id", all.x = T, all.y = F)

ggplot(data = t1, aes(x = e_value, y = align_length)) +
  geom_jitter() +
  scale_x_log10() +
  #  geom_hline(yintercept = 140) +
  #  geom_vline(xintercept = 0.001) +
  theme_classic()

ggplot(data = t5, aes(x = e_value, y = align_length, colour = inclusion)) +
  geom_point() +
  scale_x_log10() +
  #  geom_hline(yintercept = 140) +
  #  geom_vline(xintercept = 0.001) +
  theme_classic() +
  facet_wrap(vars(ncbi_id))


length(unique(missing_taxa$ncbi_id)) # a couple NCBI IDs have multiple taxonomies.


ggplot(data = t1, aes(x = e_value, y = align_length)) +
  geom_jitter() +
  scale_x_log10() +
#  geom_hline(yintercept = 140) +
#  geom_vline(xintercept = 0.001) +
  theme_classic()
  
ggplot(data = t1, aes(x = e_value, y = percent_identity)) +
  geom_point() 


#match against taxa ids going into blast to see for which taxa no blast result was found
#take this and match back to top 500 NCBI IDs to give a Y/N solution

