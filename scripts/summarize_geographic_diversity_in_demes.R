# Investigating possible clusters in data

require(dplyr)

METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/selected_sample_metadata.txt"

metadata <- read.delim(file = METADATA)
geographic_dist <- metadata %>% group_by(deme) %>%
  summarise(n_divisions = length(unique(division)),
            n_seqs = n())

View(metadata[metadata$deme == "Germany", ])
