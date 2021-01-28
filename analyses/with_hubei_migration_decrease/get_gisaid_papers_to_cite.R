require(dplyr)

main_analysis_seqs <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-03-20_to_2020-05-18_european_origins_bdmm/2020-04-19_europe_bdmm_prime/selected_sample_metadata.txt"
a2a_analysis_seqs <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-03-20_to_2020-05-18_european_origins_bdmm/2020-04-19_europe_bdmm_prime/a2a_global_bdmm/selected_sample_metadata.txt"
more_other_european_analysis_seqs <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-03-20_to_2020-05-18_european_origins_bdmm/2020-05-04_bdmm_european_origins/more_other_european_seqs/selected_sample_metadata.txt"
# gisaid_master_dir <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/data/gisaid_2020-05-22"
# ns_master_metadata <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/nextstrain_ncov/data/global_metadata.tsv"
ns_master_metadata <- "/Users/nadeaus/Downloads/metadata_2020-11-10_10-09.tsv"

# Load data
main_analysis <- read.delim(file = main_analysis_seqs)
a2a_analysis <- read.delim(file = a2a_analysis_seqs)
more_other_european_analysis <- read.delim(file = more_other_european_analysis_seqs)

ns_metadata_full <- read.delim(file = ns_master_metadata)

ns_metadata_study <- ns_metadata_full %>% 
  mutate(
    in_main = ifelse(
      test = gisaid_epi_isl %in% main_analysis$accession_id,
      yes = "a", no = ""),
    in_more_eur = ifelse(
      test = gisaid_epi_isl %in% more_other_european_analysis$accession_id,
      yes = "b", no = ""),
    in_a2a = ifelse(
      test = gisaid_epi_isl %in% a2a_analysis$accession_id,
      yes = "c", no = "")) %>%
  unite(col = "analyses", in_main, in_a2a, in_more_eur, sep = ", ", remove = F) %>%
  mutate(analyses = gsub(analyses, pattern = "^, ", replacement = "")) %>%
  mutate(analyses = gsub(analyses, pattern = "^,* ", replacement = "")) %>%
  mutate(analyses = gsub(analyses, pattern = ", , ", replacement = ", ")) %>%
  mutate(analyses = gsub(analyses, pattern = ", $", replacement = "")) %>%
  filter(analyses != "") 

# Check all seqs included
sum(ns_metadata_study$in_main != "") == nrow(main_analysis)
sum(ns_metadata_study$in_a2a != "") == nrow(a2a_analysis)
a2a_analysis$accession_id[!(a2a_analysis$accession_id %in% ns_metadata_full$gisaid_epi_isl)]
sum(ns_metadata_study$in_more_eur != "") == nrow(more_other_european_analysis)
more_other_european_analysis$accession_id[!(more_other_european_analysis$accession_id %in% ns_metadata_full$gisaid_epi_isl)]

# Summarize papers cited
temp <- ns_metadata_study %>% 
  group_by(title) %>%
  summarize(
    n_seqs = n(), 
    authors = paste0(unique(authors), collapse = "; "),
    analyses = paste0(unique(analyses), collapse = "; "))
write.table(x = temp, file = "~/Downloads/papers_cited.txt", quote = F, row.names = F, sep = "\t")  
