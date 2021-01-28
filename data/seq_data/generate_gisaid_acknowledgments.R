main_analysis_seqs <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/selected_sample_metadata.txt"
a2a_analysis_seqs <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/a2a_global_bdmm/selected_sample_metadata.txt"
more_other_european_analysis_seqs <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-05-04_bdmm_european_origins/more_other_european_seqs/selected_sample_metadata.txt"
# gisaid_master_dir <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/data/gisaid_2020-05-22"
ns_master_metadata <- "/Users/nadeaus/Documents/2019-ncov-data/CH_sequencing/nextstrain_ncov/data/global_metadata.tsv"

OUTPUT <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-05-18_european_origins/gisaid_acknowledgments.txt"

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
  filter(analyses != "") %>%
  select(in_main, in_a2a, in_more_eur, analyses, strain, virus, gisaid_epi_isl, date, country, division, location, country_exposure, originating_lab, submitting_lab)

# Check all seqs included
sum(ns_metadata_study$in_main != "") == nrow(main_analysis)
sum(ns_metadata_study$in_a2a != "") == nrow(a2a_analysis)
sum(ns_metadata_study$in_more_eur != "") == nrow(more_other_european_analysis)

ns_metadata_study <- ns_metadata_study %>%
  select(strain, gisaid_epi_isl, date, country, originating_lab, submitting_lab, analyses) %>%
  rename(
    "Analyses used in" = analyses,
    "Nextstrain sample name" = strain,
    "GISAID Accession ID" = gisaid_epi_isl,
    "Date" = date,
    "Country" = country,
    "Originating lab" = originating_lab,
    "Submitting lab" = submitting_lab)

write.table(
  x = ns_metadata_study,
  file = OUTPUT,
  row.names = F, col.names = T,
  quote = F, sep = "\t")

# first <- T
# for (file in list.files(path = gisaid_master_dir, full.names = T, pattern = ".tsv")) {
#   print(file)
#   if (first) {
#     first <- F
#     gisaid_metadata_full <- read.delim(file, skip = 2)
#   } else {
#     gisaid_metadata_full <- rbind(gisaid_metadata_full, read.delim(file, skip = 2))
#   }
# }

# gisaid_metadata_study <- gisaid_metadata_full %>% 
#   mutate(
#     in_main = ifelse(
#       test = `Accession.ID` %in% main_analysis$accession_id, 
#       yes = "Main analysis", no = ""),
#     in_a2a = ifelse(
#       test = `Accession.ID` %in% a2a_analysis$accession_id, 
#       yes = "A2a clade global", no = ""),
#     in_more_eur = ifelse(
#       test = `Accession.ID` %in% more_other_european_analysis$accession_id,
#       yes = "Death delay other European downsampling", no = "")) %>%
#   unite(col = "Analyses", in_main, in_a2a, in_more_eur, sep = ", ") %>%
#   mutate(Analyses = gsub(Analyses, pattern = "^, ", replacement = "")) %>%
#   mutate(Analyses = gsub(Analyses, pattern = "^,* ", replacement = "")) %>%
#   mutate(Analyses = gsub(Analyses, pattern = ", , ", replacement = ", ")) %>%
#   mutate(Analyses = gsub(Analyses, pattern = ", $", replacement = "")) %>%
#   filter(Analyses != "") %>%
#   select(`Virus.name`, `Accession.ID`, `Collection.date`, Location, )
