# Make table of number of introductions along each migration path, compare to 
# known introductions from metadata

# TODO: use tracerer package rather than tracer table

require(tidyr)
require(dplyr)
require(tracerer)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/with_hubei_migration_decrease"
OUTPUT <- paste(WORKDIR, "figures/n_introductions.txt", sep = "/")
LOGFILE <- paste(WORKDIR, "processed_results/combined_chains.log", sep = "/")  # combined logfile
# BURNIN <- 0  # 10% already excluded by logcombiner
TRACERTABLE <- paste(WORKDIR, "processed_results/tracer_table.txt", sep = "/")
METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-03-20_to_2020-05-18_european_origins_bdmm/2020-04-19_europe_bdmm_prime/selected_sample_metadata.txt"

# Load data
tracer_table <- read.delim(TRACERTABLE)
metadata <- read.delim(METADATA)

# Reformat tracertable data
param_labels <- c(
  "France_to_Germany",
  "France_to_Italy",
  "France_to_OtherEuropean",
  "Germany_to_France",
  "Germany_to_Italy",
  "Germany_to_OtherEuropean",
  "Italy_to_France",
  "Italy_to_Germany",
  "Italy_to_OtherEuropean",
  "OtherEuropean_to_France",
  "OtherEuropean_to_Germany",
  "OtherEuropean_to_Italy",
  "Hubei_to_France",
  "Hubei_to_Germany",
  "Hubei_to_Italy",
  "Hubei_to_OtherEuropean")
params <- paste("typeMappedTree.count_", param_labels, sep = "")

# Table median n_introductions
migration_res <- tracer_table[c("Summary.Statistic", params)]

migration_res_long <- tidyr::pivot_longer(
  data = migration_res,
  cols = params,
  names_to = "migration_path")
migration_res_wide <- tidyr::pivot_wider(
  data = migration_res_long,
  names_from = "Summary.Statistic",
  values_from = "value")

migration_res_wide$migration_path <- factor(
  x = migration_res_wide$migration_path,
  levels = params,
  labels = param_labels)

migration_res_wide <- tidyr::separate(
  data = migration_res_wide,
  col = migration_path,
  into = c("from", "to"),
  sep = "_to_")

migration_res_wide$result <- apply(
  X = migration_res_wide[c("median", "95% HPD interval")],
  MARGIN = 1, 
  FUN = paste0, collapse = ": ")

migration_res_table <- tidyr::pivot_wider(
  data = migration_res_wide[c("from", "to", "result")],
  names_from = "to",
  values_from = "result")
migration_res_table <- migration_res_table[c("from", "France", "Germany", "Italy", "OtherEuropean")]

write.table(
  x = migration_res_table,
  file = OUTPUT,
  sep = "\t",
  quote = F, row.names = F, col.names = T)

# Check against known travel history
metadata[c("country", "country_exposure")] <- apply(
  X = metadata[c("country", "country_exposure")],
  FUN = as.character,
  MARGIN = 2)

travel_cases <- metadata %>%
  filter(country != country_exposure) %>% 
  mutate(from = country_exposure, 
         to = country) %>% 
  select(c(from, to))
