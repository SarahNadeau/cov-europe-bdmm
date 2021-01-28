# Make table of number of introductions along each migration path, compare to 
# known introductions from metadata

require(bdskytools)
require(ggdistribute)
require(coda)
require(dplyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/a2a_global_bdmm/figures/"
OUTPUT <- paste(WORKDIR, "n_introductions.txt", sep = "/")
TRACERTABLE <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/a2a_global_bdmm/processed_results/tracer_table.txt"  # manually saved from tracer
NS_METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/data/metadata/nextstrain_metadata_2020-04-01.tsv"

# Load data
tracer_table <- read.delim(TRACERTABLE)
metadata <- read.delim(NS_METADATA)

# Reformat logfile data
paths <- c(
  "Africa_to_AsiaOceania",               
  "Africa_to_Europe",                     
  "Africa_to_NorthAmerica",               
  "Africa_to_SouthCentralAmerica",        
  "AsiaOceania_to_Africa",                
  "AsiaOceania_to_Europe",                
  "AsiaOceania_to_NorthAmerica",          
  "AsiaOceania_to_SouthCentralAmerica",   
  "Europe_to_Africa",                     
  "Europe_to_AsiaOceania",                
  "Europe_to_NorthAmerica",               
  "Europe_to_SouthCentralAmerica",        
  "NorthAmerica_to_Africa",               
  "NorthAmerica_to_AsiaOceania",          
  "NorthAmerica_to_Europe",               
  "NorthAmerica_to_SouthCentralAmerica",  
  "SouthCentralAmerica_to_Africa",        
  "SouthCentralAmerica_to_AsiaOceania",   
  "SouthCentralAmerica_to_Europe",        
  "SouthCentralAmerica_to_NorthAmerica")
params <- paste("typeMappedTree.count_", paths, sep = "")

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

# Clean up deme names for plotting
param_labels <- paths

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
migration_res_table <- migration_res_table[c("from", "Africa", "AsiaOceania", "Europe", "NorthAmerica", "SouthCentralAmerica")]

write.table(
  x = migration_res_table,
  file = OUTPUT,
  sep = "\t",
  quote = F, row.names = F, col.names = T)
