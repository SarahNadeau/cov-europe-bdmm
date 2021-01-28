# Break "European" clade into global region demes

require(argparse)
require(ape)
require(dplyr)
require(rjson)

# Hardcoded things
# ALIGNMENT <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/results/masked_full_name.fasta"
# METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/data/metadata.tsv"
# CLADES_JSON <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/results/clades.json"
CASE_DATA <- "/Users/nadeaus/Documents/2019-ncov-data/data/metadata/CSSE_covid_19_daily_reports_03-08-2020.csv"
N_EUROPEAN_SEQS <- 100
SEED <- 14617
set.seed(seed = SEED)

# Functions
get_strain_from_sample_name <- function(sample_name) {
  return(unlist(strsplit(sample_name, split = "\\|"))[1])
}
get_id_from_sample_name <- function(sample_name) {
  return(unlist(strsplit(sample_name, split = "\\|"))[2])
}
get_loc_from_sample_name <- function(sample_name) {
  return(unlist(strsplit(sample_name, split = "/"))[1])
}
get_date_from_sample_name <- function(sample_name) {
  return(unlist(strsplit(sample_name, split = "\\|"))[3])
}
get_clade <- function(node) {
  return(node[[1]])
}

summarize_locations <- function(locations) {
  loc_table <- table(locations)
  loc_names <- names(loc_table)
  loc_n_seqs <- unlist(unname(loc_table))
  string <- ""
  for (i in 1:length(loc_names)) {
    if (i < length(loc_names)) {
      string <- paste(string, loc_names[i], " (", loc_n_seqs[i], "), ", sep = "")
    } else {
      string <- paste(string, loc_names[i], " (", loc_n_seqs[i], ")", sep = "")
    }
  }
  return(string)
}

parser <- argparse::ArgumentParser()
parser$add_argument("--alignment", type="character", 
                    help="")
parser$add_argument("--metadata", type="character", 
                    help="")
parser$add_argument("--clades", type = "character",
                    help = "")
parser$add_argument("--output", type = "character",
                    help = "")
parser$add_argument("--selectedseqtable", type = "character",
                    help = "")
parser$add_argument("--demesummarytable", type = "character",
                    help = "")
parser$add_argument("--othereuropesummarytable", type = "character",
                    help = "")
args <- parser$parse_args()

ALIGNMENT <- args$alignment
CLADES_JSON <- args$clades
METADATA <- args$metadata
OUTPUT <- args$output
TABLE_OUTPUT <- args$selectedseqtable
TABLE_OUTPUT2 <- args$demesummarytable
TABLE_OUTPUT3 <- args$othereuropesummarytable

print(paste("alignment:", ALIGNMENT))
print(paste("metadata:", METADATA))
print(paste("clade info:", CLADES_JSON))
print(paste("alignment output:", OUTPUT))
print(paste("table output:", TABLE_OUTPUT))

alignment <- ape::read.FASTA(file = ALIGNMENT)
metadata <- read.delim(file = METADATA)
case_data <- read.delim(file = CASE_DATA, sep = ",")
clade_data <- rjson::fromJSON(file = CLADES_JSON)

clade_data <- data.frame(
  strain = names(clade_data$nodes),
  clade = unlist(lapply(X = clade_data$nodes, FUN = get_clade)))
clade_data <- clade_data[!(grepl(x = clade_data$strain, pattern = "NODE")), ]
# setdiff(seq_info$strain, clade_data$strain) 
# Seqs "Italy/UniMI02/2020"   "Malaysia/186197/2020" are in clade data but not the tree (these leaves were pruned during tree building, see log)
metadata <- merge(x = metadata, y = clade_data, by = "strain")

seq_names <- names(alignment)
seq_info <- data.frame(
  full_seq_name = seq_names,
  strain = unlist(lapply(X = seq_names, FUN = get_strain_from_sample_name)),
  accession_id = unlist(lapply(X = seq_names, FUN = get_id_from_sample_name)),
  location = unlist(lapply(X = seq_names, FUN = get_loc_from_sample_name)),
  collection_date = unlist(lapply(X = seq_names, FUN = get_date_from_sample_name)))
seq_info$collection_date <- as.Date(seq_info$collection_date)
seq_info_merged <- merge(
  x = seq_info, y = metadata[c("gisaid_epi_isl", "region", "age", "country", "division", "country_exposure", "division_exposure", "clade")], 
  by.x = "accession_id", by.y = "gisaid_epi_isl",
  all.x = T, all.y = F)

# Overwrite location from seq name (sometimes this is actually city etc.) with 
# nextstrain country where available
to_replace_filter <- !(is.na(seq_info_merged$country))
seq_info_merged$location <- as.character(seq_info_merged$location)
seq_info_merged$country <- as.character(seq_info_merged$country)
seq_info_merged[to_replace_filter, "location"] <- seq_info_merged$country[to_replace_filter]

# Check that clade A2a corresponds with ORF1b P314L
# aa_alignment <- ape::trans(alignment)
# orf1b_314 <- as.character(as.matrix(aa_alignment))[, 4803]
# orf1b_314_df <- data.frame(
#   full_seq_name = names(orf1b_314),
#   orf1b_314 = orf1b_314)
# seq_info_merged <- merge(x = seq_info_merged, y = orf1b_314_df, by = "full_seq_name")
# table(seq_info_merged[c("clade", "orf1b_314")])
# incongruous1 <- which(seq_info_merged$orf1b_314 == "L" & seq_info_merged$clade != "A2a")
# incongruous2 <- which(seq_info_merged$orf1b_314 == "P" & seq_info_merged$clade == "A2a")
# to_investigate <- seq_info_merged[c(incongruous1, incongruous2), ]
# Netherlands/NA_35/2020 (A1a, L) falls clearly in A1a clade, A1a mutation overrules A2a mutation apparently
# Spain/Madrid_LP10_12/2020 and USA/NY-NYUMC22/2020 fall clearly in A2a clade even w/o the A2a mutations. 
# All these incongruent samples collected after 08-03 so don't have to discuss in the analysis

# Filter to A2a clade before Lombardy lockdown
seq_info_filtered <- seq_info_merged %>% filter(
  clade == "A2a", collection_date <= as.Date("2020-03-08"))

# Assign demes
seq_info_filtered$region <- as.character(seq_info_filtered$region)
seq_info_filtered <- seq_info_filtered %>% mutate(
  deme = case_when(
    region %in% c("Central America", "South America") ~ "SouthCentralAmerica",
    region %in% c("Asia", "Oceania") ~ "AsiaOceania",
    region == "North America" ~ "NorthAmerica",
    T ~ region))

# Manipulate case data so that I can downsample Europe based on # deaths
case_data$Country.Region <- as.character(case_data$Country.Region)
case_data[case_data$Country.Region == "UK", "Country.Region"] <- "United Kingdom"

european_countries <- as.character(unique(
  seq_info_filtered[seq_info_filtered$deme == "Europe", "location"]))

case_data_european <- case_data %>% 
  group_by(Country.Region) %>%
  summarise(n_deaths = sum(Deaths), n_cases = sum(Confirmed)) %>%
  filter(Country.Region %in% european_countries)
all(european_countries %in% case_data_european$Country.Region)

european_countries_n_seqs <- seq_info_filtered %>% 
  group_by(country) %>% 
  filter(country %in% european_countries) %>%
  summarise(n_seqs = n())

case_data_european <- merge(
  x = case_data_european, y = european_countries_n_seqs,
  by.x = "Country.Region", by.y = "country")
case_data_european$padded_n_deaths <- case_data_european$n_deaths + 1

# Take min(# deaths + 1, # available sequences) from each European country
# Italy vastly undersampled, but unfortunately # seqs is limited
selected_samples <- seq_info_filtered %>% filter(deme != "Europe")
for (i in 1:nrow(case_data_european)) {
  country <- case_data_european[i, "Country.Region"]
  n_seqs <- case_data_european[i, "n_seqs"]
  padded_n_deaths <- case_data_european[i, "padded_n_deaths"]
  n_samples <- min(padded_n_deaths, n_seqs)
  candidate_samples <- seq_info_filtered[seq_info_filtered$country == country, ]
  selected_idx <- sample(x = 1:nrow(candidate_samples), size = n_samples)
  selected_samples <- rbind(selected_samples, candidate_samples[selected_idx, ])
}

write.table(
  x = selected_samples,
  file = TABLE_OUTPUT,
  quote = F, row.names = F, sep = "\t")

# Summarize data for each candidate deme
data_summary <- selected_samples %>% 
  group_by(deme) %>%
  summarise(
    "n_seqs" = n(),
    "locations" = summarize_locations(location),
    "first_seq" = min(collection_date),
    "last_seq" = max(collection_date))

write.table(
  x = data_summary,
  file = TABLE_OUTPUT2,
  quote = F, row.names = F, sep = "\t")

data_summary2 <- selected_samples %>% 
  filter(deme == "Europe") %>%
  group_by(country) %>%
  summarise(
    "n_selected_seqs" = n())
data_summary2 <- merge(
  x = data_summary2, y = case_data_european,
  by.x = "country", by.y = "Country.Region", all.y = T)

write.table(
  x = data_summary2,
  file = TABLE_OUTPUT3,
  quote = F, row.names = F, sep = "\t")

print(paste(length(unique(selected_samples$deme)), "demes in total"))
print(paste("Oldest sample that could have been included:", min(seq_info_filtered$collection_date)))
print(paste("Oldest sample included:", min(selected_samples$collection_date)))
print(paste("Youngest sample:", max(selected_samples$collection_date)))

# Create new sample names for analysis: accession_id|deme|collection_date
selected_samples$sample_name <- apply(
  X = selected_samples[, c("accession_id", "deme", "collection_date")],
  MARGIN = 1, 
  FUN = paste0, 
  collapse = "|")
# Remove spaces
selected_samples$sample_name <- gsub(
  x = selected_samples$sample_name,
  pattern = " ", replace = "")

alignment_to_analyse <- subset(
  x = alignment,
  subset = names(alignment) %in% selected_samples$full_seq_name)

old_names <- names(alignment_to_analyse)
new_names <- selected_samples$sample_name[match(old_names, selected_samples$full_seq_name)]
names(alignment_to_analyse) <- new_names

ape::write.FASTA(x = alignment_to_analyse, file = OUTPUT)
