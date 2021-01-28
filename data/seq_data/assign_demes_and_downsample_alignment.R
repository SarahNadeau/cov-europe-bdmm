CASE_DATA <- "/Users/nadeaus/Documents/2019-ncov-data/data/metadata/CSSE_covid_19_daily_reports_03-08-2020.csv"
DEMES_OF_INTEREST <- c("Italy", "Germany", "France")
MAX_SEQS_PER_DEME <- Inf
# N_OTHER_EUROPE_SEQS <- 50
N_HUBEI_SEQS <- 10
MIN_SEQS <- 0
SEED <- 127

set.seed(seed = SEED)

# ALIGNMENT <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/results/masked_full_name.fasta"
# METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/data/metadata.tsv"
# CLADES_JSON <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/results/clades.json"

require(argparse)
require(ape)
require(dplyr)
require(rjson)

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
metadata <- merge(x = metadata, y = clade_data, by = "strain")

seq_names <- names(alignment)
seq_info <- data.frame(
  full_seq_name = seq_names,
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

# Filter to only samples from China and Europe
# Take only Hubei samples before Jan 23 lockdown of Wuhan 
# Take only European samples before Mar 8 lockdown of Lombardy
# (so that can assume single Re regimes)
seq_info_filtered <- seq_info_merged %>%
  filter(
    (division == "Hubei" &  collection_date <= as.Date("2020-01-23")) | 
      (region == "Europe") & collection_date <= as.Date("2020-03-08"))

# Filter to countries with at least MIN_SEQS available
countries_to_keep <- seq_info_filtered %>%
  group_by(location) %>%
  summarise("n_seqs" = n()) %>%
  filter(n_seqs >= MIN_SEQS)
countries_to_keep <- countries_to_keep$location
seq_info_filtered <- seq_info_filtered[seq_info_filtered$location %in% countries_to_keep, ]

# Make sure all sequences have an assigned deme
print(paste(sum(is.na(seq_info_filtered$location)), "sequences aren't assigned a location"))

# Assign demes
seq_info_filtered <- seq_info_filtered %>%
  mutate(deme = case_when(
    division == "Hubei" ~ "Hubei",
    location %in% DEMES_OF_INTEREST ~ location,
    T ~ "OtherEuropean"))

# Manipulate case data so that I can downsample other europe proportional to # deaths
case_data$Country.Region <- as.character(case_data$Country.Region)
case_data[case_data$Country.Region == "UK", "Country.Region"] <- "United Kingdom"

other_european_countries <- as.character(unique(seq_info_filtered[seq_info_filtered$deme == "OtherEuropean", "location"]))

case_data_other_european <- case_data %>% 
  group_by(Country.Region) %>%
  summarise(n_deaths = sum(Deaths), n_cases = sum(Confirmed)) %>%
  filter(Country.Region %in% other_european_countries)

other_european_countries_n_seqs <- seq_info_filtered %>% 
  group_by(country) %>% 
  filter(country %in% other_european_countries) %>%
  summarise(n_seqs = n())

case_data_other_european <- merge(
  x = case_data_other_european, y = other_european_countries_n_seqs,
  by.x = "Country.Region", by.y = "country")

# Pad n_deaths by 1 so all countries represented
case_data_other_european$padded_n_deaths <- case_data_other_european$n_deaths + 1

# Make sure I can take as many seqs as there were deaths
if (any(case_data_other_european$n_seqs < case_data_other_european$padded_n_deaths)) {
  stop("Not enough sequences from at least 1 country to select as many as # deaths")
}

# Downsample valid sequences so that BDMM runs in reasonable time
selected_samples <- matrix(nrow = 0, ncol = ncol(seq_info_filtered))

for (deme in unique(seq_info_filtered$deme)) {
  deme_samples <- seq_info_filtered[seq_info_filtered$deme == deme, ]
  count = 0
  if (deme == "Hubei") {
    # Subsample
    ids <- sample(deme_samples$accession_id, size = N_HUBEI_SEQS)
    selected_samples <- rbind(selected_samples, seq_info_filtered[seq_info_filtered$accession_id %in% ids, ])
  } else if (deme == "OtherEuropean") {
    # Sample 1 seq for each death in a country (1 seq for countries w/ 0 deaths)
    for (country in unique(deme_samples$country)) {
      country_samples <- deme_samples[deme_samples$country == country, ]
      n_samples_to_select <- case_data_other_european[case_data_other_european$Country.Region == country, "padded_n_deaths"]
      ids <- sample(country_samples$accession_id, size = n_samples_to_select)
      selected_samples <- rbind(selected_samples, seq_info_filtered[seq_info_filtered$accession_id %in% ids, ])
    }
  } else {
    # Take all seqs from demes of interest
    selected_samples <- rbind(selected_samples, deme_samples) 
  }
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
    "locations" = paste0(unique(location), collapse = ", "),
    "first_seq" = min(collection_date),
    "last_seq" = max(collection_date))

write.table(
  x = data_summary,
  file = TABLE_OUTPUT2,
  quote = F, row.names = F, sep = "\t")

data_summary2 <- selected_samples %>% 
  group_by(country) %>%
  summarise(
    "n_selected_seqs" = n())
data_summary2 <- merge(
  x = data_summary2, y = case_data_other_european,
  by.x = "country", by.y = "Country.Region")

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
