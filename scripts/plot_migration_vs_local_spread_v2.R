# NOTE: ECDC #s negative for some days, resolved by using CSSE data for publication
# NOTE: CSSE data is cumulative cases, so take difference in # cases between days
# NOTE: ECDC data is new cases, ECDC has only China not Hubei specifically

require(dplyr)
require(ggplot2)
require(tidyr)
require(tracerer)

# Hardcoded things
WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/with_hubei_migration_decrease"
SEQ_METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/old_seqs/selected_sample_metadata.txt"
CASE_DATA_FILE_CSSE <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/case_data/csse_time_series_covid19_confirmed_global_2020-06-02.csv"
LOGFILE <- paste(WORKDIR, "processed_results/combined_chains.log", sep = "/")  # combined logfile
OUTDIR <- paste(WORKDIR, "figures", sep = "/")
TRANSMISSION_TO_CONFIRMATION_DELAY <- 5  # days
EUROPEAN_DEMES <- c("France", "Germany", "Italy", "OtherEuropean")
EUROPEAN_DEME_LABELS <- c("France", "Germany", "Italy", "other European")
IS_MIGRATION_RATE_DECREASE <- T  # expect logfile with param prefixes "migrationRateSMEpi.i0_", "migrationRateSMEpi.i1_"
PARAM_PREFIX <- "migrationRateSMEpi" 
MIGRATION_RATE_CHANGETIME <- as.Date("2020-01-23")
date_cutoff <- as.Date("2020-03-08")
colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
names(colors) <- c("Hubei", "France", "Italy", "Germany", "other European")
FIG_PREFIX <- ""

# Load & Manipulate case data
case_data_csse <- read.delim(
  file = CASE_DATA_FILE_CSSE, 
  stringsAsFactors = F, sep = ",",  check.names = FALSE)
seq_metadata <- read.delim(
  file = SEQ_METADATA, sep = "\t", stringsAsFactors = F)

# Massage data
case_data_csse <- case_data_csse %>% mutate(
  country_recoded = recode(
    `Country/Region`, 
    `Czechia` = "Czech Republic"))

# Make sure country names standardized
countries <- unique(seq_metadata$country)
if (length(countries[!(countries %in% unique(case_data_csse$country_recoded))]) > 0) {
  stop("Some country name doesn't match between seq metadata and case data.")
}

case_data_csse_filtered <- case_data_csse %>% 
  filter(
    (country_recoded %in% countries & country_recoded != "China") | 
      `Province/State` == "Hubei") %>% 
  tidyr::pivot_longer(
    cols = colnames(case_data_csse)[5:(ncol(case_data_csse) - 1)],
    names_to = "date",
    values_to = "n_conf_cases") %>% 
  group_by(country_recoded, date) %>%
  summarize(n_conf_cases = sum(n_conf_cases)) %>%
  mutate(deme = case_when(
    country_recoded %in% EUROPEAN_DEMES ~ country_recoded,
    country_recoded == "China" ~ "Hubei",
    T ~ "OtherEuropean"))

# Add time-shift to account for delay between infectivity and confirmation 
case_data_csse_filtered$date <- as.Date(
  x = case_data_csse_filtered$date, format = "%m/%d/%y") - TRANSMISSION_TO_CONFIRMATION_DELAY
min_date <- min(case_data_csse_filtered$date)

# Transform to daily case counts and generated smoothed counts
case_data_csse_demes <- case_data_csse_filtered %>%
  filter(date <= as.Date(date_cutoff)) %>%
  complete(date = seq.Date(min_date, date_cutoff, by = "1 day")) %>%
  group_by(deme, date) %>%
  summarise(n_conf_cases = sum(n_conf_cases)) %>%
  group_by(deme) %>%
  arrange(date) %>%
  mutate(
    n_conf_cases_daily = diff(c(0, n_conf_cases), lag = 1),
    n_conf_cases_daily_smoothed = zoo::rollapply(
      n_conf_cases_daily, 
      width = 7, fill = "extend",
      FUN = mean))

logfile <- tracerer::parse_beast_log(filename = LOGFILE)

tracer_table <- tracerer::calc_summary_stats(
  traces = logfile, sample_interval = 1000)

get_r0 <- function(deme, type) {
  r0_param <- paste("R0SVEpi.", deme, sep = "")
  if (type == "low_95_hpd") {
    return(tracer_table[r0_param, "hpd_interval_low"])
  } else if (type == "high_95_hpd") {
    return(tracer_table[r0_param, "hpd_interval_high"])
  } else {
    return(tracer_table[r0_param, "median"])
  }
}

sink_data <- case_data_csse_demes %>% 
  mutate(
    r0_median = get_r0(deme, type = "median"),
    t =  get_r0(deme, type = "median") * 36.5,
    t_low = get_r0(deme, type = "low_95_hpd") * 36.5,
    t_high = get_r0(deme, type = "high_95_hpd") * 36.5) %>%
  rename("sink_deme" = deme)  %>%
  mutate(
    median = n_conf_cases_daily_smoothed * t,
    hpd_interval_high = n_conf_cases_daily_smoothed * t_high,
    hpd_interval_low = n_conf_cases_daily_smoothed * t_low)

# Calculate 95% bounds for new cases from migration along each path for each date
source_demes <- c("Hubei", EUROPEAN_DEMES)
sink_demes <- EUROPEAN_DEMES
data <- as.data.frame(matrix(nrow = 0, ncol = 7))

get_hpd <- function(samples, incidence) {
  # return daily HPD stats for distribution over daily migration cases
  import_samples <- lapply(X = incidence, samples, FUN = get_hpd_helper)
  return(import_samples)
}

get_hpd_helper <- function(incidence, samples) {
  # multiply scalar daily incidence value by posterior migration rate samples
  # return HPD stats for resulting distribution over daily migration cases
  hpd_interval <- tracerer::calc_hpd_interval(trace = incidence * samples, proportion = 0.95)
  median <- median(x = incidence * samples)
  return(c(hpd_interval, median))
}

for (model_type in c("constant_rates", "hubei_rate_breakpoint")) {
  for (source_deme in source_demes) {
    for (sink_deme in sink_demes) {
      if (source_deme == sink_deme) {
        next
      } else {
        # Get posterior samples over migration rate in boh epochs
        migration_rate_param_i0 <- paste("migrationRateSMEpi.i0_", source_deme, "_to_", sink_deme, sep = "")
        migration_rate_param_i1 <- paste("migrationRateSMEpi.i1_", source_deme, "_to_", sink_deme, sep = "")
        migration_rate_samples_i0 <- logfile[migration_rate_param_i0]
        if (model_type == "constant_rates") {
          migration_rate_samples_i1 <- logfile[migration_rate_param_i0]
        } else {
          migration_rate_samples_i1 <- logfile[migration_rate_param_i1]
        }
        
        # Get daily incidence data in source location for each epoch
        incidence_i0 <- pull(
          case_data_csse_demes[case_data_csse_demes$deme == source_deme & case_data_csse_demes$date <= MIGRATION_RATE_CHANGETIME, ],
          n_conf_cases_daily_smoothed)
        incidence_i1 <- pull(
          case_data_csse_demes[case_data_csse_demes$deme == source_deme & case_data_csse_demes$date > MIGRATION_RATE_CHANGETIME, ],
          n_conf_cases_daily_smoothed)
        
        # Get HPD stats for posterior over daily migration cases
        import_samples_i0 <- as.data.frame(apply(
          X = migration_rate_samples_i0,
          incidence_i0,
          FUN = get_hpd,
          MARGIN = 2))
        import_samples_i1 <- as.data.frame(apply(
          X = migration_rate_samples_i1,
          incidence_i1,
          FUN = get_hpd,
          MARGIN = 2))
        import_samples <- cbind(import_samples_i0, import_samples_i1)
        
        # Re-format HPD stats over posterior of import cases and add to master data frame
        migration_data_temp <- data.frame(
          date = pull(case_data_csse_demes[case_data_csse_demes$deme == source_deme, ], date),
          hpd_interval_low = unlist(import_samples[1, ]),
          hpd_interval_high = unlist(import_samples[2, ]),
          median = unlist(import_samples[3, ]),
          source_deme = source_deme,
          sink_deme = sink_deme,
          model_type = model_type)
        data <- rbind(data, migration_data_temp)
      }
    }
  }
}

migration_data_long <- tidyr::pivot_longer(
  data = data,
  cols = c("hpd_interval_low", "hpd_interval_high", "median"),
  names_to = "interval",
  values_to = "rate") %>% 
  mutate(transmission_type = "Migration")

local_data_long <- sink_data %>% 
  mutate(
    transmission_type = "Within-region transmission",
    interval = "all",
    rate = median)

# Make data frame for regions to shade for having some cases detected


migration_data_long$sink_deme <- factor(
  x = migration_data_long$sink_deme,
  levels = c("Hubei", EUROPEAN_DEMES),
  labels = c("Hubei", EUROPEAN_DEME_LABELS))
migration_data_long$source_deme <- factor(
  x = migration_data_long$source_deme,
  levels = c("Hubei", EUROPEAN_DEMES),
  labels = c("Hubei", EUROPEAN_DEME_LABELS))
local_data_long$sink_deme <- factor(
  x = local_data_long$sink_deme,
  levels = c("Hubei", EUROPEAN_DEMES),
  labels = c("Hubei", EUROPEAN_DEME_LABELS))
migration_data_long$interval <- factor(
  x = migration_data_long$interval,
  levels = c("hpd_interval_high", "median", "hpd_interval_low", "all"),
  labels = c("upper 95% HPD", "median", "lower 95% HPD", "median with 95% HPD error bars"))
local_data_long$interval <- factor(
  x = local_data_long$interval,
  levels = c("hpd_interval_high", "median", "hpd_interval_low", "all"),
  labels = c("upper 95% HPD", "median", "lower 95% HPD", "median with 95% HPD error bars"))

# Make supplemental figures with upper, lower 95% HPD & median
p_supplement <- ggplot( 
  data = migration_data_long %>% filter(model_type == "hubei_rate_breakpoint"),
  aes(x = date, y = rate, fill = source_deme)) +
  geom_bar(stat = "identity") +
  geom_bar(
    data = local_data_long, 
    inherit.aes = F,
    aes(x = date, y = rate, fill = sink_deme),
    stat = "identity") + 
  geom_errorbar(
    data = local_data_long,
    inherit.aes = F,
    alpha = 0.5,
    aes(x = date, ymin = hpd_interval_low, ymax = hpd_interval_high)) +
  facet_grid(
    transmission_type + interval ~ sink_deme, switch = "y", 
    labeller = label_wrap_gen(10),
    scales = "free_y") + 
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.margin = margin(t=0, r=0, b=0, l=0, unit="cm")) +
  scale_x_date(
    breaks = "1 week", date_labels = "%b-%d",
    limits = c(min_date, date_cutoff), expand = c(0, 0)) + 
  labs(x = element_blank(), y = "Rate of new cases") + 
  scale_y_continuous(position="right")
show(p_supplement)

get_date_ranges_with_cases <- function(dates, daily_cases) {
  # Takes a complete date series and the corresponding daily case counts
  # Returns a data frame with columns "start_date" and "end_date" of ranges with cases
  ranges <- matrix(nrow = 0, ncol = 2)
  ranges <- as.data.frame(ranges, stringsAsFactors = F)
  colnames(ranges) <- c("start_date", "end_date")
  prev <- F
  for (i in 1:length(dates)) {
    day_has_cases <- daily_cases[i] > 0
    is_series_end <- i == length(dates)
    if (is_series_end & prev & day_has_cases) {
      end <- dates[i]
      ranges <- rbind(ranges, data.frame(start_date = start, end_date = end))
    } else if (is_series_end & !prev & day_has_cases) {
      date <- dates[i] 
      ranges <- rbind(ranges, data.frame(start_date = date, end_date = date))
    } else if (!prev & day_has_cases) {
      start <- dates[i]
      prev <- T
    } else if (prev & !day_has_cases) {
      end <- dates[i - 1]
      ranges <- rbind(ranges, data.frame(start_date = start, end_date = end))
      prev <- F
    }
  }
  return(ranges)
}

is_first <- T
for (deme in unique(local_data_long$sink_deme)) {
  deme_data <- local_data_long %>% filter(sink_deme == deme) %>% arrange(date)
  deme_ranges <- get_date_ranges_with_cases(
    dates = deme_data$date, daily_cases = deme_data$n_conf_cases_daily)
  deme_ranges$sink_deme <- deme
  deme_ranges$transmission_type <- "Within-region transmission"
  if (is_first) {
    conf_case_range_data <- deme_ranges
    is_first <- F
  } else {
    conf_case_range_data <- rbind(conf_case_range_data, deme_ranges)
  }
}
conf_case_range_data$sink_deme <- factor(
  x = conf_case_range_data$sink_deme,
  levels = c("Hubei", EUROPEAN_DEME_LABELS))

# Make main text figure with median only
p_main_text <- ggplot(
  data = migration_data_long %>% filter(interval == "median", model_type == "hubei_rate_breakpoint"),
  aes(x = date, y = rate, fill = source_deme)) +
  geom_bar(stat = "identity") +
  geom_bar(
    data = local_data_long,
    aes(x = date, y = median, fill = sink_deme),
    stat = "identity") +
  facet_grid(transmission_type ~ sink_deme, switch = "y", 
             labeller = label_wrap_gen(10),
             scales = "free_y") + 
  scale_fill_manual(values = colors) + 
  theme_bw() + 
  geom_rect(
    data = conf_case_range_data,
    inherit.aes = F,
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    alpha = 0.2,
    fill = "black") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.margin = margin(t=0, r=0, b=0, l=0, unit="cm"),
    panel.grid.minor = element_line(size = 0.3), 
    panel.grid.major = element_line(size = 0.3)) +
  scale_x_date(
    breaks = "1 week", date_labels = "%b-%d",
    limits = c(min_date, date_cutoff), expand = c(0, 0)) + 
  labs(x = element_blank(), y = "Rate of new cases") +
  scale_y_continuous(position="right", expand = expand_scale(mult = c(0, .1)))
show(p_main_text)

# Save plots
pdf(
  file = paste(OUTDIR, paste(FIG_PREFIX, "transmission_rate_comparison_main.pdf", sep = ""), sep = "/"),
  width = 7,
  height = 3.5)
show(p_main_text)
dev.off()
pdf(
  file = paste(OUTDIR, paste(FIG_PREFIX, "transmission_rate_comparison_supplemental.pdf", sep = ""), sep = "/"),
  width = 7,
  height = 4.5)
show(p_supplement + theme(legend.position = "none"))
dev.off()
