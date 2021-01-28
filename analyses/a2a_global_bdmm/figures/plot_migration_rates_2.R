# Plot BDMM migration rates

require(bdskytools)
require(ggdistribute)
require(coda)
require(dplyr)
require(ggplot2)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/a2a_global_bdmm"
OUTPUT1 <- paste(WORKDIR, "figures/log_migration_rates_2.png", sep = "/")
OUTPUT2 <- paste(WORKDIR, "figures/migration_rates_2.png", sep = "/")
LOGFILE <- paste(WORKDIR, "processed_results/combined_chains.log", sep = "/")
BURNIN <- 0  # 10% already excluded by logcombiner
TRACERTABLE <- paste(WORKDIR, "processed_results/tracer_table.txt", sep = "/")
LOGNORM_PRIOR_MEAN <- 0
LOGNORM_PRIOR_SD <- 1

# Load data
logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

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
params <- paste("migrationRateSMEpi.", paths, sep = "")

logfile_long <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "migration_path",
  values_to = "migration_rate")

# Get ESS for each path
ess_data <- logfile_long %>%
  group_by(migration_path) %>%
  summarize(
    "ESS_coda" = coda::effectiveSize(migration_rate))
ess_data$ESS_tracer <- 0
for (path in ess_data$migration_path) {
  ess_tracer = as.numeric(as.character(tracer_table[10, path]))
  ess_data[ess_data$migration_path == path, "ESS_tracer"] <- ess_tracer
}

logfile_long$migration_path <- factor(
  x = logfile_long$migration_path,
  levels = params,
  labels = paths)

ess_data$migration_path <- factor(
  x = ess_data$migration_path,
  levels = params,
  labels = paths)

logfile_long <- tidyr::separate(
  data = logfile_long, 
  col = migration_path, 
  into = c("from", "to"),
  sep = "_to_",
  remove = F)

ess_data <- tidyr::separate(
  data = ess_data,
  col = migration_path,
  into = c("from", "to"),
  sep = "_to_",
  remove = F)

deme_levels <- c("Africa", "AsiaOceania", "Europe", "NorthAmerica", "SouthCentralAmerica")
deme_labels <- c("Africa", "Asia &\nOceania", "Europe", "North America", "South &\nCentral America")

logfile_long$to <- factor(
  x = logfile_long$to, levels = deme_levels, labels = paste("sink:", deme_labels))
logfile_long$from <- factor(
  x = logfile_long$from, levels = deme_levels, labels = paste("source:", deme_labels))

ess_data$to <- factor(
  x = ess_data$to, levels = deme_levels, labels = paste("sink:", deme_labels))
ess_data$from <- factor(
  x = ess_data$from, levels = deme_levels, labels = paste("source:", deme_labels))

logfile_long <- logfile_long[with(logfile_long, order(from, to)), ]
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), 
    meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD),
  from = logfile_long$from,
  to = logfile_long$to)

log_migration_rates_plot <- ggplot(
  data = logfile_long,
  aes(x = migration_rate)) +
  ggdistribute::geom_posterior(
    data = samples_under_prior, 
    aes(x = prior),
    draw_ci = F,
    draw_sd = F,
    midline = NA) + 
  ggdistribute::geom_posterior(
    draw_ci = T,
    draw_sd = F,
    alpha = 0.5, fill = "#E69F00",
    interval_type = "hdi",
    ci_width = 0.95,
    brighten = 0) +
  geom_text(
    data = ess_data,
    aes(
      x = 0, y = Inf, hjust = -0.04, vjust = 1.4,
      label = paste("ESS:", round(ESS_tracer, digits = 0))),
    size = 3) + 
  facet_grid(
    from ~ to, 
    scales = "free_y", 
    labeller = label_wrap_gen(10), 
    switch = "y") +  
  theme_bw() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous("Count", position="right") + 
  labs(x = expression(log[10]*" migration rate"))

png(
  file = OUTPUT1,
  width = 6.5,
  height = 5,
  units = "in",
  res = 300)
show(log_migration_rates_plot)
dev.off()

migration_rates_plot <- ggplot(
  data = logfile_long,
  aes(x = migration_rate)) +
  ggdistribute::geom_posterior(
    draw_ci = T,
    draw_sd = F,
    fill = "blue",
    # midline = "black",
    interval_type = "hdi",
    ci_width = 0.95) +
  geom_text(
    data = ess_data,
    aes(
      x = -Inf, y = Inf, hjust = -0.04, vjust = 1.4,
      label = paste("ESS:", round(ESS_tracer, digits = 0)))) + 
  facet_grid(to ~ from, scales = "free_y", labeller = label_both) + 
  theme_bw() +
  labs(x = paste("migration_rate")) + 
  ggdistribute::geom_posterior(
    data = samples_under_prior, 
    aes(x = prior),
    draw_ci = F,
    draw_sd = F,
    midline = "red",
    interval_type = "hdi",
    ci_width = 0.95,
    alpha = 0.3, fill = "red")

png(
  file = OUTPUT2,
  width = 8,
  height = 6,
  units = "in",
  res = 300)
show(migration_rates_plot)
dev.off()

