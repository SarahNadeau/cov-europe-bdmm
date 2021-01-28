# Plot BDMM migration rates
# OKAY: what is the warning?
#   Seems like duplicate x-values are dropped when ggdistribute interpolates segments to draw the curve of the distribution
#   https://stackoverflow.com/questions/56861001/how-to-suppress-warnings-from-statsregularize-values
# DONE: add color legend
# OKAY: why are distributions all different sizes?
#   ggdistribute was doing something funny in representing the tails of the distribution 
#   I've switched to just plotting using geom_density

require(bdskytools)
require(ggdistribute)
require(dplyr)
require(ggplot2)
require(coda)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/with_hubei_migration_decrease"
LOGFILE <- paste(WORKDIR, "processed_results/combined_chains.log", sep = "/")  # combined logfile
# TRACERTABLE <- paste(WORKDIR, "processed_results/chain_2_tracer_table.txt", sep = "/")  # exported from tracer
BURNIN <- 0  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0
LOGNORM_PRIOR_SD <- 1
OUTDIR <- paste(WORKDIR, "figures", sep = "/")
EUROPEAN_DEMES <- c("France", "Germany", "Italy", "OtherEuropean")
EUROPEAN_DEME_LABELS <- c("France", "Germany", "Italy", "other European")
IS_MIGRATION_RATE_DECREASE <- T  # expect logfile with param prefixes "migrationRateSMEpi.i0_", "migrationRateSMEpi.i1_"
PARAM_PREFIX <- "migrationRateSMEpi" 

system(command = paste("mkdir -p", OUTDIR))

# Load data
logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
# tracer_table <- read.delim(TRACERTABLE)

# Get migration paths
migration_paths <- matrix(nrow = length(EUROPEAN_DEMES) + 1, ncol = length(EUROPEAN_DEMES))
colnames(migration_paths) <- EUROPEAN_DEMES
rownames(migration_paths) <- c("Hubei", EUROPEAN_DEMES) 
for (i in 1:nrow(migration_paths)) {
  for (j in 1:ncol(migration_paths)) {
    source <- c("Hubei", EUROPEAN_DEMES)[i]
    sink <- EUROPEAN_DEMES[j]
    if (source != sink) {
      migration_paths[[i, j]] <- paste(source, sink, sep = "_to_")
    }
  }
}
paths <- c(migration_paths)[!is.na(c(migration_paths))] 

if (IS_MIGRATION_RATE_DECREASE) {
  params <- c(
    paste(PARAM_PREFIX, ".i0_", paths, sep = ""),  # all before Hubei lockdown
    paste(PARAM_PREFIX, ".i1_", paths[grepl(x = paths, pattern = "Hubei_to")], sep = ""))  # from Hubei after lockdown
} else {
  params <- paste(PARAM_PREFIX, paths, sep = "")
}

# Reformat logfile
logfile_long <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "param",
  values_to = "sampled_value")

# Get ESS for each path
# ess_data <- logfile_long %>% 
#   group_by(param) %>%
#   summarize(
#     "ESS_coda" = coda::effectiveSize(sampled_value))
# ess_data$ESS_tracer <- 0
# for (path in ess_data$param) {
#   ess_tracer = as.numeric(as.character(tracer_table[10, path]))
#   ess_data[ess_data$param == path, "ESS_tracer"] <- ess_tracer
# }

get_path_from_param <- function(param, param_prefixes) {
  for (prefix in param_prefixes) {
    if (grepl(x = param, pattern = prefix)) {
      return(gsub(x = param, pattern = prefix, replacement = ""))
    }
  }
  stop(paste("Can't find path from", param))
}

if (IS_MIGRATION_RATE_DECREASE) {
  param_prefixes <- paste(PARAM_PREFIX, c(".i0_", ".i1_"), sep = "")
} else {
  param_prefixes <- PARAM_PREFIX
}
  
logfile_long$path <- unlist(lapply(
  X = logfile_long$param, 
  FUN = get_path_from_param, param_prefixes = param_prefixes))

# ess_data$path <- unlist(lapply(
#   X = ess_data$param,
#   FUN = get_path_from_param, param_prefixes = param_prefixes))

logfile_long <- tidyr::separate(
  data = logfile_long, 
  col = path, 
  into = c("from", "to"),
  sep = "_to_",
  remove = F)

# ess_data <- tidyr::separate(
#   data = ess_data, 
#   col = path, 
#   into = c("from", "to"),
#   sep = "_to_",
#   remove = F)

get_prefix_from_param <- function(param, param_prefixes) {
  for (prefix in param_prefixes) {
    if (grepl(x = param, pattern = prefix)) {
      return(prefix)
    }
  }
  stop(paste("Can't find prefix from", param))
}

logfile_long$prefix <- unlist(lapply(
  X = logfile_long$param,
  FUN = get_prefix_from_param, param_prefixes = param_prefixes))

deme_levels <- c("Hubei", EUROPEAN_DEMES) 
deme_labels <- c("Hubei", EUROPEAN_DEME_LABELS)

logfile_long$to <- factor(
  x = logfile_long$to, levels = deme_levels, labels = paste("sink:", deme_labels))
logfile_long$from <- factor(
  x = logfile_long$from, levels = deme_levels, labels = paste("source:", deme_labels))

# ess_data$to <- factor(
#   x = ess_data$to, levels = deme_levels, labels = paste("sink:", deme_labels))
# ess_data$from <- factor(
#   x = ess_data$from, levels = deme_levels, labels = paste("source:", deme_labels))

logfile_long <- logfile_long[with(logfile_long, order(from, to)), ]

filtered_logfile_long <- logfile_long %>% filter(prefix == param_prefixes[1]) # to make sure I'm generating correct # samples
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), 
    meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD),
  from = filtered_logfile_long$from,
  to = filtered_logfile_long$to)

log_migration_rates_plot <- ggplot() + 
  geom_density(
    data = samples_under_prior,
    aes(x = prior, fill = "prior", alpha = "prior"),
    size = 0) +
  geom_density(
    data = logfile_long %>% filter(prefix == param_prefixes[1]),
    aes(x = sampled_value, fill = "posterior", alpha = "posterior"),
    size = 0) +
  facet_grid(
    from ~ to, 
    scales = "free_y", 
    labeller = label_wrap_gen(10), 
    switch = "y") + 
  theme_bw() + 
  scale_x_continuous(trans='log10', breaks = c(0.1, 1, 10), limits = c(0.05, 10)) + 
  labs(x = "Migration rate") + 
  scale_fill_manual(
    name = element_blank(),
    values = c("prior" = "lightgrey", "posterior" = "#E69F00")) + 
  scale_alpha_manual(
    name = element_blank(),
    values = c("prior" = 0.8, "posterior" = 0.5)) +
  theme(
    legend.position = "bottom", 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())

if (IS_MIGRATION_RATE_DECREASE) {
  log_migration_rates_plot <- log_migration_rates_plot + 
    geom_density(
    data = logfile_long %>% filter(prefix == param_prefixes[2]),
    aes(x = sampled_value, 
        fill = "posterior before Hubei lockdown",
        alpha= "posterior before Hubei lockdown"),
    size = 0) + 
    scale_fill_manual(
      name = element_blank(),
      values = c(
        "prior" = "lightgrey", 
        "posterior" = "#E69F00",
        "posterior before Hubei lockdown" = "#7570B3")) +
    scale_alpha_manual(
      name = element_blank(),
      values = c(
        "prior" = 0.8, 
        "posterior" = 0.5,
        "posterior before Hubei lockdown" = 0.5))
}

png(
  file = paste(OUTDIR, "log_migration_rates.png", sep = "/"),
  width = 6.5,
  height = 5,
  units = "in",
  res = 300)
show(log_migration_rates_plot)
dev.off()

migration_rates_plot <- log_migration_rates_plot + 
  scale_x_continuous(limits = c(0, 10))

png(
  file = paste(OUTDIR, "migration_rates.png", sep = "/"),
  width = 8,
  height = 6,
  units = "in",
  res = 300)
show(migration_rates_plot)
dev.off()

# log_migration_rates_plot <- ggplot() + 
#   ggdistribute::geom_posterior(
#     data = samples_under_prior, 
#     aes(x = prior),
#     draw_ci = F,
#     draw_sd = F,
#     midline = NA) + 
#   ggdistribute::geom_posterior(
#     data = logfile_long %>% filter(prefix == param_prefixes[1]),
#     aes(x = sampled_value),
#     draw_ci = T,
#     draw_sd = F,
#     alpha = 0.5, fill = "#E69F00",
#     interval_type = "hdi",
#     ci_width = 0.95,
#     brighten = 0) +
#   facet_grid(
#     from ~ to, 
#     scales = "free_y", 
#     labeller = label_wrap_gen(10), 
#     switch = "y") + 
#   theme_bw() + 
#   scale_x_continuous(trans='log10') + 
#   scale_y_continuous("Count", position="right") + 
#   labs(x = "Migration rate")
# 
# if (IS_MIGRATION_RATE_DECREASE) {
#   log_migration_rates_plot <- log_migration_rates_plot + 
#     ggdistribute::geom_posterior(
#       data = logfile_long %>% filter(prefix == param_prefixes[2]),
#       aes(x = sampled_value),
#       draw_ci = T,
#       draw_sd = F,
#       alpha = 0.5, fill = "#7570B3",
#       interval_type = "hdi",
#       ci_width = 0.95,
#       brighten = 0)
# }
