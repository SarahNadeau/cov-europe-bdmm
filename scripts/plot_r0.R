require(bdskytools)
require(ggdistribute)
require(coda)
require(dplyr)
require(ggplot2)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/death_delay_more_other_eur_seqs"
OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "processed_results/combined_chains.log", sep = "/")
TRACERTABLE <- paste(WORKDIR, "processed_results/tracer_table.txt", sep = "/")
BURNIN <- 0  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Hubei", "France", "Germany", "Italy", "OtherEuropean")
deme_labels <- c("Hubei", "France", "Germany", "Italy", "other European")
params <- paste("R0SVEpi.", demes, sep = "")                               

logfile_long <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Get ESS for each deme
ess_data <- logfile_long %>% 
  group_by(Deme) %>%
  summarize(
    "ESS_coda" = coda::effectiveSize(Re))
ess_data$ESS_tracer <- 0
for (deme in ess_data$Deme) {
    ess_tracer = as.numeric(as.character(tracer_table[10, deme]))
    ess_data[ess_data$Deme == deme, "ESS_tracer"] <- ess_tracer
}

# Clean up deme names for plotting
logfile_long$Deme <- factor(
  x = logfile_long$Deme,
  levels = params,
  labels = deme_labels)

ess_data$Deme <- factor(
  x = ess_data$Deme,
  levels = params,
  labels = deme_labels)

colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
names(colors) <- c("Hubei", "France", "Italy", "Germany", "other European")

# Plot R0 posteriors by deme
re_across_demes_plot <- ggplot(
  data = logfile_long,
  aes(x = Re, fill = Deme)) +
  ggdistribute::geom_posterior(
    draw_ci = T,
    draw_sd = F,
    interval_type = "hdi",
    ci_width = 0.95, 
    mirror = T,
    brighten = 0) + 
  theme_bw() +
  coord_flip() + 
  facet_wrap(Deme ~ ., nrow = 1, scale = "free_x") +  
  theme(
    text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none") +
  scale_fill_manual(values = colors) + 
  labs(x = "R", y = "Density")

# show(re_across_demes_plot)

pdf(
  file = OUTPUT,
  width = 7,
  height = 2.5)
show(re_across_demes_plot)
dev.off()

