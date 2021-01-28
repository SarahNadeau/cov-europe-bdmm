# Plot BDMM sampling proportions

require(bdskytools)
require(ggdistribute)
require(dplyr)
require(argparse)
require(ggplot2)
require(tidyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/with_hubei_migration_decrease"
OUTPUT1 <- paste(WORKDIR, "figures/sampling_proportion.pdf", sep = "/")
OUTPUT2 <- paste(WORKDIR, "figures/prevelence.pdf", sep = "/")
OUTPUT3 <- paste(WORKDIR, "figures/prevelence.txt", sep = "/")
LOGFILE <- paste(WORKDIR, "processed_results/combined_chains.log", sep = "/")
TRACERTABLE <- paste(WORKDIR, "processed_results/tracer_table.txt", sep = "/")
SAMPLING_DATA <- paste(WORKDIR, "sampling_information.txt", sep = "/")
BURNIN <- 0  # 10% already excluded by logcombiner

demes <- c("Hubei", "France", "Germany", "Italy", "OtherEuropean")
deme_labels <- c("Hubei", "France", "Germany", "Italy", "other European")

# Load data
logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)
sampling_data <- read.delim(file = SAMPLING_DATA)
sampling_prop_upper_bounds <- as.list(sampling_data$sampling_prop_upper_bound)
names(sampling_prop_upper_bounds) <- sampling_data$deme
n_samples <- as.list(sampling_data$n_samples)
names(n_samples) <- sampling_data$deme

# Reformat logfile data
params <- paste("samplingProportionSVEpi.i1_", demes, sep = "")
logfile_long <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "sampling_proportion")

# Generate samples from under the prior for each deme
logfile_long <- logfile_long %>%
  group_by(Deme) %>%
  mutate(Prior = case_when(
    Deme == "samplingProportionSVEpi.i1_France" ~
      runif(n = n(), min = 0, max = sampling_prop_upper_bounds$France),
    Deme == "samplingProportionSVEpi.i1_Germany" ~
      runif(n = n(), min = 0, max = sampling_prop_upper_bounds$Germany),
    Deme == "samplingProportionSVEpi.i1_Hubei" ~
      runif(n = n(), min = 0, max = sampling_prop_upper_bounds$Hubei),
    Deme == "samplingProportionSVEpi.i1_Italy" ~
      runif(n = n(), min = 0, max = sampling_prop_upper_bounds$Italy),
    Deme == "samplingProportionSVEpi.i1_OtherEuropean" ~
      runif(n = n(), min = 0, max = sampling_prop_upper_bounds$OtherEuropean)))

# Transform sampling proportion to prevalence
logfile_long <- logfile_long %>% 
  group_by(Deme) %>%
  mutate(
    prevalence = case_when(
      Deme == "samplingProportionSVEpi.i1_France" ~ 
        n_samples$France / sampling_proportion,
      Deme == "samplingProportionSVEpi.i1_Germany" ~ 
        n_samples$Germany / sampling_proportion,
      Deme == "samplingProportionSVEpi.i1_Hubei" ~ 
        n_samples$Hubei / sampling_proportion,
      Deme == "samplingProportionSVEpi.i1_Italy" ~ 
        n_samples$Italy / sampling_proportion,
      Deme == "samplingProportionSVEpi.i1_OtherEuropean" ~ 
        n_samples$OtherEuropean / sampling_proportion),
    prevalence_prior = case_when(
      Deme == "samplingProportionSVEpi.i1_France" ~
        n_samples$France / Prior,
      Deme == "samplingProportionSVEpi.i1_Germany" ~
        n_samples$Germany / Prior,
      Deme == "samplingProportionSVEpi.i1_Hubei" ~
        n_samples$Hubei / Prior,
      Deme == "samplingProportionSVEpi.i1_Italy" ~
        n_samples$Italy / Prior,
      Deme == "samplingProportionSVEpi.i1_OtherEuropean" ~
        n_samples$OtherEuropean / Prior))

# Clean up deme names for plotting
logfile_long$Deme <- factor(
  x = logfile_long$Deme,
  levels = c(params, "Prior"),
  labels = c(deme_labels, "Prior"))

colnames(sampling_data)[colnames(sampling_data) == "deme"] <- "Deme"
sampling_data$Deme <- factor(
  x = sampling_data$Deme,
  levels = demes,
  labels = deme_labels)
print(sampling_data)

colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
names(colors) <- c("Hubei", "France", "Italy", "Germany", "other European")

# Plot posteriors by Deme
sampling_proportion_plot <- ggplot(
  data = logfile_long,
  aes(x = sampling_proportion, fill = Deme)) +
  ggdistribute::geom_posterior(
    draw_ci = T,
    draw_sd = F,
    interval_type = "hdi",
    ci_width = 0.95,
    brighten = 0,
    cut = 0) +
  geom_vline(
    data = sampling_data, 
    linetype = "dashed",
    aes(xintercept = sampling_prop_upper_bound)) + 
  facet_wrap(Deme ~ ., nrow = 1, scales = "free") +
  theme_bw() + 
  scale_fill_manual(values = colors) +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Sampling proportion", y = "Density")

# sampling_proportion_plot

# pdf(
#   file = OUTPUT1,
#   width = 7,
#   height = 2.5)
# show(sampling_proportion_plot)
# dev.off()

# Add last sample day to deme label
deme_date_labels <- mapply(
  FUN = paste, deme_labels, sampling_data$last_sample_day, sep = "\n")

logfile_long$Deme_last_sample_day <- factor(
  x = logfile_long$Deme,
  levels = c(deme_labels, "Prior"),
  labels = c(deme_date_labels, "Prior"))

colnames(sampling_data)[colnames(sampling_data) == "deme"] <- "Deme"
sampling_data$Deme_last_sample_day <- factor(
    x = sampling_data$Deme,
    levels = deme_labels,
    labels = deme_date_labels)

print(sampling_data)

# Plot posteriors by migration path
prevalence_plot <- ggplot(
  data = logfile_long,
  aes(x = prevalence, fill = Deme)) +
  ggdistribute::geom_posterior(
    draw_ci = T,
    draw_sd = F,
    interval_type = "hdi",
    ci_width = 0.95,
    brighten = 0) +
  # ggdistribute::geom_posterior(
  #   aes(x = prevalence_prior),
  #   fill = "grey", alpha = 0.3,
  #   draw_ci = F,
  #   draw_sd = F,
  #   midline = NA,
  #   brighten = 0) +
  facet_wrap(Deme_last_sample_day ~ ., nrow = 1, scales = "free") + 
  theme_bw() + 
  theme(
    text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90)) +
  geom_vline(
    data = sampling_data, 
    linetype = "dashed",
    aes(xintercept = confirmed_cases_on_last_sample_day)) + 
  labs(x = "Prevalence", y = "Density") + 
  scale_x_continuous(trans='log10') + 
  scale_fill_manual(values = colors)

# prevalence_plot 

# pdf(
#   file = OUTPUT2,
#     width = 7,
#     height = 2.5)
# show(prevalence_plot)
# dev.off()

# Write out 95% HPD for prevalence
prevalence_table <- tracer_table %>%
  select(paste("samplingProportionSVEpi.i1", demes, sep = "_"))
rownames(prevalence_table) <- tracer_table$Summary.Statistic

prevalence_hpd_table <- as.data.frame(t(prevalence_table["95% HPD interval", ])) %>% 
  separate(col = "95% HPD interval", into = c("low", "high"), sep = ", ")
prevalence_hpd_table$low <- unlist(lapply(
  FUN = gsub, 
  X = prevalence_hpd_table$low,
  pattern = "\\[", replacement = ""))
prevalence_hpd_table$high <- unlist(lapply(
  FUN = gsub, 
  X = prevalence_hpd_table$high,
  pattern = "\\]", replacement = ""))
prevalence_hpd_table$low <- as.numeric(prevalence_hpd_table$low)
prevalence_hpd_table$high <- as.numeric(prevalence_hpd_table$high)
prevalence_hpd_table$Deme <- unlist(lapply(
  FUN = gsub,
  X = rownames(prevalence_hpd_table),
  pattern = "samplingProportionSVEpi.i1_", replacement = ""))

prevalence_hpd_table$Deme[prevalence_hpd_table$Deme == "OtherEuropean"] <- "other European"

prevalence_hpd_table_2 <- merge(x = prevalence_hpd_table, y = sampling_data)
prevalence_hpd_table_2 <- prevalence_hpd_table_2 %>% mutate(
  prevalence_high = n_samples / low,
  prevalence_low = n_samples / high)

write.table(
  file = OUTPUT3,
  row.names = F, quote = F, sep = "\t",
  x = prevalence_hpd_table_2[c("Deme", "prevalence_low", "prevalence_high")])
