# TODO: fix A2_data <- tree2_data %>% filter(A2a_clade == "A2a" | accession_id == "EPI_ISL_406862") to take GISAID_EPI_ISL_COLNAME (lazyeval puzzle)

require(treeio)
require(ggimage)
require(ggtree)
require(ggplot2)
require(ape)
require(tidyr)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/with_hubei_migration_decrease_low_migration_prior"
TREEFILE <- paste(WORKDIR, "processed_results/combined_chains_mcc.typed.node.tree", sep = "/")
ALIGNMENT <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/old_seqs/europe_demes.fasta"
METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/old_seqs/selected_sample_metadata.txt"
OUTDIR <- paste(WORKDIR, "figures", sep = "/")
NS_CLADES <- "" # paste(WORKDIR, "clades/clades.txt", sep = "/")
GISAID_EPI_ISL_COLNAME <- "accession_id"

OUTPUT <- paste(OUTDIR, "mcc_tree_2.pdf", sep = "/")
OUTPUT2 <- paste(OUTDIR, "a2a_mrca.pdf", sep = "/")
OUTPUT3 <- paste(OUTDIR, "a2_mrca.pdf", sep = "/")
OUTPUT4 <- paste(OUTDIR, "a2a_a2_mrca.txt", sep = "/")

N_DEMES <- 5
DEME_NAMES <- c("Hubei", "France", "Italy", "Germany", "OtherEuropean")
DEME_LABELS <- c("Hubei", "France", "Italy", "Germany", "other European")
DEME_COLORS <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
names(DEME_COLORS) <- DEME_LABELS

system(command = paste("mkdir -p", OUTDIR))

# Load data
tree <- treeio::read.beast(file = TREEFILE)
tree_phylo <- treeio::as.phylo(tree)
metadata <- read.delim(file = METADATA, stringsAsFactors = F)
alignment <- ape::read.FASTA(file = ALIGNMENT)
alignment_matrix <- as.character(as.matrix(alignment))
if (file.exists(NS_CLADES)) {
  plot_clades <- T
  clade_data <- read.delim(NS_CLADES)
}

# Get A2a clade based on ORF1b:P314L == C14408T 
is_A2a <- function(seq, A2a_pos = 14408, A2a_nuc = "T") {
  if (toupper(seq[14408]) == "T") {
    return(T)
  } else {
    return(F)
  }
}

clades <- apply(
  X = alignment_matrix,
  FUN = is_A2a,
  MARGIN = 1)
clades_df <- data.frame(
  seq_name = names(clades),
  A2a_clade = ifelse(test = clades, yes = "A2a", no = "other"))
clades_df <- tidyr::separate(
  data = clades_df,
  col = seq_name,
  into = c(GISAID_EPI_ISL_COLNAME, "deme", "date"), 
  sep = "\\|")
metadata2 <- merge(x = metadata, y = clades_df)

# Make dataframe with a row for each tip and tip data in columns
tip_data <- tidyr::separate(
  data = data.frame(label = tree_phylo$tip.label), 
  col = "label", 
  into = c(GISAID_EPI_ISL_COLNAME, "deme", "date"),
  sep = "\\|",
  remove = F)
tip_data$date <- as.Date(tip_data$date)
tip_data <- merge(
  x = tip_data, y = metadata2[c(GISAID_EPI_ISL_COLNAME, "country", "A2a_clade")])
tip_data <- tip_data[c("label", GISAID_EPI_ISL_COLNAME, "deme", "date", "country", "A2a_clade")]

if (plot_clades) {
  tip_data <- merge(x = tip_data, y = clade_data, by.x = "label", by.y = "name")
}

# Get node number for MRCA of A2a samples, create group to differentiate the clade
tree2 <- full_join(tree, tip_data, by = 'label')
tree2_data <- treeio::get.data(tree2)
A2a_data <- tree2_data[tree2_data$A2a_clade == "A2a" & !is.na(tree2_data$A2a_clade), ]
A2a_mrca <- MRCA(tree2, A2a_data$node)
tree2 <- groupClade(tree2, .node = A2a_mrca)

# Generate dataframe with 1 row per node, 1 column per type containing type probabilities 
node_data <- treeio::get.data(tree)

# Figure out which order deme (type) probabilities are given in
first_node_w_all_types <- which(
  unlist(lapply(X = node_data$type.set, FUN = length)) == N_DEMES)[1]
type_order <- node_data[first_node_w_all_types, "type.set"][[1]][[1]]
type_to_index <- as.list(1:N_DEMES)
names(type_to_index) <- type_order

node_data_2 <- matrix(nrow = 0, ncol = length(type_to_index) + 1)
for (i in 1:nrow(node_data)) {
  node <- node_data[i, "node"]
  types <- unlist(c(node_data[[i, "type.set"]]))
  probs <- unlist(c(node_data[[i, "type.set.prob"]]))
  full_probs <- rep(0, length(type_to_index))
  for (j in 1:length(types)) {
    full_probs[type_to_index[[types[j]]]] <- probs[j]
  }
  data <- c(node, full_probs)
  node_data_2 <- rbind(node_data_2, data)
}
node_data_2 <- as.data.frame(node_data_2)
colnames(node_data_2) <- c("node", names(type_to_index))

node_data_2$node <- as.integer(node_data_2$node)
n_demes <- length(type_to_index)
node_data_2[, 2:(n_demes + 1)] <- apply(X = node_data_2[, 2:(n_demes + 1)], MARGIN = 2, FUN = as.numeric)

# Fix deme order and labels
colnames(node_data_2)[colnames(node_data_2) == "OtherEuropean"] <- "other European"
colnames(node_data)[colnames(node_data) == "OtherEuropean"] <- "other European"
tip_data$deme <- factor(
  x = tip_data$deme,
  levels = DEME_NAMES,
  labels = DEME_LABELS)

# Get node age uncertainty
numeric_date_range_to_date <- function(range, mrsd) {
  date_max <- mrsd - (range[1] * 365)
  date_min <- mrsd - (range[2] * 365)
  return(c(date_min, date_max))
}
node_data$CAheight_0.95_HPD_dates <- lapply(
  X = node_data$CAheight_0.95_HPD,
  FUN = numeric_date_range_to_date,
  mrsd = max(tip_data$date))

# Make pie charts for internal nodes of tree
pies <- ggtree::nodepie(node_data_2, cols=2:(n_demes + 1))
n_tips <- nrow(tip_data)
n_internal_nodes <- n_tips - 1
internal_pies <- pies[(n_tips + 1):(n_tips + n_internal_nodes)]
internal_pies <- lapply(
  X = internal_pies, 
  FUN = function(g) g + scale_fill_manual(values = DEME_COLORS))

# Plot tree
p_big <- ggtree(
  tr = tree2, 
  aes(linetype = group), 
  mrsd = max(tip_data$date), as.Date = T) %<+% tip_data %<+% node_data + 
  geom_range(
    range='CAheight_0.95_HPD_dates', 
    color='red', 
    alpha=.2, 
    size=1) +
  geom_tippoint(aes(color = deme)) + 
  scale_color_manual(values = DEME_COLORS) +
  scale_x_date(
    date_labels = "%b-%d", 
    limits = c(as.Date("2019-11-13"), as.Date("2020-03-20"))) + 
  guides(linetype = F) + 
  theme_tree2() + 
  geom_tiplab(
    aes_string(label = GISAID_EPI_ISL_COLNAME), 
    geom = "text", 
    size = 1.6, hjust = -0.1) + 
  geom_inset(
    insets = internal_pies, 
    width = 0.07, height = 0.07,
    hjust = 0.4, vjust = 0.4) + 
  geom_nodelab(
    aes(x = branch, 
        label = round(as.numeric(posterior), 2)), 
    geom = "text", size = 1.8, vjust = -0.4) + 
  theme(axis.text.x = element_text(size=12),
        legend.text = element_text(size=12)) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)))

if (plot_clades) {
  p_big <- p_big + geom_tippoint(aes(shape = clade))
}

pdf(file = OUTPUT, width = 7, height = 9)
show(p_big)
dev.off()

# Plot A2a and A2 MRCA pie charts larger for annotation of figure 
A2_data <- tree2_data %>% filter(A2a_clade == "A2a" | accession_id == "EPI_ISL_406862")
A2_mrca <- MRCA(tree2, A2_data$node)

pdf(file = OUTPUT2, width = 1, height = 1)
internal_pies[[as.character(A2a_mrca)]]
dev.off()

pdf(file = OUTPUT3, width = 1, height = 1)
internal_pies[[as.character(A2_mrca)]]
dev.off()

# Write out node type probabilities for A2a, A2 MRCA
node_data_of_interest <- rbind(
  node_data_2[node_data_2$node == A2a_mrca, ],
  node_data_2[node_data_2$node == A2_mrca, ])
node_data_of_interest$node_desc <- c("a2a", "a2")
write.table(
  x = node_data_of_interest,
  file = OUTPUT4,
  col.names = T, row.names = F, quote = F, sep = "\t")
