# Plot MCC tree with pie charts at internal nodes representing type 
# probabilities.

require(ape)
require(tidyr)
require(ggtree)
require(treeio)
require(ggplot2)
require(dplyr)
require(ggimage)
require(gridExtra)

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-10-21_european_origins_bdmm_pnas_revision/a2a_global_bdmm"
TREEFILE <- paste(WORKDIR, "processed_results/combined_chains_mcc.typed.node.tree", sep = "/")
ALIGNMENT <- paste(WORKDIR, "a2a_demes.fasta", sep = "/")
METADATA <- paste(WORKDIR, "selected_sample_metadata.txt", sep = "/")
OUTPUT <- paste(WORKDIR, "figures/types_on_tree_2.pdf", sep = "/")
OUTPUT2 <- paste(WORKDIR, "figures/root_type_probabilities.txt", sep = "/")

N_DEMES <- 5

demes <- c("Africa", "AsiaOceania", "Europe", "NorthAmerica", "SouthCentralAmerica")
deme_labels <- c("Africa", "Asia &\nOceania", "Europe", "North\nAmerica", "South &\nCentral\nAmerica")

colors <- RColorBrewer::brewer.pal(n = N_DEMES, name = "Set2")
names(colors) <- c("Africa", "Asia &\nOceania", "Europe", "North\nAmerica", "South &\nCentral\nAmerica")

# Load data
tree <- treeio::read.beast(file = TREEFILE)
tree_phylo <- treeio::as.phylo(tree)
metadata <- read.delim(file = METADATA)
alignment <- ape::read.FASTA(file = ALIGNMENT)
alignment_matrix <- as.character(as.matrix(alignment))

# Confirm all sequences are from A2a clade based on ORF1b:P314L == C14408T 
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
  into = c("accession_id", "deme", "date"), 
  sep = "\\|")
metadata2 <- merge(x = metadata, y = clades_df)

# Make dataframe with a row for each tip and tip data in columns
tip_data <- tidyr::separate(
  data = data.frame(label = tree_phylo$tip.label), 
  col = "label", 
  into = c("accession_id", "deme", "date"),
  sep = "\\|",
  remove = F)
tip_data$date <- as.Date(tip_data$date)
tip_data <- merge(
  x = tip_data, y = metadata2[c("accession_id", "country_exposure", "A2a_clade", "country")])
tip_data <- tip_data[c("label", "accession_id", "deme", "date", "country_exposure", "A2a_clade", "country")]

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
  types <- unlist(node_data[[i, "type.set"]])
  probs <- unlist(node_data[[i, "type.set.prob"]])
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
tip_data$deme <- factor(
  x = tip_data$deme,
  levels = demes,
  labels = deme_labels)

colnames(node_data_2)[2:ncol(node_data_2)] <- 
  deme_labels[match(colnames(node_data_2)[2:ncol(node_data_2)], demes)]

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
  FUN = function(g) g + scale_fill_manual(values = colors))



# Plot tree
p_big <- ggtree(
  tr = tree, 
  mrsd = max(tip_data$date), as.Date = T) %<+% tip_data %<+% node_data + 
  geom_range(
    range='CAheight_0.95_HPD_dates', 
    color='red', 
    alpha=.2, 
    size=1) +
  geom_tippoint(aes(color = deme)) + 
  scale_color_manual(values = colors) +
  scale_x_date(
    date_labels = "%b-%d", 
    limits = c(as.Date("2020-01-10"), as.Date("2020-03-20"))) + 
  guides(linetype = F) + 
  theme_tree2() + 
  geom_tiplab(
    aes(label = country), 
    geom = "text", 
    size = 1.6, hjust = -0.1) + 
  geom_inset(
    insets = internal_pies, 
    width = 0.15, height = 0.15,
    hjust = 0.1, vjust = 0.1) + 
  geom_nodelab(
    aes(x = branch, 
        label = round(as.numeric(posterior), 2)), 
    geom = "text", size = 1.8, vjust = -0.4) + 
  theme(axis.text.x = element_text(size=12),
        legend.text = element_text(size=10)) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p2 <- pies[[as.character(root)]] + scale_fill_manual(values = colors)

pdf(file = OUTPUT, width = 6.5, height = 8)
multiplot(
  p2 + ggtitle("Root type posterior"), 
  p_big + ggtitle("Maximum clade credibility tree"),
  ncol=2, widths = c(1.5, 4))
dev.off()

# Get type probabilites at root
root <- MRCA(tree, tip_data$label)
print(paste("Root is node:", root))
print("Make sure correct pie is being selected - this is hardcoded!!")

# Write out node type probabilities at the root
node_data_of_interest <- node_data_2[node_data_2$node == root, ]
node_data_of_interest$node_desc <- c("root")
write.table(
  x = node_data_of_interest,
  file = OUTPUT2,
  col.names = T, row.names = F, quote = F, sep = "\t")