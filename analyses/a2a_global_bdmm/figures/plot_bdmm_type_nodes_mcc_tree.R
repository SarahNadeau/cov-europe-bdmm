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

WORKDIR <- "/Users/nadeaus/Documents/2019-ncov-data/analyses/2020-04-19_europe_bdmm_prime/a2a_global_bdmm"
OUTPUT <- paste(WORKDIR, "figures/types_on_tree.png", sep = "/")
# NS_METADATA_LINK <- "https://raw.githubusercontent.com/nextstrain/ncov/master/data/metadata.tsv"
NS_METADATA <- "/Users/nadeaus/Documents/2019-ncov-data/data/metadata/nextstrain_metadata_2020-04-01.tsv"
TREEFILE <- paste(WORKDIR, "processed_results/combined_chains_mcc.typed.node.tree", sep = "/")
CLADES_JSON <- "/Users/nadeaus/Documents/2019-ncov-data/data/sequences/2020-04-01/ncov/results/clades.json"

get_clade <- function(node) {
  return(node[[1]])
}

# Load data
tree <- treeio::read.beast(file = TREEFILE)
tree_phylo <- treeio::as.phylo(tree)
# download.file(NS_METADATA_LINK, destfile = NS_METADATA, method = "curl")
ns_info <- read.delim(file = NS_METADATA)
clade_list <- rjson::fromJSON(file = CLADES_JSON)

# Parse data to get clades
clade_data <- data.frame(
  strain = names(clade_list$nodes),
  clade = unlist(lapply(X = clade_list$nodes, FUN = get_clade)))
ns_info <- merge(x = ns_info, y = clade_data, by = "strain")

# Make dataframe with a row for each tip and tip data in columns
tip_data <- tidyr::separate(
  data = data.frame(label = tree_phylo$tip.label), 
  col = "label", 
  into = c("accession_id", "deme", "date"),
  sep = "\\|",
  remove = F)
tip_data$date <- as.Date(tip_data$date)
tip_data <- merge(
  x = tip_data, y = ns_info[c("gisaid_epi_isl", "country", "country_exposure", "clade")],
  by.x = "accession_id", by.y = "gisaid_epi_isl",
  all.x = T, all.y = F)
tip_data <- tip_data[c("label", "accession_id", "deme", "date", "country", "country_exposure", "clade")]
tip_data[c("country", "country_exposure", "clade")] <- apply(
  X = tip_data[c("country", "country_exposure", "clade")], 
  FUN = as.character,
  MARGIN = 2)

# Generate dataframe with 1 row per node, 1 column per type containing type probabilities 
node_data <- treeio::get.data(tree)

type_to_index <- list(
  "Europe" = 1, 
  "AsiaOceania" = 2, 
  "Africa" = 3,
  "SouthCentralAmerica" = 4, 
  "NorthAmerica" = 5)

node_data_2 <- matrix(nrow = 0, ncol = length(type_to_index) + 1)
for (i in 1:nrow(node_data)) {
  node <- node_data[i, "node"]
  types <- c(node_data[[i, "type.set"]])
  probs <- c(node_data[[i, "type.set.prob"]])
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

# Get node number for MRCA of A2a samples, create group to differentiate the clade
tree2 <- full_join(tree, tip_data, by = 'label')
tree2_data <- treeio::get.data(tree2)
A2a_data <- tree2_data %>% filter(clade == "A2a")
A2a_mrca <- MRCA(tree2, A2a_data$node)
tree2 <- groupClade(tree2, .node = A2a_mrca)

# Plot basic tree, highlight A2a clade with different linetype
colors <- c("#00B0F6", "#A3A500", "#F8766D", "#00BF7D", "#E76BF3")
names(colors) <- names(type_to_index)

p <- ggtree(
  tr = tree2, 
  aes(linetype = group), 
  mrsd = max(tip_data$date), as.Date = T) %<+% tip_data + 
  scale_x_date(date_breaks = "1 week", 
               limits = c(as.Date("2020-01-15"), as.Date("2020-03-20"))) + 
  geom_tippoint(aes(color = deme)) + 
  scale_color_manual(values = colors) +
  geom_tiplab(aes(label = country), size = 2, hjust = -0.3) + 
  guides(linetype = F) + 
  theme_tree2() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
show(p)

# Get probabilities for a few nodes in particular
p <- p %<+% node_data_2
p + geom_nodelab(aes(label = node))
# View(node_data_2[node_data_2$node == A2a_mrca, ])

# Add pie charts to tips & internal nodes of tree
pies <- ggtree::nodepie(node_data_2, cols=2:(n_demes + 1))
pies <- lapply(
  X = pies, 
  FUN = function(g) g + scale_fill_manual(values = colors))
p2 <- inset(p, pies, width = 0.1, height = 0.1, vjust = +0.07, hjust = +0.07)

# Add posterior support branch labels in front of pies
p2 <- p2 %<+% node_data +
  geom_nodelab(
    aes(x = branch, 
        label = round(as.numeric(posterior), 2),
        size = 15), 
    geom = "text", size = 2, vjust = -0.5) + 
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(text = element_text(size = 15))

# Get type probabilites at root
print(A2a_mrca)
p3 <- pies$`104`
node_data_2[node_data_2$node == 104, ]

png(file = OUTPUT, width = 10, height = 12, units = "in", res = 300)
multiplot(
  p3 + ggtitle("Root type posterior"), 
  p2 + ggtitle("Maximum clade credibility tree"),
  ncol=2, widths = c(1.5, 4))
dev.off()

