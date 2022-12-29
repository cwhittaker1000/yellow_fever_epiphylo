# Loading Required Libraries
library(tidyverse); library(ape); library(seqinr); library(ggtree); library(anytime); 
library(phytools); library(ggpmisc); library(tidytree)

# Lots of the below comes from the excellent book: # https://yulab-smu.top/treedata-book/chapter2.html

# Reading in Unrooted Tree 
# tree <- treeio::read.iqtree("3_IQTREE_Unrooted_Phylogenetic_Tree/NCBIVirus_YFVSelectAlignedTrimmedSeqs.fasta.treefile")
tree <- treeio::read.iqtree("4_TempEst_Outlier_IQTREE_Repeat/NCBIVirus_YFVSelectAlignedTrimmedOutlierRemovedSeqs.fasta.treefile")

# Rooting the Tree Using Simple Midpoint Rooting
rooted_tree <- tree 
rooted_tree@phylo <- midpoint.root(tree@phylo) 
ggtree(rooted_tree@phylo)

# Calculating Branch Lengths On This Midpoint Rooted Tree

# rooted_tree@phylo$edge is a dataframe where:
#   2nd column describes the child
#   1st column describes the parent
#   Hence each node only appears once in the 2nd column, 
#   and each node can appear multiple times in 1st column (>2 implies polytomy)
#   and the tips of the tree do NOT appear in 1st column, as they are not parent to any nodes. 
tip_df <- data.frame(label = rooted_tree@phylo$tip.label, node = seq(1:length(rooted_tree@phylo$tip.label)))
oldest_internal_node <- unique(rooted_tree@phylo$edge[, 1][!(rooted_tree@phylo$edge[, 1] %in% rooted_tree@phylo$edge[, 2])])
edge_lengths <- c()
for (i in 1:(oldest_internal_node - 1)) {
  continue <- TRUE
  label <- tip_df$label[i]
  current_node_number <- which(rooted_tree@phylo$tip.label == label)
  temp_edge_lengths <- c()
  while (continue) {
    current_edge_index <- which(rooted_tree@phylo$edge[, 2] == current_node_number)
    temp_edge_lengths <- c(temp_edge_lengths, rooted_tree@phylo$edge.length[current_edge_index])
    current_node_number <- rooted_tree@phylo$edge[current_edge_index, 1]
    if (current_node_number == oldest_internal_node) {
      continue <- FALSE
    }
  }
  edge_lengths <- c(edge_lengths, sum(temp_edge_lengths))
}

# Basic Root-To-Tip Regression
dates <- do.call(rbind, stringr::str_split(rooted_tree@phylo$tip.label, "\\|"))[, 4]
dates <- anydate(dates)
mod2 <- lm(edge_lengths ~ dates) # have an intercept in here as none of the 322 actual genomes are the tMRCA and so should have branch length > 0
summary(mod2)
df <- data.frame(dates = dates, branch_length = edge_lengths)
ggplot(df, aes(dates, branch_length)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  lims(x = c(as.Date("2010-01-01"), as.Date("2022-12-01")))

# Loading in and joining metadata to the phylogenetic tree data 
metadata <- read_tsv("2_Multiple_Sequence_Alignment_Generation/Step2_MAFFT/NCBIVirus_YFVSelectSeqs_Metadata.tsv")
tibble_tree <- as.tibble(rooted_tree)
merged_metadata_tree <- full_join(tibble_tree, metadata, by = c("label" = "taxa"))
rooted_metadata_tree <- as.treedata(merged_metadata_tree)
ggtree(rooted_metadata_tree, mrsd = max(dates)) +
  geom_tippoint(aes(shape = Host, color = Country)) + 
  theme(legend.position = "right") +
  theme_tree2()

# Above is a manually rooted tree - branch lengths aren't in terms of time yet, just number of substitutions
# per site (per year??) I think. So then what we need to do is take the coefficient from the linear regression;
# (which measures root to tip divergence per day) and then multiply it by 365 to get year;
# then changing edge length gives you a rough date (I think)
# See here for explanation: https://biology.stackexchange.com/questions/60841/branch-length-in-phylogenetic-trees
# Note confusion here - when we plot the time-tree along the x-axis, do
# we plot the implied age given the clock, or the *ACTUAL* age based on date
# of sampling (which is known) - ask Filipe this.
rooted_metadata_tree@phylo$edge.length <- rooted_metadata_tree@phylo$edge.length/(coef(mod2)[2] * 365)

ggtree(rooted_metadata_tree, mrsd = max(dates)) +
  geom_tippoint(aes(shape = Host, color = Country)) + 
  theme(legend.position = "right") +
  theme_tree2()

ggtree(rooted_metadata_tree) +
  scale_x_reverse() +
  coord_flip() +
  geom_tippoint(aes(shape = Host, color = Country)) + 
  theme(legend.position = "right") 

# Read the data

# supply a most recent sampling date so you get the dates
# and add a scale bar

# I've gone through in detail and doesn't look like theme_tree2 is using
# any of the @data section
# Perhaps in @phylo??

# Note confusion here - when we plot the time-tree along the x-axis, do
# we plot the implied age given the clock, or the *ACTUAL* age based on date
# of sampling (which is known) - ask Filipe this.
tree <- treeio::read.beast("3_IQTREE_Unrooted_Phylogenetic_Tree/flu_tree_beast.tree")
tree@phylo$edge.length[1] <- 100
ggtree(tree, mrsd="2013-01-01") + 
  theme_tree2() 
# this shows that theme_tree2 is just scaling the edge length based on the most recent sampling 
# date. I.e. edge length needs to be in units of time, not substitutions. 

colnames(tree@data)
colnames(rooted_metadata_tree@data)

# Finally, add tip labels and adjust axis
ggtree(tree, mrsd="2013-01-01") + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5) + 
  xlim(1990, 2020)

##### Old and Scrap Code
# 
# # 1st 539 are the tips
# tree@phylo$tip.label
# 
# # there are 1075 edges that connect all of the 539 tips
# dim(tree@phylo$edge)
# length(tree@phylo$edge.length)
# 
# 
# set.seed(10)
# tree <- rtree(n=5,tip.label=c("A", "B", "C", "D", "E"))
# ggtree(tree, branch.length = "branch.length") +
#   geom_tree() +
#   geom_tiplab() +
#   geom_treescale()
# tips <- "C"
# 
# tip_df <- data.frame(label = tree$tip.label, node = seq(1:5))
# oldest_internal_node <- unique(tree$edge[, 1][!(tree$edge[, 1] %in% tree$edge[, 2])])
# edge_lengths <- c()
# for (i in 1:5) {
#   continue <- TRUE
#   label <- tip_df$label[i]
#   current_node_number <- which(tree$tip.label == label)
#   temp_edge_lengths <- c()
#   while (continue) {
#     current_edge_index <- which(tree$edge[, 2] == current_node_number)
#     temp_edge_lengths <- c(temp_edge_lengths, tree$edge.length[current_edge_index])
#     current_node_number <- tree$edge[current_edge_index, 1]
#     if (current_node_number == oldest_internal_node) {
#       continue <- FALSE
#     }
#   }
#   edge_lengths <- c(edge_lengths, sum(temp_edge_lengths))
#   print(i)
# }
# 
# # Convert BEAST Tree to Tibble, Add Extra Information and Convert Back to Tree
# tree_tibble <- as.tibble(tree)
# label_split <- do.call(rbind.data.frame, str_split(tree_tibble$label[!is.na(tree_tibble$label)], "\\|"))
# colnames(label_split) <- c("Accession", "Country", "Collection_Date")
# for_lab <- tibble(label = tree_tibble$label[!is.na(tree_tibble$label)],
#                   label_split)
# tree_tibble_extra_info <- full_join(tree_tibble, for_lab, by = 'label')
# new_tree <- as.treedata(tree_tibble_extra_info)
# 
# ggtree(new_tree, aes(x, y), size = 0.1, ladderize = TRUE, mrsd = "2021-01-01") +
#   geom_tree() + 
#   scale_x_continuous(limits = c(1940, 2021)) +
#   theme_tree2() +
#   geom_tippoint(size=1, colour="red") +
#   geom_tiplab(aes(y = Country), as_ylab=TRUE, color='firebrick') +
#   theme(axis.text.y = element_text(size = 5)) +
#   geom_rootedge()

## Root-To-Tip Regression for Unrooted Tree

# Calculating Branch Lengths On The Unrooted Tree
# tip_df <- data.frame(label = tree@phylo$tip.label, node = seq(1:length(tree@phylo$tip.label)))
# oldest_internal_node <- unique(tree@phylo$edge[, 1][!(tree@phylo$edge[, 1] %in% tree@phylo$edge[, 2])])
# edge_lengths <- c()
# for (i in 1:(oldest_internal_node - 1)) {
#   continue <- TRUE
#   label <- tip_df$label[i]
#   current_node_number <- which(tree@phylo$tip.label == label)
#   temp_edge_lengths <- c()
#   while (continue) {
#     current_edge_index <- which(tree@phylo$edge[, 2] == current_node_number)
#     temp_edge_lengths <- c(temp_edge_lengths, tree@phylo$edge.length[current_edge_index])
#     current_node_number <- tree@phylo$edge[current_edge_index, 1]
#     if (current_node_number == oldest_internal_node) {
#       continue <- FALSE
#     }
#   }
#   edge_lengths <- c(edge_lengths, sum(temp_edge_lengths))
#   #print(i)
# }
# 
# # Getting Dates and Plotting Tip Cumulative Branch Lengths Against Dates
# ## Note that is for an unrooted tree, and so the branch lengths are meaningless,
# ## and we're not expecting anything particularly realistic looking
# dates <- do.call(rbind, stringr::str_split(tree@phylo$tip.label, "\\|"))[, 4]
# dates <- anydate(dates)
# adj_dates <- (as.numeric(dates) - min(as.numeric(dates)))
# mod <- lm(edge_lengths ~ adj_dates)
# summary(mod)
# plot(dates, edge_lengths)
# lines(dates, predict(mod, dates))