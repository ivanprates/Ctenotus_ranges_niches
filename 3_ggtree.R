library(ape)
library(castor)
library(cowplot)
library(ggtree)
library(here)
library(magrittr)
library(tidyverse)

## PART 1: Defining putative clusters.

## First creating an empty list to populate with qmatrices for each K:
qmatrix_tree_list <- vector("list", length(minK:maxK))

## Looping over K values:
for (K in minK:maxK){

  ## Selecting inferred qmatrix for corresponding K:
  qmatrix_tree <- sNMF_results$qmatrix_list[[K-minK+1]]
  
  ## Keeping only columns of interest:
  qmatrix_tree <- qmatrix_tree[c("ID", "cluster_assign")]
  
  ## Selecting list of sNMF clusters for K:
  putative_taxa <- putative_taxa_list[[K-minK+1]]

  ## Change levels from clusters to putative taxa:
  levels(qmatrix_tree$cluster_assign) <- putative_taxa

  ## Function will return:
  qmatrix_tree_list[[K-minK+1]] <- qmatrix_tree
  
  } ## Close loop.

## Merging data frames in that resulting list using "reduce" from the purrr package:
assignments <- reduce(qmatrix_tree_list, left_join, by = "ID")
names(assignments) <- c("ID", "K13", "K14", "K15")

## Create new column to include final assigments:
assignments$decision <- NA

## Now changing manually the assigment of certain samples:
assignments[12:13, 5] <- "essingtonii 1a"
assignments[14:16, 5] <- "essingtonii 1b"
assignments[29:30, 5] <- "robustus 1"
assignments[48:61, 5] <- "spaldingi 1"
assignments[62, 5] <- "spaldingi 2a"
assignments[63:85, 5] <- "spaldingi 2b"
assignments[99:106, 5] <- "superciliaris 1"
assignments[107:126, 5] <- "superciliaris 2"
  
## For the remaining samples, replacing rows with the values from another column:
assignments$decision[is.na(assignments$decision)] <- as.character(assignments$K13)[is.na(assignments$decision)]

## Changing levels to ensure their order as in the tree (otherwise R will make it alphabetical):
assignments$decision <- factor(assignments$decision, levels = unique(assignments$decision))

## Splitting this dataframe into a list by putative taxon:
decision <- split(x = assignments, f = assignments$decision)

## PART 2: Getting node number to annotate corresponding clades in phylogeny:

## Loading phylogenetic tree:
full_tree <- read.nexus(paste0(here(),"/RAxML/ni277/RAxML_bipartitions.nex"))
full_tree <- ladderize(full_tree) ## Ladderize tree.

## Assigning colors to tree terminals. For now we're going to set them all as "black".
color_tips <- data.frame(ID = full_tree$tip.label, color = "black")

## Now listing the individuals that were not included in sNMF analyses:
tips_not_in_sNMF <- full_tree$tip.label[(full_tree$tip.label %in% qmatrix_tree_list[[1]]$ID) == FALSE]

## Change the "color" column to character, and then change color of samples not included in sNMF to red:
color_tips$color <- as.character(color_tips$color)
color_tips$color[color_tips$ID %in% tips_not_in_sNMF] <- "#FF3333"

## Now changing the color of outgroups to gray:
color_tips[1:11, "color"] <- "blue"

## Defining color scheme for clade labels on tree plot:
#clade_colors <- palette_list[[3]] ## For 15 clusters.
clade_colors <- c("#474747", "#A4A4A4", "#DFB9B0", "#EE7555", "#EE7555", "#F58F4D", "#F7B649", "#F0E33B",
                  "#A2DA3C", "#4EBF4A", "#2F9975", "#2F9975", "#275578", "#310B6C", "#6E0689", "#9D1484", "#C03265")

## Creating an empty list to populate with node numbers:
cluster_nodes <- vector("list", length(decision))

## Change name of each element in this list:
names(cluster_nodes) <- names(decision)

## Recording the MRCA of each cluster over a loop:
for (n in 1:length(decision)){ 
  
  ## Find MRCA node of all samples corresponding to a species and record in the list we created before:
  cluster_nodes[[n]] <- get_mrca_of_set(tree = full_tree, descendants = as.vector(decision[[n]]$ID))

  } ## Close loop (node values).

## PART 3: Plotting tree.

## Creating new tree plot object with ggtree:
tree_plot <- ggtree(full_tree, size = 0.35, ladderize = FALSE) ## Size controls branch line thickness.

## Editing tree tips:associated taxa
tree_plot <- tree_plot + geom_tiplab(size = 4, color = color_tips$color) ## Change tip font size, color.

## Getting rid of empty margins around plot. This is not working very well and wide margins remain...
tree_plot <- tree_plot + theme(plot.margin = unit(c(t = -2, r = 0, b = -2, l = 0), "in"))
#labs(x = NULL, y = NULL)

## Flip nodes of to approximate non-sister clades that correspond to the same sNMF cluster:
tree_plot <- flip(tree_plot, node1 = 378, node2 = 379) ## essingtonii 1b and essingtonii  2.
tree_plot <- flip(tree_plot, node1 = 314, node2 = 319) ## spaldingi 2a and rimacola.

## Adding bars to each clade over a loop:
for (n in 1:length(cluster_nodes)){ 
  
  ## Adding bars: 
  tree_plot <- tree_plot + geom_cladelabel(node = cluster_nodes[[n]], label = names(cluster_nodes)[[n]], 
                                           angle = 0, fontsize = 6, offset = 0.0015, offset.text = 0.0001,
                                           color = clade_colors[[n]], barsize = 2.5)
  
  } ## Close loop (nodes).

## Check tree plot:
tree_plot

## Saving plot:
save_plot(plot = tree_plot, base_width = 25, base_height = 50, limitsize = FALSE, 
          filename = paste0(here(), "/ggtree/RAxML_n277_putative_clusters.jpeg"))

  ## End of script.
