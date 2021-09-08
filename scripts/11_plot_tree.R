### The goals of this script are:
### To plot a phylogenetic tree showing the relationships between Ctenotus samples.
### Checked for functionality on April 23rd 2021.
### Written by Ivan Prates (ivanprates.org).

### PART 1: Getting ready ----

## Packages:
library(ggtree)
library(phytools)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/" ## Rhinellax.
setwd(paste0(path, "spheno_ddRAD/Ctenotus_species_cohesion/"))
detach("package:here", unload = TRUE)
library(here)

## PART 2: Clade information ----

## Sample information:
sample_info <- read.csv("sample_information/ddRAD_sample_information.csv", header = TRUE)

## Labels:
sample_info$UPDATED_SP <- gsub(sample_info$UPDATED_SP, pattern = "cf. ", replacement = "cf")
sample_info$ID_label <- paste0("paste(\"", sample_info$SAMPLE_ID)
sample_info$ID_label <- gsub(sample_info$ID_label, pattern = "_Er_[a-z]+", replacement = "_E.")
sample_info$ID_label <- gsub(sample_info$ID_label, pattern = "_Le_[a-z]+", replacement = "_L.")
sample_info$ID_label <- gsub(sample_info$ID_label, pattern = "_Ct_[a-z]+", replacement = "_C.")
sample_info$ID_label <- gsub(sample_info$ID_label, pattern = "_ct_[a-z]+", replacement = "_C.")
sample_info$ID_label <- gsub(sample_info$ID_label, pattern = "NA_([A-Z]+)([0-9]+)", replacement = "\\1_\\2")
sample_info$ID_label <- gsub(sample_info$ID_label, pattern = "_", replacement = " ")
sample_info$tip_label <- paste0(sample_info$ID_label, " ", sample_info$UPDATED_SP)
sample_info$tip_label <- gsub(sample_info$tip_label, pattern = "(C. [a-z]+)", replacement = ", italic\\(\" \\1\"\\)\\)")
sample_info$tip_label <- gsub(sample_info$tip_label, pattern = "(L. [a-z]+)", replacement = ", italic\\(\" \\1\"\\)\\)")
sample_info$tip_label <- gsub(sample_info$tip_label, pattern = "(E. [a-z]+)", replacement = ", italic\\(\" \\1\"\\)\\)")
sample_info$tip_label <- gsub(sample_info$tip_label, pattern = " ,", replacement = "\",")
sample_info$tip_label <- gsub(sample_info$tip_label, pattern = "cf", replacement = "cf. ")

## Read tree:
tree <- read.nexus(file = "RAxML/Ctenotus.nex")
tree <- ladderize(tree, right = TRUE)

## Change tree trips:
tree$tip.label <- sample_info$tip_label

## Tips that define clades:
t_atla_a <- c("paste(\"WAMR 146927\", italic(\" C. mimetes\"))", 
              "paste(\"WAMR 115119\", italic(\" C. australis\"))")
t_atla_b <- c("paste(\"WAMR 104349\", italic(\" C. angusticeps\"))", 
              "paste(\"SAMR 57397\", italic(\" C. atlas\"))")
t_coll <- c("paste(\"CCM 4867\", italic(\" C. cf. storri\"))", 
            "paste(\"CCM 0751\", italic(\" C. cf. halysis\"))")
t_essi_a <- c("paste(\"CCM 6499\", italic(\" C. cf. coggeri\"))", 
              "paste(\"CCM 2294\", italic(\" C. cf. coggeri\"))")
t_essi_b <- c("paste(\"ANWC R05116\", italic(\" C. brevipes\"))", 
              "paste(\"K 073\", italic(\" C. cf. essingtonii\"))")
t_inor <- c("paste(\"NTMR 20378\", italic(\" C. robustus\"))", 
            "paste(\"CUMV 14590\", italic(\" C. inornatus\"))")
t_labi_a <- c("paste(\"WAMR 117161\", italic(\" C. catenifer\"))", 
              "paste(\"WAMR 135504\", italic(\" C. youngsoni\"))")
t_labi_b <- c("paste(\"SAMR 29480\", italic(\" C. impar\"))", 
              "paste(\"WAMR 144222\", italic(\" C. impar\"))")
t_leon <- c("paste(\"SAMAR 51242\", italic(\" C. regius\"))", 
            "paste(\"WAMR 131785\", italic(\" C. leonhardii\"))")
t_pant <- c("paste(\"WAMR 146582\", italic(\" C. rubicundus\"))", 
            "paste(\"DLR 0038\", italic(\" C. cf. pantherinus\"))")
t_scho <- c("paste(\"WAMR 166437\", italic(\" C. schomburgkii\"))", 
            "paste(\"CUMV 14411\", italic(\" C. brooksi\"))")
t_taen <- c("paste(\"QM 84335\", italic(\" C. quinkan\"))", 
            "paste(\"SAMAR 33571\", italic(\" C. taeniolatus\"))")

## Clades:
c_atla_a <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_atla_a))$tip.label
c_atla_b <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_atla_b))$tip.label
c_coll <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_coll))$tip.label
c_essi_a <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_essi_a))$tip.label
c_essi_b <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_essi_b))$tip.label
c_inor <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_inor))$tip.label
c_labi_a <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_labi_a))$tip.label
c_labi_b <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_labi_b))$tip.label
c_leon <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_leon))$tip.label
c_pant <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_pant))$tip.label
c_scho <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_scho))$tip.label
c_taen <- extract.clade(tree, node = findMRCA(tree = tree, type = "node", tips = t_taen))$tip.label

## List of samples in each taxon:
clade_list <- list(c_atla_a, c_atla_b, c_coll, c_essi_a, c_essi_b, c_inor, c_labi_a, c_labi_b, c_leon, c_pant, c_scho, c_taen)

## Group tips by taxon to color tips in ggtree. 
names(clade_list) <- c("a atla a", "b atla b", "c coll", "d essi a", "e essi b", "f inor", 
                       "g labi a", "h labi b", "i leon", "j pant", "k scho", "l taen") ## To preserve order.
tree <- groupOTU(tree, clade_list)

## Group labels:
labels <- c("atlas", "colletti", "essingtonii", "inornatus", "labillardieri", "leonhardii", "pantherinus", "schomburgkii", "taeniolatus")
labels <- paste0("italic(\"", labels, "\")")
labels <- paste0("paste(", labels, ", \" group\")")

## Prepare tree to plot:
tree_plot <- tree

## PART 3: Plot tree ----

## Plot tree:
plot_tree <- ggtree(tree_plot, color = "gray20", size = 0.4, ladderize = FALSE) + ## Size = branch line thickness.
    
  ## Editing tree tips:
  geom_tiplab(size = 1.5, col = "black", offset = 0.001, parse = TRUE) +
  scale_color_manual(values = c("gray70", palette)) +
  
  ## Node support:
  #geom_nodelab(size = 1.5, col = "gray50", nudge_x = -0.01, nudge_y = -0.75) +
  
  ## Clade bars:
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[1], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_atla_a)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[1], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_atla_b)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[2], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_coll)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[3], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_essi_a)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[3], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_essi_b)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[4], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_inor)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[5], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_labi_a)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[5], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_labi_b)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[6], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_leon)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[7], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_pant)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[8], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_scho)) +
   geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.035, label = labels[9], parse = TRUE, color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_taen)) +
  
  ## Other edits:
  theme(legend.position = "none")
  
## Plot limits:
plot_t <- plot_tree + xlim(0, 0.37)

## Save:
ggsave(plot = plot_t, width = 26, height = 100, units = "cm", limitsize = FALSE, device = "pdf", dpi = 200, filename = here("plots/FigS1_test.pdf"))
  
## End of script.
