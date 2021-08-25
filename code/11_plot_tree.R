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
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## PART 2: Clade information ----

## Read tree:
tree <- read.nexus(file = here("RAxML/Ctenotus.nex"))
tree <- ladderize(tree, right = TRUE)

## Tips that define clades:
t_atla_a <- c("WAMR_146927_Ct_mime", "WAMR_115119_Ct_aust")
t_atla_b <- c("WAMR_104349_Ct_angu", "SAMR_57397_Ct_atla")
t_coll <- c("NA_CCM4867_Ct_stor", "NA_CCM0751_Ct_haly")
t_essi_a <- c("NA_CCM6499_Ct_essi", "NA_CCM2294_Ct_essi")
t_essi_b <- c("ANWC_R05116_Ct_brev", "NA_K073_Ct_essi")
t_inor <- c("NTMR_20378_Ct_robu", "CUMV_14590_Ct_hele")
t_labi_a <- c("WAMR_117161_Ct_cate", "WAMR_135504_Ct_youn")
t_labi_b <- c("SAMR_29480_Ct_impa", "WAMR_144222_Ct_impa")
t_leon <- c("SAMAR_51242_Ct_regi", "WAMR_131785_Ct_leon")
t_pant <- c("WAMR_146582_Ct_rubi", "NA_DLR0038_Ct_pant")
t_scho <- c("WAMR_166437_Ct_scho", "CUMV_14411_Ct_broo")
t_taen <- c("QM_84335_Ct_quin", "SAMAR_33571_Ct_taen")

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

## Labels:
labels <- c("atlas", "colletti", "essingtonii", "inornatus", "labillardieri", "leonhardii", "pantherinus", "schomburgkii", "taeniolatus")

## Prepare tree to plot:
tree_plot <- tree

## PART 3: Plot tree ----

## Plot tree:
plot_tree <- ggtree(tree_plot, color = "gray20", size = 0.4, ladderize = FALSE) + ## Size = branch line thickness.
    
  ## Editing tree tips:
  geom_tiplab(size = 1.5, col = "black", offset = 0.001) +
  scale_color_manual(values = c("gray70", palette)) +
  
  ## Node support:
  #geom_nodelab(size = 1.5, col = "gray50", nudge_x = -0.01, nudge_y = -0.75) +
  
  ## Clade bars:
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[1], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_atla_a)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[1], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_atla_b)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[2], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_coll)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[3], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_essi_a)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[3], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_essi_b)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[4], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_inor)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[5], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_labi_a)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[5], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_labi_b)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[6], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_leon)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[7], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_pant)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[8], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_scho)) +
  geom_cladelabel(barsize = 1, align = TRUE, fontsize = 5, offset.text = 0.005, offset = 0.03, label = paste0(labels[9], " group"), color = "gray20", node = findMRCA(tree = tree_plot, type = "node", tips = t_taen)) +
  
  ## Other edits:
  theme(legend.position = "none")
  
## Plot limits:
plot_t <- plot_tree + xlim(0, 0.37)
#plot_t <- plot_tree  

## Save:
ggsave(plot = plot_t, width = 25, height = 100, units = "cm", limitsize = FALSE, device = "pdf", dpi = 200, filename = here("plots/FigS1.pdf"))
  
## End of script.
