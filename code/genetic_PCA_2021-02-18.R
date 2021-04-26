####################
### The goals of this script are:
### To run genetic PCAs based on SNP data to support OTU delimitation.
### To make plots of the first few PC axes.
### Checked for functionality on April 23rd 2021.

## PART 1: Getting ready ----

## Packages:
library(LEA)
library(pals)
library(patchwork)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## Creating directories to save results:
dir.create(path = here("sNMF_runs/PCA"))

## PART 2: Function: Running PCA on the SNPs ----

## Function:
run_PCA <- function(group, miss_ind, miss_SNP) {
  
  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5
  
  ## Location of geno file:
  geno <- here(paste0("sNMF_runs/data/", group, "_mi", miss_ind, ".geno"))
  
  ## Run PCA analysis:
  PCA <- pca(geno)
  
  ## Displaying information on the PCA analysis:
  #summary(PCA)
  
  ## Plotting eigenvalues:
  #plot(PCA, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")
  
  ## Storing the proportion of variation represented by selected PC axes"
  PC_perc <- vector("list", 4) # We'll use the four first PC axes and save them into this list.
  for (p in c(1:4)) { PC_perc[[p]] <- round((summary(PCA)[2, p]*100), digits = 0) 
  names(PC_perc) <- c(paste0("PC", rep(1:4))) }
  
  ## Saving genetic PCA axes to plot later:
  pcadata <- as.data.frame(PCA$projections)
  pcadata <- pcadata[, 1:4]
  names(pcadata) <- c(paste0("PC", rep(1:4)))
  
  ## Loading assigments:
  assignments <- read.csv(header = TRUE, file = here(paste0("outputs/assignments_Ctenotus_OTUs.csv")))
  assignments <- assignments[assignments$group == group, ]
  
  ## Merge:
  pcadata$ID <- assignments$ID
  pcadata <- merge(pcadata, assignments, by = "ID")
  
  ## Adding labels to make plots:
  sample_info <- read.csv(here("sample_information/ddRAD_sample_information.csv"), header = TRUE)
  sample_info <- arrange(sample_info, by = SAMPLE_ID)
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% pcadata$ID, ]
  pcadata$label <- as.character(gsub(sample_info$UPDATED_SP, pattern = "(^[a-z]{4}).+", replacement = "\\1"))
  pcadata$taxon <- sample_info$UPDATED_SP
  
  ## Using a phylogeny to order samples in structure plots:
  full_tree <- read.nexus(here("RAxML/Ctenotus.nex"))
  sample_order <- as.data.frame((full_tree$tip.label[full_tree$tip.label %in% pcadata$ID]))
  colnames(sample_order) <- "ID" 
  sample_order$plot_order <- 1:nrow(sample_order) ## Order of samples in the phylogeny.
  pcadata <- merge(pcadata, sample_order, by = "ID")
  pcadata <- arrange(pcadata, by = plot_order)
  
  ## Function will return:
  return(list(pcadata = pcadata, PC_perc = PC_perc))
  
} ## End of function.

## PART 3: Function: Setting plot colors ----

## Color palette:
get_palette <- function(K) {
  colors <- stepped2()[seq(1, 18, 1.25)] 
  color_mask <- ceiling(x = seq(from = 1, to = length(colors), length.out = K))
  palette <- colors[color_mask]
  return(palette)
}

## PART 4: Function: Plotting PCA results ----

## Function: Plotting PCs:
plot_PCA <- function(group, miss_ind) {
  
  ## Testing:
  #group <- "atlas_gr" ; miss_ind <- 0.5 
  
  ## Run PCA:
  PCA_results <- run_PCA(group = group, miss_ind = miss_ind, miss_SNP = miss_SNP)
  
  ## Data we'll use:
  pcadata <- PCA_results$pcadata[!is.na(PCA_results$pcadata$cluster), ]
  pcadata <- pcadata[, c("ID", "cluster", "label", c(paste0("PC", rep(1:4)))) ]
  
  ## If needed, multiply PCs by -1 to match axes of previously made figure:
  ## This does not affect the results, only PC axis orientation, which is arbitrary.
  #for (PC in c("PC2")) {
  #  pcadata[[PC]] <- pcadata[[PC]]*-1
  #}
  
  ## Ordering factors to keep order in plot:
  pcadata$cluster <- factor(pcadata$cluster, levels = unique(pcadata$cluster))
  
  ## Best fit-K:
  assignments <- read.csv(header = TRUE, file = here(paste0("outputs/assignments_Ctenotus_OTUs.csv")))
  assignments <- assignments[assignments$group == group, ]
  K <- length(unique(assignments$cluster))
  
  ## Color palette:
  palette <- get_palette(K)
  
  ## Function:
  plot_PCaxes <- function(PCx, PCy) {
    
    ## Testing:
    #PCx <- "PC1" ; PCy <- "PC2"
    
    ## Status:
    print(paste0("Now plotting ", PCy, " vs. ", PCx, "!"))
    
    ## Plotting PCA scatterplot with ggplot2:
    plot <- ggplot(pcadata) +
      
      ## Plotting ellipses for each group:
      ggforce::geom_mark_ellipse(aes(x = .data[[PCx]], y = .data[[PCy]], fill = cluster, color = cluster), 
                                 expand = unit(2, "mm"), size = 0.25, alpha = 0.3) +
      
      ## Configuring point size and shape:
      geom_point(aes(x = .data[[PCx]], y = .data[[PCy]], fill = cluster), size = 2.5, shape = 21, 
                 alpha = 0.80, stroke = 0.15, position = position_jitter(width = 0.9, height = 0.9)) +
      
      ## Specifying ggplot theme to be used:
      theme_light() +
      
      ## Expand plot limits:
      expand_limits(y = c(min(pcadata[[PCy]])-2, max(pcadata[[PCy]])+2),
                    x = c(min(pcadata[[PCx]])-2, max(pcadata[[PCx]])+2)) +
      
      ## Axis labels:
      labs(y = paste0(PCy, " (", PCA_results$PC_perc[[PCy]], "%)"), 
           x = paste0(PCx, " (", PCA_results$PC_perc[[PCx]], "%)")) +
      
      ## Setting color scheme and, if needed, removing legend:
      scale_fill_manual(values = palette) + 
      scale_color_manual(values = palette) +
      guides(fill = "none", color = "none", label = "none") +
      
      ## A few edits:
      theme(panel.border = element_rect(size = 0.75, colour = "gray30"),
            axis.text.x = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, r = 0, b = 0, l = 0)),
            axis.text.y = element_text(size = 12, hjust = 0.5, margin = margin(t = 0, r = 5, b = 0, l = 0)),
            axis.title.x = element_text(size = 14, hjust = 0.5, margin = margin(t = 5, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 14, hjust = 0.5, margin = margin(t = 0, r = 5, b = 0, l = 0)),
            panel.grid = element_blank(), ## Removing background grid.
            axis.ticks = element_line(size = 0.5, color = "gray30"), ## Making axis thicks thicker.
            axis.ticks.length = unit(0.25, "cm")) ## Making axis ticks longer.
   
    ## Function will return:
    return(plot)
    
  } ## End of function.
  
  ## Plot for PC combinations while getting rid of redundant axis labels:
  plot_A <- plot_PCaxes(PCx = "PC1", PCy = "PC2") + theme(axis.title.x = element_blank())
  plot_B <- plot_PCaxes(PCx = "PC3", PCy = "PC2") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  plot_C <- plot_PCaxes(PCx = "PC1", PCy = "PC4") 
  plot_D <- plot_PCaxes(PCx = "PC3", PCy = "PC4") + theme(axis.title.y = element_blank())
  
  # Status:
  print("Now combining plots!")
  
  ## Combining plots using the patchwork package:
  plots_combined <- ((plot_A  | plot_B ) / (plot_C | plot_D )) +
    
  ## Combining the legends of individual plots:
  plot_layout(guides = "collect") + 
    
    ## Setting common title:
    plot_annotation(title = paste0("Ctenotus ", gsub(x = group, pattern = "_gr", replacement = ""), " group"),
                    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "italic", margin = margin(t = 10, r = 0, b = 10, l = 0))))
  
  ## Saving plot (with cowplot):
  ggsave(filename = here(paste0("sNMF_runs/PCA/", group, "_mi", miss_ind, ".pdf")), 
          plot = plots_combined, width = 10, height = 8, dpi = 200, units = "in")
  
  ## Status: 
  print("Combined plot ready!")
  
  ## Remove pca-related files:
  unlink(here("*.pcaProject"))
  unlink(here("sNMF_runs/data/*.lfmm"))
  unlink(here("sNMF_runs/data/*.pca"), recursive = TRUE)
  
  ## Return:
  return(plot)
  
} ## End of function.

## Set parameters:
miss_ind <- 0.5
groups <- c("atlas_gr", "colletti_gr", "inornatus_gr", "leonhardii_gr", "essingtonii_gr", "pantherinus_gr", "schomburgkii_gr", "taeniolatus_gr")

## Run for all groups:
map(groups, plot_PCA, miss_ind = miss_ind)

## End of script.
