####
### R script by Ivan Prates (ivanprates.org).
### University of Michigan, Ann Arbor, MI, USA.
### May 2020.

### The goals of this script are:
### To run PCA analyses based on the unlinked SNP data.
### To make scatter plots based on the PC axes.

## Installing packages we'll need:
#install.packages("cowplot")
#install.packages("here")
#install.packages("LEA")
#install.packages("tidyverse")

## Loading packages:
library(cowplot)
library(here)
library(LEA)
library(patchwork)
library(tidyverse)

## PART 1: Performing PCA analyses: ----

## Function: Running PCA:
run_PCA <- function(K){

  ## Selecting qmatrix for focal K value:
  qmatrix <- sNMF_results$qmatrix_list[[K-minK+1]]
  
  ## Order qmatrix by sample ID:
  qmatrix <- arrange(qmatrix, SAMPLE_ID)
  
  ## Location of geno file:
  geno <- paste0(here(), "/LEA/", ipyrad_run, "/mac", mac, "/sNMF_analyses/mid", mid, "_ms", ms, "_a", a, ".geno")
  
  ## Run PCA analysis:
  PCA <- pca(geno, center = TRUE, scale = FALSE)
  
  ## Displaying information on the PCA analysis:
  #summary(PCA)
  
  ## Plotting eigenvalues:
  #plot(PCA, lwd = 5, col = "blue", cex = 0.7, xlab = "Factors", ylab = "Eigenvalues")
  
  ## Plotting standard deviations:
  #plot(PCA$sdev)
  
  ## Performing Tracy-Widom tests for all eigenvalues:
  ## Create file: genotypes.tracyWidom - tracy-widom test information, in the directory genotypes.pca/
  #tw <- tracy.widom(PCA)
  
  ## Plotting the percentage of variance explained by each component:
  #plot(tw$percentage)
  
  ## Plot p-values for the Tracy-Widom test:
  #plot(tw$pvalues)
  
  ## Saving genetic PCA axes to plot later:
  pcadata <- PCA$projections # This is what we'll plot (the pc axes).
  pcadata <- as.data.frame(pcadata)
  
  ## Storing the proportion of variation represented by selected PC axes. 
  # Creating an empty list that has the length of the number of PCs we want to keep.
  perc_PCs <- vector("list", 5)
  
  ## Looping over PC axes:
  for (p in c(1:5)) { # We'll use the five first PC axes.
    temp <- round((summary(PCA)[2, p]*100), digits = 1)
    perc_PCs[[p]] <- temp ## Creating temporary object.
    rm(temp) ## Removing temporary variable.
    rm(p) ## Removing temporary variable.
    } ## Close loop. 
  
  ## Grouping PCA results based on SNMF clusters:
  pcadata$cluster = qmatrix$cluster_assign
  
  ## Changing PC names:
  colnames(pcadata)[1:length(qmatrix$ID)] = c(paste0(rep("PC", length(qmatrix$ID)), 1:length(qmatrix$ID)))

  ## Function will return:
  return(list(pcadata = pcadata,
  perc_PCs = perc_PCs))
  
  } ## End of function.

## PART 2: Plotting PCA results: ----
  
## Function: Plotting PCs:
plot_PCA <- function(K, PCx, PCy){
  
  ## Data we'll use:
  pcadata <- PCA_results$pcadata
  perc_PCs <- PCA_results$perc_PCs
  
  ## For each K, selecting color palette for plots:
  cluster_palette <- palette_list[[K-minK+1]]
  
  ## For each K, selecting list of putative taxa for plot legends:
  putative_taxa <- putative_taxa_list[[K-minK+1]]
  
  ## Defining arguments to name the axes on plots:     
  PCx_plot = pcadata[, PCx]
  label_x = paste0("PC", PCx, " (", perc_PCs[[PCx]], "%)")
  PCy_plot = pcadata[, PCy]
  label_y = paste0("PC", PCy, " (", perc_PCs[[PCy]], "%)")

  ## Plotting PCA scatterplot with ggplot2:
  PCA_plot <- ggplot(pcadata, aes(x = PCx_plot, y = PCy_plot, fill = cluster)) +
  
    ## Adding individual labels (for verification purposes only; makes the plot very busy):
    #geom_text_repel(aes(x = PCx_plot, y = PCy_plot, label = qmatrix$SAMPLE_ID), size = 5) + 
  
    ## Plotting ellipses around the points corresponding to each dewlap:
    #stat_ellipse(aes(x = PCx_plot, y = PCy_plot, color = cluster), level = 0.99, size = 1) +
    
    ## Configuring point size and shape:
    geom_point(size = 6, shape = 21) +
  
    ## Specifying ggplot theme to be used:
    theme_light() +
  
    ## Adding labels to the graph:
    labs(y = label_y, x = label_x) + #, title = "Ctenotus inornatus group: PCA on 1537 unlinked SNPs") +
        
    ## Setting color scheme and, if needed, adding legend:
    scale_fill_manual(name = "Genetic clusters:", values = cluster_palette, labels = putative_taxa) + 
    
    ## In case elipses are used, set colors:
    #scale_color_manual(values = cluster_palette, guide = "none") + 
    
    ## Changing text:
    theme(
    
          ## Changing thickness of lines around plot:
          panel.border = element_rect(size = 3, colour = "gray30"),
              
          ## Changing title:  
          plot.title = element_text(size = 40, hjust = 0.5, margin = margin(t = 10, r = 0, b = 20, l = 0)),
              
          ## Changing axes' text:
          axis.text.x = element_text(size = 32, hjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 32, hjust = 0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
              
          ## Changing axes' labels:
          axis.title.x = element_text(size = 36, hjust = 0.5, margin = margin(t = 20, r = 0, b = 10, l = 0)),
          axis.title.y = element_text(size = 36, hjust = 0.5, margin = margin(t = 0, r = 20, b = 0, l = 10)),
              
          ## Changing legend title and text:
          legend.text = element_text(size = 30, margin = margin(t = 0, r = 30, b = 0, l = 0)), 
          legend.title = element_text(size = 32, margin = margin(t = 0, r = 30, b = 0, l = 0)),
          
          ## Changing vertical space between elements of the legend:
          legend.key.size = unit(2.5, 'lines'),
          
          ## A few more changes:
          panel.grid = element_blank(), ## Removing background grid.
          axis.ticks = element_line(size = 2, color = "gray30"), ## Making axis thicks thicker.
          axis.ticks.length = unit(.25, "cm") ## Making axis ticks longer.
          ) ## Close "theme" settings.

    ## Saving plot (with cowplot):
    save_plot(filename = paste0(here(), "/LEA/", ipyrad_run, "/mac", mac, "/PCA/K", K, "/K", K, "_PC", PCy, "-PC", PCx, ".jpeg"), plot = PCA_plot, base_width = 25, base_height = 20)
        
    ## Function will return:
    return(PCA_plot = PCA_plot)
    
  } ## End of function.

## PART 3: Running the functions we created:
  
## If needed, setting some parameters again:
#ipyrad_run <- "inornatus_gr_c90_ni281_mi070"
#mac <- 2
#mid <- 0.8
#ms <- 0.35
#a <- 500
#minK <- 13
#maxK <- 15

## If needed, remove existing folder to start again from scratch:
#unlink(paste0(here(), "/LEA/", ipyrad_run, "/mac", mac, "/PCA"), recursive = TRUE)

## Creating a directory to place PCA results:
dir.create(paste0(here(), "/LEA/", ipyrad_run,  "/mac", mac, "/PCA"))

## Setting focal K:
for (K in minK:maxK){ 

  ## Creating a directory to place PCA results:
  dir.create(paste0(here(), "/LEA/", ipyrad_run,  "/mac", mac, "/PCA/K", K))
  
  ## Run PCA analysis:
  PCA_results <- run_PCA(K)

  ## Plot PCA results :
  PCA_plot_A <- plot_PCA(K, PCx = 1, PCy = 2)
  PCA_plot_B <- plot_PCA(K, PCx = 3, PCy = 2)
  PCA_plot_C <- plot_PCA(K, PCx = 1, PCy = 4)
  PCA_plot_D <- plot_PCA(K, PCx = 3, PCy = 4)
  
  ## Combining plots using the patchwork package:
  ## We'll also get rid of redundant axis labels:
  PCA_plot_combined <- (PCA_plot_A + theme(axis.title.x = element_blank()) |
                        PCA_plot_B + theme(axis.title.x = element_blank(),
                                           axis.title.y = element_blank()) )/
                       (PCA_plot_C |
                        PCA_plot_D + theme(axis.title.y = element_blank()) )+
                          
  ## Combining the legends of individual plots:
  plot_layout(guides = "collect") + 
  
  ## Setting common title:
  plot_annotation(title = expression(paste(italic("Ctenotus inornatus"), " species group: PCA on 1537 unlinked SNPs")),
                  theme = theme(plot.title = element_text(
                    size = 42, hjust = 0.5, ## Title size and justification.
                    margin = margin(t = 30, r = 0, b = 30, l = 0)))) ## Margins around title.
                                
  ## Saving plot (with cowplot):
  save_plot(filename = paste0(here(), "/LEA/", ipyrad_run, "/mac", mac, "/PCA/K", K, "/K", K, "_PCA_plots_combined.jpeg"), 
            plot = PCA_plot_combined, base_width = 25, base_height = 20)

  } ## Close loop (K values).

## PART 4: Wrapping up: ----

## Saving our progress as an R workspace image:
save.image(file = paste0(here(), "/LEA/", ipyrad_run, "/mac" , mac, "/sNMF_workspaces/mid", mid, "_ms", ms, "_a", a, ".RData"))

## End of script.
