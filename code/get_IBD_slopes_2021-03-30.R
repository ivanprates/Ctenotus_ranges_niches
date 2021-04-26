####################
### The goals of this script are:
### To estimate IBD for Ctenotus OTUs and correponding taxa.
### To plot the slope of IBD vs. geographic range size and climatic niche breadth.
### To test for significant relationships between these variables.
### Checked for functionality on April 22nd 2021.

## Useful resources:
## https://popgen.nescent.org/2015-05-18-Dist-SNP.html
## https://www.rdocumentation.org/packages/hierfstat/versions/0.5-7/topics/genet.dist
## https://www.rdocumentation.org/packages/poppr/versions/2.8.5/topics/nei.dist
## https://pubmed.ncbi.nlm.nih.gov/20067366/ ## Study comparing metrics of genetic distance.

## PART 1: Getting ready ----

## Loading packages:
library(ape)
library(BEDASSLE)
library(caper)
library(fossil)
library(LEA)
library(patchwork)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## Creating directories to save results:
dir.create(path = here("IBD"))
dir.create(path = here("IBD/data"))
dir.create(path = here("IBD/dist_matrices"))
dir.create(path = here("IBD/MMRR_results"))
dir.create(path = here("IBD/plots"))

## Loading environmental and OTU information:
environ_data <- read.csv(file = here("climatic_data_Ctenotus_OTUs/extracted_climatic_data.csv"), header = TRUE)
assignments <- read.csv(file = here("outputs/assignments_Ctenotus_OTUs.csv"), header = TRUE)
environ_data <- merge(environ_data, assignments, by = "ID")

## PART 2: Function: Preparing genetic and environmental data ----

## Function:
get_data <- function(group) {

  ## Testing:
  #group <- "inornatus_gr"
  
  ## Status:
  print(paste0("Now preparing data for ", group, "!"))
    
  ## Keeping only focal group:
  environ_data <- environ_data[environ_data$group == group, ] 
  
  ## Loading genetic and sample information:
  gendata <- as.data.frame(read.geno(input.file = here(paste0("sNMF_runs/data/", group, "_mi0.5.geno"))))
  ID <- read.table(file = here(paste0("sNMF_runs/samples/samples_", group, "_mi0.5.txt")), header = FALSE)
  gendata$ID <- ID$V1
  
  ## What samples are present in both the genetic and environmental dataset?
  overlapping <- intersect(environ_data$ID, gendata$ID)
  
  ## Keep only overlapping samples in both datasets. We'll also adjust factor levels:
  gendata <- gendata[gendata$ID %in% overlapping, ]
  gendata$ID <- factor(gendata$ID, levels = gendata$ID)
  environ_data <- environ_data[environ_data$ID %in% overlapping, ]
  environ_data$ID <- factor(environ_data$ID, levels = environ_data$ID)

  ## Check: Is the order of samples the same in the genetic and environmental data?
  table(environ_data$ID == gendata$ID)
  
  ## Function will return:
  return(list(gendata = gendata, environ_data = environ_data))
  
} ## End of function.

## PART 3: Function: Estimating IBD and IBE ----

## Function:
get_dist <- function(group, group_var, min_samples){
  
  ## Testing:
  #group <- "colletti_gr" ; group_var <- "cluster" ; min_samples <- 2
  
  ## First get environmental and genetic data:
  datasets <- get_data(group = group)
  
  ## Select environmental data:
  environ_data <- datasets$environ_data
  
  ## Get rid of unassigned samples:
  environ_data <- environ_data[(!is.na(environ_data[[group_var]])), ]
  
  ## Get rid of OTUs with less than a minimum of samples:
  keep <- names(table(environ_data[[group_var]])[table(environ_data[[group_var]]) >= min_samples])
  environ_data <- environ_data[environ_data[[group_var]] %in% keep, ]
  
  ## Function:
  calc_dist <- function(OTU) {
  
    ## Testing: 
    #OTU = "inornatus_gr_cluster_1"
 
    ## Print status:
    print(paste0("Estimating Fst for ", OTU, " in the ", group_var, " scheme!"))
    
    ## Data for focal OTU:
    environ_df <- environ_data[environ_data[[group_var]] == OTU, ]
    
    ## We're gonna consider individuals separately:
    environ_df$site <- environ_df$ID
    
    ## Estimating a geographic distance matrix:
    geo_dist <- earth.dist(environ_df[, c("Longitude", "Latitude")], dist = FALSE)
    colnames(geo_dist) = environ_df$site
    rownames(geo_dist) = environ_df$site
    diag(geo_dist) <- NA ## Replacing diagonals with NA.
    
    ## Estimating an environmental distance matrix:
    environ_var <- environ_df[, c("Annual_mean_temperature", "Temperature_seasonality", "Annual_precipitation", "Precipitation_seasonality")]
    eco_dist <- as.matrix(dist(environ_var, method = "euclidean"))
    colnames(eco_dist) = environ_df$site
    rownames(eco_dist) = environ_df$site
    diag(eco_dist) <- NA ## Replacing diagonals with NA.
    
    ## Genetic data:
    gendata <- datasets$gendata
    
    ## Keep only OTU samples:
    gendata <- gendata[gendata$ID %in% environ_df$ID, ]
    
    ## Check: Is the order of samples the same in the genetic and environmental data?
    table(environ_df$ID == gendata$ID)
    
    ## Add rownames:
    rownames(gendata) <- gendata$ID ## Use IDs as pop names.
    
    ## Remove ID column:
    SNP_df <- gendata[, names(gendata) != "ID"]
    
    ## Estimating Fst using BEDASSLE:
    ## Let's first prepare the data format.
    ## Convert to matrix with numeric column values:
    allele_counts <- apply(SNP_df, 2, as.numeric) 
      
    ## Add row names again:
    row.names(allele_counts) <- row.names(SNP_df)
      
    ## Create a sample size object and reserve:
    sample_sizes <- allele_counts
      
    ## Change 9 (missing data) for zeros in the allele counts matrix:
    allele_counts[allele_counts == 9] <- 0     
      
    ## Removing invariable sites:
    SNPs_to_keep <- apply(allele_counts, 2, function(x) length(unique(x)) > 1)
    allele_counts <- allele_counts[, SNPs_to_keep]
      
    ## Now, in the sample size object, keep only SNPs present in the final allele count matrix:
    sample_sizes <- sample_sizes[, colnames(allele_counts)]
      
    ## Since organism is diploid, replace with total number of alleles per sample:
    sample_sizes[sample_sizes == 1 ] <- 2
    sample_sizes[sample_sizes == 0 ] <- 2
      
    ## Change NAs for zeros in the sample size matrix:
    sample_sizes[sample_sizes == 9] <- 0
      
    ## Estimate Fst using BEDASSLE:
    Fst_dist <- calculate.all.pairwise.Fst(allele.counts = allele_counts, sample.sizes = sample_sizes)
      
    ## Linearizing Fst:
    Fst_dist[Fst_dist == 1] <- 0.99 ## If needed, transforming 1 to 0.99 or transformation will result in Inf.
    gen_dist <- Fst_dist/(1-Fst_dist)
    
    ## Convert to matrix:
    gen_dist <- as.matrix(gen_dist)
    
    # ## Check: Is matrix symetric? Must be "TRUE":
    isSymmetric(gen_dist)
    
    ## Change column and row names, set diagonal as NA:
    colnames(gen_dist) <- row.names(SNP_df)
    rownames(gen_dist) <- row.names(SNP_df)
    diag(gen_dist) <- NA
    
    ## Check: is data paired across matrices? Must all be "TRUE".
    table(colnames(gen_dist) == colnames(geo_dist))
    table(colnames(gen_dist) == colnames(eco_dist))
    
    ## Save matrices:
    write.csv(x = gen_dist, file = here(paste0("IBD/dist_matrices/gen_dist_", group_var, "_", OTU, ".csv")), quote = FALSE)
    write.csv(x = eco_dist, file = here(paste0("IBD/dist_matrices/eco_dist_", group_var, "_", OTU, ".csv")), quote = FALSE)
    write.csv(x = geo_dist, file = here(paste0("IBD/dist_matrices/geo_dist_", group_var, "_", OTU, ".csv")), quote = FALSE)
    
    ## Turn dist matrices into vectors:
    eco_vector <- eco_dist[lower.tri(eco_dist, diag = FALSE)]
    gen_vector <- gen_dist[lower.tri(gen_dist, diag = FALSE)]
    geo_vector <- geo_dist[lower.tri(geo_dist, diag = FALSE)]
    
    ## Return:
    return(data.frame(OTU = OTU, Climatic_distance = eco_vector, Genetic_distance = gen_vector, Geographic_distance = geo_vector))
    
    ## Print status:
    print(paste0("Done processing ", OTU, "!"))
    
  } ## End of function.
  
  ## Apply:
  all_dist_df <- map_df(unique(environ_data[[group_var]]), calc_dist)
  
  ## Saving:  
  write.csv(x = all_dist_df, row.names = FALSE, file = here(paste0("IBD/data/distances_", group, "_", group_var, ".csv")), quote = FALSE)
  
  ## Return:
  return(all_dist_df)
  
} ## End of function.

## PART 4: Function: Combine IBD for different groups ----

## Function:
get_dist_all <- function(groups, group_var) {
  
  ## Testing:
  #group_var <- "cluster"
  #groups <- c("atlas_gr", "colletti_gr", "inornatus_gr", "leonhardii_gr", "pantherinus_gr", "schomburgkii_gr", "taeniolatus_gr")
              
  ## Function: Get data for each group:
  get_data <- function(group) {
    print(paste0("Now running group ", group, "!"))
    dist_data <- read.csv(header = TRUE, file = here(paste0("IBD/data/distances_", group, "_", group_var, ".csv")))
    dist_data$OTU <- as.character(dist_data$OTU)
    dist_data$group <- group
    return(dist_data)
  } ## End of function.

  ## Apply:
  dist_df <- map_df(groups, get_data)
  
  ## Save:
  write.csv(dist_df, row.names = FALSE, quote = FALSE, file = here(paste0("IBD/data/distances_all_", group_var, ".csv")))
  
  ## Return:
  return(dist_df)
  
} ## End of function.

## PART 5: Function: Extracting slope of IBD and IBE ----

## Function:
get_slope <- function(group_var) {
  
  ## Testing:
  #group_var <- "cluster"
  
  ## Import data from a previous run:
  dist_df <- read.csv(file = here(paste0("IBD/data/distances_all_", group_var, ".csv")), header = TRUE)
  
  ## Function:
  calc_slope <- function(OTU) {
    
    ## Testing:
    #OTU <- "atlas_gr_cluster_7"
    
    ## Data for OTU:
    OTU_df <- dist_df[dist_df$OTU == OTU, ]
    
    ## Perform linear model to get slope of IBD.
    ## Note: lm cannot be used to test significance for distance matrices. For this, use MMRR.
    lm <- lm(Genetic_distance ~ Geographic_distance, data = OTU_df)
      
    ## Extracting the slope and group and saving:
    slope_OTU <- data.frame(OTU = OTU, slope = lm$coefficients[2])
    
    ## Return:
    return(slope_OTU)
    
  } ## End of function.
  
  ## Apply:
  slope_df <- map_df(unique(dist_df$OTU), calc_slope)
  
  ## Save:
  if (group_var == "cluster") { group_label <- "OTUs" }
  if (group_var == "taxon") { group_label <- "taxa" }
  write.csv(slope_df, file = here(paste0("outputs/IBD_slopes_Ctenotus_", group_label, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Function will return:
  return(slope_df)
  
} ## End of function.

## PART 6: Function: Plotting IBD ----

## Function:
plot_IBD <- function(group, group_var, min_comparisons) {
  
  ## Testing:
  #group <- "inornatus_gr" ; group_var <- "cluster" ; min_comparisons <- 2
  
  ## Import distances from a previous run:
  dist_df <- read.csv(file = here(paste0("IBD/data/distances_", group, "_", group_var, ".csv")), header = TRUE, stringsAsFactors = FALSE)
  
  ## Import data on niche breadths and range sizes:
  if (group_var == "cluster") { group_label <- "OTUs" }
  if (group_var == "taxon") { group_label <- "taxa" }
  nb_rs <- read.csv(file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_", group_label, ".csv")), header = TRUE)
  
  ## Merge with the distance data:
  dist_df <- merge(dist_df, nb_rs, by.x = "OTU", by.y = "group")
  
  ## What groups have a minimum number of comparisons?
  keep_groups <- names(table(dist_df$OTU) >= min_comparisons)[table(dist_df$OTU) >= min_comparisons]
  
  ## Keep only groups with enough samples:
  dist_df <- dist_df[dist_df$OTU %in% keep_groups, ]
  
  ## Setting facet labeller:
  label_df <- read.csv(file = here(paste0("sNMF_runs/samples/labels_", group, "_mi0.5_delimited_OTUs.csv")), header = TRUE, stringsAsFactors = FALSE)
  label_df$OTU <- paste0(group, "_", label_df$cluster)
  
  ## Labels:
  if (group_var == "cluster") {
    dist_df$label <- factor(dist_df$OTU, levels = label_df$OTU, labels = label_df$label) 
    group_label <- "OTUs" }
  if (group_var == "taxon") { 
    dist_df$label <- dist_df$OTU 
    group_label <- "taxa" }
  
  ## Plot:
  plot <- ggplot(data = dist_df) + 
      
      ## Add points:
      geom_point(aes(x = Geographic_distance, y = Genetic_distance), color = "gray50", size = 2, shape = 19) +
      
      ## Adding regression line:
      geom_smooth(aes(x = Geographic_distance, y = Genetic_distance), method = "lm", color = "blue", size = 0.75, se = FALSE) +
      
      ## Setting axis labels:
      ylab("Genetic distances (Fst/1-Fst)") +
      xlab("Geographic distances (Km)") +
      labs(title = paste0("Ctenotus ", gsub(x = group, pattern = "_gr", replacement = ""), " group")) +  
    
      ## Facet by group:
      facet_wrap(~label, ncol = 4) +

      ## Changing overall theme:
      theme_bw() +
      
      ## Also setting values on axes:
      scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      
      ## Other params:
      theme(
        plot.margin = margin(t = 10, b = 10, l = 10, r = 20),
        strip.text = element_text(size = 18, face = "italic"),
        axis.title.x = element_text(size = 18, margin = margin(t = 10, b = 0, l = 0, r = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, b = 0, l = 0, r = 10)),
        axis.text = element_text(size = 14),
        #plot.title = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5, face = "italic", margin = margin(t = 10, r = 0, b = 10, l = 0)),
        plot.subtitle = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 1, colour = "gray20"), ## Changing box around plot.
        axis.ticks = element_line(size = 0.75, color = "gray20"), ## Making axis thicks thicker.
        panel.grid.minor = element_blank()) ## Removing background grid.
    
    ## Plot configuration and size:
    bh <- 1+(ceiling(length(unique(dist_df$OTU))/4)*2.5)
    if (group == "colletti_gr") { 
      bw <- 6
    } else if (group == "taeniolatus_gr") { 
      bw <- 4 
    } else if (group == "pantherinus_gr" & group_var == "taxon") { 
      bw <- 4 
    } else { 
      bw <- 12 
    }
    
    ## Save plot in pdf format:
    ggsave(plot = plot, height = bh, width = bw, device = "pdf", units = "in", dpi = 500,
              filename = here(paste0("IBD/plots/", group, "_", group_label, "_IBD.pdf")))
    
} ## End of function.

## PART 7: Function: Plot slope of IBD vs. niche breadth and range size ----

## Function:
plot_slope <- function() {
  
  ## Testing:
  #group_var <- "cluster" ; xvar <- "log_rs" ; facet_var <- "no"
  
  ## Function:
  plot_IBD_correlate <- function(xvar, group_var) {
  
    ## Status:
    print(paste0("Now plotting IBD slopes for the ", group_var, " scheme!"))
    
    ## Load IBD slopes:
    if (group_var == "cluster") { group_label <- "OTUs" }
    if (group_var == "taxon") { group_label <- "taxa" }
    slope_df <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_", group_label, ".csv")))
    
    ## Load range sizes and climatic niches and merge:
    rs_nb_df <- read.csv(header = TRUE, file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_", group_label, ".csv")))
    slope_df <- merge(slope_df, rs_nb_df, by.x = "OTU", by.y = "group")
    
    ## Remove NAs:
    slope_df <- slope_df[!is.na(slope_df$range_size), ]
    
    ## Define plot labels:
    if (group_var == "taxon") { ylabel <- "Slope of IBD (taxa)" }
    if (group_var == "cluster") { ylabel <- "Slope of IBD (OTUs)" }
    if (xvar == "log_rs") { xlabel <- "Range size"
                            xmin <- 2
                            xmax <- 7 } 
    if (xvar == "log_nb") { xlabel <- "Climatic niche breadth"
                            xmin <- -4
                            xmax <- 1 } 
    
    ## Plot:
    plot <- ggplot(data = slope_df) +
          
      ## Scatter plot:    
      geom_point(aes(x = .data[[xvar]], y = slope), size = 2.5, color = "gray50") +
          
      ## Adding regression line:
      geom_smooth(aes(x = .data[[xvar]], y = slope), method = "lm", color = "blue", 
                  size = 0.75, se = FALSE) +
          
      ## Setting title:
      ylab(label = ylabel) +
      xlab(label = xlabel) +
          
      ## Changing overall theme:
      theme_bw() +
          
      ## Also setting values on axes:
      scale_x_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(xmin, xmax)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
      
      ## Other params:
      theme(
        plot.margin = margin(t = 10, b = 10, l = 10, r = 10),
        axis.title.x = element_text(size = 18, margin = margin(t = 10, b = 0, l = 0, r = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, b = 0, l = 0, r = 5)),
        axis.text = element_text(size = 16),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        panel.border = element_rect(size = 1, colour = "gray20"), ## Changing box around plot.
        axis.ticks = element_line(size = 0.75, color = "gray20"), ## Making axis thicks thicker.
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) ## Removing background grid.
        
    
    ## Return:
    return(plot)
    
  } ## End of function.
  
  ## Plot:
  plot_crs <- plot_IBD_correlate(xvar = "log_rs", group_var = "cluster")
  plot_cnb <- plot_IBD_correlate(xvar = "log_nb", group_var = "cluster")
  plot_trs <- plot_IBD_correlate(xvar = "log_rs", group_var = "taxon")
  plot_tnb <- plot_IBD_correlate(xvar = "log_nb", group_var = "taxon")
  
  ## Combine and save for OTUs:
  plot_c <- ( plot_crs ) + ( plot_cnb + theme(axis.title.y = element_blank(), axis.text.y = element_blank()) ) +
    plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(size = 24),
            plot.tag.position = c(0.95, 0.95))
  ggsave(plot = plot_c, height = 5, width = 11, device = "jpeg", units = "in", dpi = 500,
            filename = here(paste0("plots/Fig6.jpeg")))
  
  ## Combine and save for taxa:            
  plot_t <- ( plot_trs | plot_tnb + theme(axis.title.y = element_blank(), axis.text.y = element_blank()) ) +
     plot_annotation(tag_levels = "A") &
       theme(plot.tag = element_text(size = 24),
             plot.tag.position = c(0.95, 0.95))
   ggsave(plot = plot_t, height = 5, width = 11, device = "jpeg", units = "in", dpi = 500, 
          filename = here(paste0("IBD/plots/slope_of_IBD_taxa.jpeg")))
  
  ## Return:
  return(plot_c)
  
} ## End of function.
    
## PART 8: Run it all ----

## Setting a few parameters:
schemes <- c("cluster", "taxon")
groups <- c("atlas_gr", "colletti_gr", "essingtonii_gr", "inornatus_gr", "leonhardii_gr", "pantherinus_gr", "schomburgkii_gr", "taeniolatus_gr")
min_samples <- 2
min_comparisons <- 2
  
## Get IBD:
for (group_var in schemes) { 
  for (group in groups) { 
    get_dist(group, group_var, min_samples) 
  }
}

## Get slope of IBD:
for (group_var in schemes) { 
  get_dist_all(groups = groups, group_var) 
  get_slope(group_var) 
}

## Plot:
for (group_var in schemes) {
  for (group in groups) { 
    plot_IBD(group, group_var, min_comparisons)
  }
}

## Fig 6:
plot_slope()

## PART 9: Function: Getting IBD significance for OTUs ----

## Function:
get_MMRR_p <- function(group_var) {

  ## Testing:
  #group_var <- "cluster"
  
  ## Load data:
  IBD_df <- read.csv(header = TRUE, file = here(paste0("IBD/data/distances_all_", group_var, ".csv")))
  
  ## Samples with two or more comparisons:
  keep <- names(table(IBD_df$OTU))[table(IBD_df$OTU) >= 2]
  
  ## Filter:
  IBD_df <- IBD_df[IBD_df$OTU %in% keep, ]
  
  ## Function:
  test_MMRR <- function(OTU) {
  
    ## Testing:
    #OTU <- "atlas_gr_cluster_6"
    
    ## Status:
    print(paste0("Now running MMRR for ", OTU, "!"))
  
    ## Which group is this OTU from?
    OTU_group <- unique(IBD_df$group[IBD_df$OTU == OTU])
    
    ## Read gen dist data for group:
    gen_dist <- as.matrix(read.csv(header = TRUE, row.names = 1, file = here(paste0("IBD/dist_matrices/gen_dist_", group_var, "_", OTU, ".csv"))))
    
    ## Read geo dist data for group:
    geo_dist <- as.matrix(read.csv(header = TRUE, row.names = 1, file = here(paste0("IBD/dist_matrices/geo_dist_", group_var, "_", OTU, ".csv"))))
    
    ## Fitting an MMRR model:
    MMRR_test <- MMRR(Y = gen_dist, X = list(geo_dist = geo_dist), nperm = 1000)
    
    ## Results:
    MMRR_df <- data.frame(OTU = OTU,
                          group = OTU_group,
                          n_samples = nrow(gen_dist),
                          p_value = MMRR_test$tpvalue[2],
                          r_squared = MMRR_test$r.squared,
                          slope = MMRR_test$coefficients[2], 
                          row.names = "")
                            
    ## Return
    return(MMRR_df)
    
  } ## End of function.

  ## Run:
  MMRR_df_all <- purrr::map_df(unique(IBD_df$OTU), test_MMRR)
  
  ## Keeping only OTUs with range size and climatic niche estimates:
  if (group_var == "cluster") { group_label <- "OTUs" }
  if (group_var == "taxon") { group_label <- "taxa" }
  slope_df <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_", group_label, ".csv")))
  MMRR_df_all <- MMRR_df_all[MMRR_df_all$OTU %in% slope_df$OTU, ]
  
  ## Save:
  write.csv(MMRR_df_all, file = here(paste0("IBD/MMRR_results/MMRR_", group_label, ".csv")))

  ## Number of OTUs with p <= 0.05 and positive IBD slope):
  n_sig <- length(MMRR_df_all$OTU[MMRR_df_all$p_value <= 0.05])
  n_neg <- length(MMRR_df_all$OTU[MMRR_df_all$slope > 0])
  min(MMRR_df_all$slope)
  max(MMRR_df_all$slope)
  
  ## Return:
  return(list(n_sig = n_sig, n_neg = n_neg))
  
} ## End of function.

## Run:
MMRR_res <- get_MMRR_p(group_var = "cluster")
MMRR_res$n_sig
MMRR_res$n_neg

## PART 10: Testing significance ----

## Load data:
slope_df <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_OTUs.csv")))
rs_nb_df <- read.csv(header = TRUE, file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_OTUs.csv")))
slope_df <- merge(slope_df, rs_nb_df, by.x = "OTU", by.y = "group")

## How many samples?
length(unique(slope_df$OTU))

## Summary of data:
summary(slope_df)

## Histograms:
hist(slope_df$log_nb)
hist(slope_df$log_rs)
hist(slope_df$slope)

## Scatter plot:
plot(log_rs ~ log_nb, data = slope_df)
plot(slope ~ log_rs, data = slope_df)
plot(slope ~ log_nb, data = slope_df)

## Linear regression:
lm_res <- lm(log_rs ~ log_nb, data = slope_df) ; summary(lm_res)
lm_res <- lm(slope ~ log_rs, data = slope_df) ; summary(lm_res)
lm_res <- lm(slope ~ log_nb, data = slope_df) ; summary(lm_res)

## PART 11: Phylogenetic regressions ----

## Load data:
slope_df <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_OTUs.csv")))
assignments <- read.csv(file = here("outputs/assignments_Ctenotus_OTUs.csv"), header = TRUE)
assignments <- assignments[assignments$cluster %in% slope_df$OTU, ]
rs_nb_df <- read.csv(header = TRUE, file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_OTUs.csv")))
slope_df <- merge(slope_df, rs_nb_df, by.x = "OTU", by.y = "group")

## Making tree ultrametric:
tree <- read.nexus(here("RAxML/Ctenotus.nex"))
tree <- chronopl(phy = tree, lambda = 1)

## Keeping one random tip per OTU:
keep_samples <- assignments %>% group_by(cluster) %>% sample_n(size = 1)
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep_samples$ID])

## Replace tree tips:
keep_samples$ID <- factor(keep_samples$ID, levels = tree$tip.label)
keep_samples <- arrange(keep_samples, by = ID)
tree$tip.label <- as.character(keep_samples$cluster)

## Prepare data:
row.names(slope_df) <- slope_df$OTU
geiger::name.check(tree, slope_df) ## Data match?

## Run PGLS using caper:
compdata <- comparative.data(tree, slope_df, names.col = "OTU", vcv.dim = 2, warn.dropped = TRUE)

## For range size:
pglsModel <- pgls(slope ~ log_rs, data = compdata)
summary(pglsModel)

## Plot:
plot(slope_df$slope ~ slope_df$log_rs)
abline(lm(slope_df$slope ~ slope_df$log_rs), col = "blue")
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], col = "red")

## For niche breadth:
pglsModel <- pgls(slope ~ log_nb, data = compdata)
summary(pglsModel)

## Plot:
plot(slope_df$slope ~ slope_df$log_nb)
abline(lm(slope_df$slope ~ slope_df$log_nb), col = "blue")
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], col = "red")

## Range size vs. niche breadth:
pglsModel <- pgls(log_rs ~ log_nb, data = compdata) ## Lambda = 1.
summary(pglsModel)

## Plot:
plot(slope_df$log_rs ~ slope_df$log_nb)
abline(lm(slope_df$log_rs ~ slope_df$log_nb), col = "blue")
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2], col = "red")

## PART 12: Checking the effect of biome on IBD slope ----

## Load data:
slope_df <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_OTUs.csv")))
rs_nb_df <- read.csv(header = TRUE, file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_OTUs.csv")))
slope_df <- merge(slope_df, rs_nb_df, by.x = "OTU", by.y = "group")
biome_df <- read.csv(header = TRUE, file = here(paste0("outputs/biomes_Ctenotus_OTUs.csv")))     
biome_df <- merge(slope_df, biome_df, by = "OTU")

## Keep only most frequent biomes?
#biome_df <- biome_df[biome_df$biome %in% c("desert", "trop_grassland"), ]

## ANCOVAs, biome as factor:
ancova1 <- aov(slope ~ log_rs * biome, data = biome_df) ; summary(ancova1) 
ancova1 <- aov(slope ~ log_rs + biome, data = biome_df) ; summary(ancova1) ## No range size vs. biome interaction.
ancova2 <- aov(slope ~ log_nb * biome, data = biome_df) ; summary(ancova2) 
ancova2 <- aov(slope ~ log_nb + biome, data = biome_df) ; summary(ancova2) ## No niche breadth vs. biome interaction.

## Prepare column and biome names to plot:
names(biome_df)[names(biome_df) == "log_rs"] <- "Range size"
names(biome_df)[names(biome_df) == "log_nb"] <- "Climatic niche breadth"
names(biome_df)[names(biome_df) == "slope"] <- "Slope of IBD"
biome_df$No_biomes <- factor(biome_df$n_biomes)
biome_df$biome <- gsub(x = biome_df$biome, pattern = "desert", replacement = "Desert")
biome_df$biome <- gsub(x = biome_df$biome, pattern = "mediter_woodland", replacement = "Mediterranean woodlands")
biome_df$biome <- gsub(x = biome_df$biome, pattern = "temp_forest", replacement = "Temperate forests")
biome_df$biome <- gsub(x = biome_df$biome, pattern = "trop_grassland", replacement = "Tropical grasslands")

## Function:
plot_by_biome <- function(xvar, yvar) {
  
  ## Plot:
  plot <- ggplot(data = biome_df) +
    
    ## Points:
    geom_point(aes(x = .data[[xvar]], y = .data[[yvar]], fill = biome, size = No_biomes), color = "black", shape = 21, stroke = 0.5, alpha = 0.7, show.legend = TRUE) +
    
    ## Curves:
    geom_smooth(data = biome_df %>% filter(biome == "Desert"), aes(x = .data[[xvar]], y = .data[[yvar]]), method = "lm", se = FALSE, size = 0.5, color = "red") +
    geom_smooth(data = biome_df %>% filter(biome == "Tropical grasslands"), aes(x = .data[[xvar]], y = .data[[yvar]]), method = "lm", se = FALSE, size = 0.5, color = "blue") +
    
    ## Editing legends:
    guides(fill = guide_legend(title = "Most frequent biome for OTU:")) +
    scale_fill_manual(values = c("red", "green", "yellow", "blue")) +
    scale_size_manual(values = c(2,3,4,5,6)) +
    
    ## Editing legends and axis labels:
    scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    
    ## Other params:
    theme_bw() +
    theme(plot.title = element_blank(),
          axis.title = element_text(size = 14), #margin = margin(t = 5, l = 5, b = 5, r = 5)),
          axis.text = element_text(size = 12), #margin = margin(t = 5, l = 5, b = 5, r = 5)),
          panel.border = element_rect(size = 1, colour = "gray20"), ## Changing thickness of lines around plot.
          axis.ticks = element_line(size = 0.75, color = "gray20"), ## Ticks thickness.
          panel.grid = element_blank(),
          
          ## Legend:
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key = element_rect(fill = "transparent", colour = "transparent"), ## Getting rid of boxes around legend elements.
          legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80")) ## Box around legend.
  
  ## Save:
  #ggsave(plot = plot, filename = here(paste0("IBD/plots/by_biome_", yvar, "_vs_", xvar, ".jpeg")), dpi = 150, width = 8, height = 5)
  
  ## Return:
  return(plot)
  
}

## Plot:
A <- plot_by_biome(xvar = "Range size", yvar = "Slope of IBD")
B <- plot_by_biome(xvar = "Climatic niche breadth", yvar = "Slope of IBD")
C <- plot_by_biome(xvar = "Climatic niche breadth", yvar = "Range size")

## Combine:
ABC <- C + A + B +  plot_layout(guides = 'collect')
#ggsave(plot = ABC, filename = here(paste0("plots/by_biome_ABC.jpeg")), dpi = 300, width = 20, height = 5)

AB <- A + B +  plot_layout(guides = 'collect')
ggsave(plot = AB, filename = here(paste0("plots/FigS4.jpeg")), dpi = 300, width = 14, height = 5)

## End of script.
