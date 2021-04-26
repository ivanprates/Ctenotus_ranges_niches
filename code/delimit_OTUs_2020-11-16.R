####################
### The goals of this script are:
### To run genetic clustering analyses using sNMF.
### To make bar plots based on the ancestry coefficients.
### To make maps based on cluster assignments.
### To extract and organize cluster assignments for downstream analyses.
### Checked for functionality on April 22nd 2021.

## PART 1: Getting ready ----

## Packages:
library(ggtree)
library(LEA)
library(magrittr)
library(pals)
library(patchwork)
library(phytools)
library(tidyverse)
library(rgeos)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

## If needed, clearing working space:
rm(list = ls())

## If needed, remove previous runs:
#unlink("/home/ivan/sNMF_runs/", recursive = TRUE)

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## rhinellax.
runs_path <- "/home/ivan/sNMF_runs/" ## rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## Creating directories:
dir.create(runs_path)
dir.create(here("sNMF_runs"))

## Creating directories to save results:
dir.create(path = here("sNMF_runs/assignments"))
dir.create(path = here("sNMF_runs/data"))
dir.create(path = here("sNMF_runs/miss_data"))
dir.create(path = here("sNMF_runs/plots"))
dir.create(path = here("sNMF_runs/qmatrices"))
dir.create(path = here("sNMF_runs/samples"))

## PART 2: Function: Calculating missing data in the SNP dataset ----

## Function:
calc_miss <- function(group) {

  ## Testing:
  #group = "inornatus_gr"
  
  ## Reading genetic data:  
  ## Missing sites should be coded as 9! Don't forget to replace -1 to 9 if using VCF Tools.
  gendata <- read.table(here("genetic_data_filtered_vcftools/", group, "/usnps.012"), sep = "\t", row.names = 1)
  dim(gendata)
  
  ## Testing:
  #gendata <- gendata[, 1:1000]

  ## Let's now estimate levels of missing data per SNP:
  print("Now estimating missing data for SNPs!")
  
  ## Function:
  calc_miss_SNP <- function(SNP) {
  
    ## Let's now estimate levels of missing data per sample:
    print(paste0("Now processing SNP ", SNP))
    
    ## How many sites have no data?
    miss <- table(gendata[, SNP] == 9)["TRUE"]
    
    ## What proportion is this from the total?
    if (is.na(miss) == TRUE) { 
      SNP_miss <- 0 
    } else { 
      SNP_miss <- round(miss/(dim(gendata)[1]), digits = 2) }
    
    ## Storing:
    SNP_miss <- data.frame(SNP = SNP, miss_prop = SNP_miss)
    SNP_miss$SNP <- as.character(SNP_miss$SNP)
    
    ## Return:
    return(SNP_miss)
    
  } ## End of function.
    
  ## Apply:
  SNP_df <- map_df(names(gendata), calc_miss_SNP)
  
  ## Save missing data information:
  write.csv(SNP_df, file = here(paste0("sNMF_runs/miss_data/SNP_df_", group, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Let's now estimate levels of missing data per sample:
  print("Now estimating missing data for samples!")
  
  ## Import sample IDs:
  sample_IDs <- read.table(here("genetic_data_filtered_vcftools/", group, "/usnps.012.indv"), header = FALSE)
  names(sample_IDs) <- "SAMPLE_ID"  
  gendata$SAMPLE_ID <- sample_IDs$SAMPLE_ID
  
  ## Function:
  calc_miss_ind <- function(sample) {
      
    ## Testing:
    #sample <- "WAMR_166392_Ct_hele"
    
    ## Let's now estimate levels of missing data per sample:
    print(paste0("Now processing sample ", sample))
      
    ## Get SNPs for sample:
    sample_data <- gendata[gendata$SAMPLE_ID == sample, ]

    ## What proportion is this from the total number of loci? 
    miss <- table(sample_data == 9)["TRUE"]
    ind_miss <- round(miss/(dim(gendata)[2]), digits = 2)
    
    ## Storing:
    ind_miss <- data.frame(SAMPLE_ID = sample, miss_prop = ind_miss)
    ind_miss$SAMPLE_ID <- as.character(ind_miss$SAMPLE_ID)
    
    ## Return:
    return(ind_miss)
    
  } ## End of function.
        
  ## Apply:
  ind_df <- map_df(gendata$SAMPLE_ID, calc_miss_ind)
  
  ## Save this information:
  write.csv(ind_df, file = here(paste0("sNMF_runs/miss_data/ind_df_", group, ".csv")), row.names = FALSE, quote = FALSE)
  
} ## End of function.

## PART 3: Function: Filtering SNP data by max missing data ----

## Function:
prep_data <- function(group, miss_ind) {
  
  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5
  
  ## Print status:
  print(paste0("Now preparing genetic data for ", group, "!"))
  
  ## Reading genetic data:  
  ## Missing sites should be coded as 9! Don't forget to replace -1 to 9 if using VCF Tools.
  gendata <- read.table(here("genetic_data_filtered_vcftools/", group, "/usnps.012"), sep = "\t", row.names = 1)
  dim(gendata)
  
  ## Import sample IDs:
  sample_IDs <- read.table(here("genetic_data_filtered_vcftools/", group, "/usnps.012.indv"), header = FALSE)
  names(sample_IDs) <- "SAMPLE_ID"  
  
  ## Add sample IDs:
  gendata$SAMPLE_ID <- sample_IDs$SAMPLE_ID
  
  ## If needed, exclude poorly sampled taxa or missidentified samples:
  exclude <- c(
               ## atlas clade:
               #"NA_ABTC125975_Ct_allo", ## Too few samples for taxon.
	             "WAMR_104349_Ct_angu", ## Too few samples for taxon.
               "WAMR_110581_Ct_pian", "WAMR_119711_Ct_pian", ## Now rhabdotus. Too few samples for taxon.
               "NA_CCM1470_Ct_inor", ## Missidentified. rhabdotus?
               "UMMZ_242650_Ct_serv", "UMMZ_242649_Ct_serv", ## Too few samples for taxon.
               #"SAMR_48698_Ct_gran", ## Missidentified. hanloni?
               #"WAMR_158380_Ct_hanl", "WAMR_157583_Ct_hanl", ## Too few samples for taxon.
               #"WAMR_102472_Ct_hele", ## Missidentified. hanloni?
               #"WAMR_140983_Ct_mary", "WAMR_158428_Ct_mary", ## Too few samples for taxon.
               #"AMR_123129_Ct_zast", "UMMZ_242655_Ct_zast", "UMMZ_242654_Ct_zast", ## Too few samples for taxon.
               
	             ## colletti clade:
	             #"NA_CCM1013_Ct_ehma", "NA_CCM1006_Ct_ehma", ## Too few samples for taxon. Also, what species is this?
	             "NA_CCM6484_Ct_deca", ## Missidentified. 
	             "NA_CCM4867_Ct_stor", ## Different species?
	             "WAMR_151014_Ct_tant", ## Missidentified? striaticeps.   
	             "WAMR_151833_Ct_tant", "WAMR_165940_Ct_tant", ## Too few samples per taxon.
	             "NTMR_22441_Ct_stri", ## Too few samples for taxon. 
	             "NA_K014_Ct_stor", "NA_K015_Ct_stor", ## Different species?
	             #"WAMR_102658_Ct_nasu", ## Too few samples for taxon.
	             #"CUMV_14845_Ct_calu", ## Missidentified?
	             #"WAMR_157327_Ct_rufe", ## Too few samples for taxon.
	             #"WAMR_139415_Ct_nigr", "WAMR_139414_Ct_nigr") ## Too few samples for taxon.
	             
               ## essingtonii clade:
               "NTMR_20364_Ct_arnh", "ANWC_R05116_Ct_brev", ## Too few samples for taxon.
               "NA_CCM6499_Ct_essi", ## Different species?
	             "NA_YMR74_Ct_essi", "NA_YMR77_Ct_essi", ## Different species?

               ## inornatus clade:
               "WAMR_126010_Ct_rima", "WAMR_126015_Ct_rima", ## Too few samples for taxon.
               "NA_CCM3739_Ct_inor", ## Different species? Island sample in NT.
               "SAMAR_55731_Ct_robu", "SAMR_55746_Ct_robu", 
	             "NA_CCM0044_Ct_euta", ## Different species?
               "AMSR_111493_Ct_spal", "AMSR_111494_Ct_spal", ## Different species?
	             "QM_82086_Ct_spal", ## Different species?
               #"SAMAR_55799_Ct_late", "NA_CCM0090_Ct_late", ## Different species?
               #"NTMR_20378_Ct_robu", "NTMR_22166_Ct_robu", ## Different species?
               
               ## leonhardii clade:
               "WAMR_146927_Ct_mime", ## Too few samples for taxon.
               "UMMZ_242653_Ct_tana", ## Too few samples for taxon.
               "NTMR_22387_Ct_mili", ## Too few samples for taxon.
               "SAMR_46967_Ct_sept", "SAMAR_44726_Ct_sept", ## Too few samples for taxon.
               "WAMR_110738_Ct_ruti", ## Too few samples for taxon.
               #"WAMR_139532_Ct_iape", "WAMR_157205_Ct_iape", ## Too few samples for taxon.
               "WAMR_139164_Ct_serv", ## Missidentified. Clusters with leonhardii but messes distribution.
               "CUMV_14379_Ct_uber", "WAMR_145646_Ct_uber", "WAMR_154802_Ct_uber", ## Too few samples for taxon.   
               "SAMR_54025_Ct_pulc", "SAMR_54030_Ct_pulc", ## Different species?

               ## pantherinus clade:
               "NA_PMO179_Ct_pant", ## Different species?
               "WAMR_146582_Ct_rubi", ## Too few samples for taxon.
               
               ## schomburgkii clade:
               "WAMR_166437_Ct_scho", ## Different species? 
               "WAMR_158345_Ct_scho", ## Different species?
	             "CUMV_14858_Ct_calu", ## Different species?
               "QM_87405_Ct_rosa", ## Too few samples for taxon.
               
               ## taeniolatus clade:
               "SAMR_65418_Ct_agre", ## Too few samples per taxon.
               "QM_59674_Ct_ingr", ## Too few samples per taxon.
               "QM_84335_Ct_quin") ## Too few samples per taxon.

  ## Exclude:
  exclude <- as.character(gendata$SAMPLE_ID[gendata$SAMPLE_ID %in% exclude])
  gendata <- gendata[!gendata$SAMPLE_ID %in% exclude, ]
   
  ## Read missing data information for individuals:
  ind_df <- read.csv(file = here(paste0("sNMF_runs/miss_data/ind_df_", group, ".csv")), header = TRUE)
  
  ## What samples have less than max_miss missing data?
  inds_to_keep <- ind_df[ind_df$miss_prop <= miss_ind, ]
  
  ## Keep samples with less than max_miss missing data:
  gendata <- gendata[gendata$SAMPLE_ID %in% inds_to_keep$SAMPLE_ID, ]
  
  ## Read sample information:
  sample_info <- read.csv(file = here("sample_information/ddRAD_sample_information.csv"), header = TRUE)
  
  ## Let's only keep the individuals included in the genetic dataset:
  sample_info <- sample_info[sample_info$SAMPLE_ID %in% gendata$SAMPLE_ID, ]
  
  ## How many samples per site?
  site_list <- group_by(sample_info, UPDATED_SP, LAT, LON) %>% group_split()
  n_per_site <- purrr::map(.x = site_list, .f = nrow)
   
  ## We're going to enforce a maximum number of samples per taxon.
  ## This is because sampling imbalance can impair clustering analyses.
  max_n <- 5
  
  ## Sites that have less or equal the maximum number of samples:
  norm_sampled <- site_list[n_per_site <= max_n]
  norm_sampled <- dplyr::bind_rows(norm_sampled)
   
  ## Oversampled sites:
  over_sampled <- site_list[n_per_site > max_n]
  
  ## Rarify oversampled sites:
  over_sampled_filtered <- purrr::map_dfr(over_sampled, sample_n, size = max_n)
   
  ## Combine "normal" and filtered oversampled sites:
  sample_info <- rbind(norm_sampled, over_sampled_filtered)
  
  ## Let's only keep genetic data for the individuals that passed all filters:
  gendata <- gendata[gendata$SAMPLE_ID %in% sample_info$SAMPLE_ID, ]
  
  ## Adding row names:
  row.names(gendata) <- gendata$SAMPLE_ID
  
  ## Remove sample ID column:
  gendata <- subset(gendata, select = -SAMPLE_ID)
  
  ## Saving IDs:
  write.table(row.names(gendata), file = here(paste0("sNMF_runs/samples/samples_", group, "_mi", miss_ind, ".txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ## Writing down genetic data in geno format.
  ## Because we want to compare different alpha values in sNMF, we have to save the geno data once per alpha.
  ## This is due to a limitation of sNMF: saving projects with the same name as the input data files.
  alphas <- c(400, 500)
  purrr::map(alphas, function(a) write.geno(gendata, paste0(runs_path, group, "_mi", miss_ind, "_a", a, ".geno")))
  
  ## Without alpha for downstream analyses:
  write.geno(gendata, here(paste0("sNMF_runs/data/", group, "_mi", miss_ind, ".geno")))

} ## End of function.
  
## PART 4: Function: Run clustering analyses using sNMF ----

## Function:  
run_sNMF <- function(a, cpu, group, maxK, minK, miss_ind, project, rep) {
  
  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5 ; a <- 500 ; minK <- 1 ; maxK <- 20 ; rep <- 3 ; cpu = 4 ; project <- "previous"
  
  ## Running sNMF:
  if (project == "new") {
    project.snmf <- snmf(input.file = paste0(runs_path, "/", group, "_mi", miss_ind, "_a", a, ".geno"), 
                         entropy = TRUE, ploidy = 2, project = "new", #seed = 2020,   
                         CPU = cpu, K = minK:maxK, alpha = a, repetitions = rep) }
    
  ## Loading previously saved sNMF project:
  if (project == "previous") {
    project.snmf <- load.snmfProject(paste0(runs_path, "/", group, "_mi", miss_ind, "_a", a, ".snmfProject")) }
  
  ## Showing summary of project results:
  snmf_summary <- summary(project.snmf)
  snmf_summary
  
  ## Selecting criterion to determine the best K:
  crossEntropy <- snmf_summary$crossEntropy[1,] # Using min cross entropy across runs.
  #crossEntropy <- snmf_summary$crossEntropy[2,] # Using mean cross entropy across runs.
  names(crossEntropy) %<>% gsub(pattern = "K = ", replacement = "", .)
  
  ## Formatting:
  crossEntropy_df <- data.frame(K = names(crossEntropy), ce = crossEntropy)
  crossEntropy_df <- arrange(crossEntropy_df, by = ce)
  
  ## Selecting best K based on minimum (or mean) cross entropy among runs:
  bestK <- as.numeric(as.character(crossEntropy_df$K[which.min(crossEntropy_df$ce)]))
  K <- bestK
 
  ## Plotting mean entropy scores:
  png(filename = here(paste0("sNMF_runs/plots/", group, "_mi", miss_ind, "_a", a, "_cross_entropy.png")))
  plot(y = crossEntropy, x = names(crossEntropy), cex = 1.2, col = "blue", pch = 19, xlab = paste0("K values"), ylab = "Cross-entropy")
  lines(y = crossEntropy, x = names(crossEntropy), col = "blue")
  dev.off() # Saving plot.
  
  ## Getting the cross entropy of all runs for K:
  ce <- cross.entropy(project.snmf, K = K)
    
  ## A few edits:
  ce_df <- as.data.frame(ce)
  ce_df$run <- 1:rep
  ce_df <- ce_df[order(ce_df[, 1]), ]
  names(ce_df)[1] <- "ce"
    
  ## Selecting the run with the lowest cross-entropy for K = best K:
  bestrun <- which.min(ce)
    
  ## Getting the Q matrix:
  qmatrix <- as.data.frame(Q(project.snmf, K = K, run = bestrun))
  
  ## Replace column names:
  names(qmatrix) <- gsub(names(qmatrix), pattern = "V", replacement = "cluster_\\1")
      
  ## Adding individual ID names for plotting:
  individuals <- read.table(file = here(paste0("sNMF_runs/samples/samples_", group, "_mi", miss_ind, ".txt")), header = FALSE)
  qmatrix$ID <- individuals$V1

  ## Save "raw" qmatrix:
  write.csv(qmatrix, file = here(paste0("sNMF_runs/qmatrices/", group, "_mi", miss_ind, "_a", a, "_raw.csv")), row.names = FALSE)
  
} ## End of function

## PART 5: Function: Get cluster assignments ----

## Function:
get_assigments <- function(a, group, miss_ind) {

  ## Testing:
  #group <- "atlas_gr" ; miss_ind <- 0.5 ; a <- 200
  
  ## Read "raw" qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF_runs/qmatrices/", group, "_mi", miss_ind, "_a", a, "_raw.csv")), header = TRUE)
  
  ## Best K:
  qscore_col <- names(qmatrix)[grep(x = names(qmatrix), pattern = "cluster_[0-9]+")]
  K <- length(qscore_col)
  
  ## "Melt" dataframe using to assign samples to clusters based on max qscores:
  qmatrix_melt <- gather(qmatrix, key = cluster_assign, value = coeff, all_of(qscore_col))
  
  ## Assign specimens to cluster based on the highest qscore (coeff) values:
  cluster_assign <- qmatrix_melt %>% group_by(ID) %>% top_n(n = 1, wt = coeff)
    
  ## Combine qmatrix and cluster assignments:
  qmatrix_c <- merge(qmatrix, cluster_assign, by = "ID")
    
  ## Now, let's add locality information (including lat-longs) to the qmatrix:
  sample_info_all <- read.csv(file = here("sample_information/ddRAD_sample_information.csv"), header = TRUE)
  
  ## Let's only keep information for the individuals included in sNMF analyses.
  sample_info <- (sample_info_all[sample_info_all$SAMPLE_ID %in% qmatrix_c$ID, ])
  names(sample_info)[names(sample_info) == "SAMPLE_ID"] <- "ID" ## Renaming ID column for consistency across files.
    
  ## Merging:
  qmatrix_c <- merge(qmatrix_c, sample_info, by = "ID")
  
  ## Save qmatrix:
  write.csv(qmatrix_c, file = here(paste0("sNMF_runs/qmatrices/", group, "_mi", miss_ind, "_a", a, ".csv")), row.names = FALSE)
  
  ## Keeping only columns of interest:
  assignments <- data.frame(ID = qmatrix_c$ID, cluster_assign = qmatrix_c$cluster_assign, coeff = qmatrix_c$coeff)
  names(assignments)[2] <- "cluster"
  
  ## Save assignments:
  write.csv(assignments, file = here(paste0("sNMF_runs/assignments/", group, "_mi", miss_ind, "_a", a, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Using a phylogeny to order samples in structure plots:
  tree <- read.nexus(here("RAxML/Ctenotus.nex"))
  
  ## Using tip labels to guide sample order in plots:
  phylo_order <- as.character(tree$tip.label)
  assignments$ID <- factor(assignments$ID, levels = phylo_order)
  assignments <- arrange(assignments, ID)
  
  ## Saving new sample order: 
  assignments$plot_order <- 1:nrow(assignments)
  
  ## Keeping columns of interest:
  sample_order <- assignments[, c("ID", "plot_order")]
  
  ## Save plot order:
  write.csv(sample_order, file = here(paste0("sNMF_runs/samples/order_", group, "_mi", miss_ind, "_a", a, ".csv")), row.names = TRUE)
  
} ## End of function.

## PART 6: Function: Color palette for plots ----

## Color palette:
get_palette <- function(K) {
  colors <- stepped2()[seq(1, 18, 1.25)] 
  color_mask <- ceiling(x = seq(from = 1, to = length(colors), length.out = K))
  palette <- colors[color_mask]
  return(palette)
}

## PART 7: Function: Plotting ancestry coefficients ----

## Creating function to plot with ggplot:
plot_sNMF_bars <- function(a, group, minK, miss_ind) {

  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5 ; a <- 500
  
  ## Print status:
  print(paste0("Now making barplot for ", group, "!"))
   
  ## Read corresponding qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF_runs/qmatrices/", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
  
  ## Best-fit K:
  K <- length(unique(qmatrix$cluster_assign))
  
  ## Load plot order:
  plot_order <- read.csv(file = here(paste0("sNMF_runs/samples/order_", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
  
  ## Merge to plot order:
  qmatrix <- merge(qmatrix, plot_order, by = "ID")

  ## Arrange:
  qmatrix <- arrange(qmatrix, plot_order)
  
  ## Set factors levels to preserve order in plots:
  qmatrix$ID <- factor(qmatrix$ID, levels = unique(qmatrix$ID))

  ## Exclude admixed samples to set factor levels:
  qmatrix_l <- qmatrix[qmatrix$coeff > 0.95, ]

  ## "Melt":
  qmatrix_m <- gather(qmatrix, key = sNMF_cluster, value = qscores, 2:(K+1))
        
  ## Set factors levels to preserve order in plots:
  qmatrix_m$sNMF_cluster <- factor(qmatrix_m$sNMF_cluster, levels = unique(qmatrix_l$cluster_assign))
    
  ## Rounding qscores:
  qmatrix_m$qscores <- round(qmatrix_m$qscores, digits = 2)
    
  ## Color palette:
  palette <- get_palette(K)

  ## Finally, creating stacked bar plot of ancestry coefficients:
  plot <- ggplot(data = qmatrix_m, aes(y = ID)) +
      
    ## Adding bars that represent ancestry coefficients:
    geom_bar(aes(x = qscores, fill = sNMF_cluster), 
             color = "white", ## Color of bar lines.
             size = 0.01, ## Thickness of bar lines.
             stat = "identity", position = "fill",
	           show.legend = FALSE) +
        
    ## Filling bars by cluster assigment:
    scale_fill_manual(values = palette) + 
          
    ## Adjusting labels:
    labs(x = "Ancestry proportions", y = "") +
          
    ## Adjusting limits of the x axis:
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = round(c(0, 1), 0)) +
            
    ## Changing theme:
    theme_minimal() +
  
    ## Theme parameters:
    theme(
          axis.text.y = element_text(color = "gray30", angle = 0, vjust = 0.5, hjust = 1, size = 5),
          #axis.text.y = element_blank(), ## Removing IDs from plot.
          axis.title.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid = element_blank(), 
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 0, margin = margin(0, 0, 0, 0)),
          plot.margin = margin(r = 15, l = 10, t = 0, b = 0),
          plot.title = element_blank())
	    
  ## Return:
  return(plot)
  
} ## End of function.

## PART 8: Function: Plotting maps for the resulting clusters ----

## Function:
plot_maps <- function(a, group, minK, miss_ind) {
  
  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5 ; a <- 500
  
  ## Print status:
  print(paste0("Now making maps for ", group, "!"))
  
  ## Let's start by setting up a few map features. Importing world map:
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  ## Selecting the part that represents Australia:
  AUS <- filter(world, admin == "Australia")
  
  ## Restrict geographic area using a box:
  bbox <- st_bbox(AUS)
  bbox[["xmin"]] = 110.5
  bbox[["xmax"]] = 156
  bbox[["ymin"]] = -46
  bbox[["ymax"]] = -8.5
  
  ## Read corresponding qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF_runs/qmatrices/", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
  
  ## Best-fit K:
  K <- length(unique(qmatrix$cluster_assign))
  
  ## Load plot order:
  plot_order <- read.csv(file = here(paste0("sNMF_runs/samples/order_", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
  
  ## Merge to plot order:
  qmatrix <- merge(qmatrix, plot_order, by = "ID")

  ## Arrange:
  qmatrix <- arrange(qmatrix, desc(plot_order))

  ## Exclude admixed samples to set factor levels:
  qmatrix_l <- qmatrix[qmatrix$coeff >= 0.95, ]

  ## Set order of factors to preserve sample order in plots:
  qmatrix$cluster_assign <- factor(qmatrix$cluster_assign, levels = unique(qmatrix_l$cluster_assign))
      
  ## Change order of columns to keep cluster color order in map plots:
  qmatrix[c(paste0("cluster_", 1:K))] <- qmatrix[c(paste0("cluster_", 1:K))][levels(qmatrix$cluster_assign)]
  
  ## Make qscores numeric:
  qmatrix[c(paste0("cluster_", 1:K))] <- apply(qmatrix[c(paste0("cluster_", 1:K))], 2, as.numeric)
  
  ## Creating labels for map facets:
  label_df <- data.frame(cluster = qmatrix$cluster_assign, taxon = qmatrix$UPDATED_SP)
  #label_df$taxon <- paste0("C. ",  label_df$taxon) 
  table_lb <- data.frame(table(label_df))
  table_lb <- arrange(table_lb, cluster)
  label_df <- table_lb %>% group_by(cluster) %>% top_n(Freq, n = 1)
  label_df <- label_df %>% group_by(cluster) %>% sample_n(size = 1)
  
  ## Create labels for duplicated taxon name across clusters:
  label_df$label <- as.character(label_df$taxon)
  duplicated <- names(table(label_df$taxon))[table(label_df$taxon) >= 2]
  for (taxon in duplicated) {
    n_taxon <- table(label_df$taxon)[taxon]
    new_labels <- paste0(taxon, " (", 1:n_taxon, ")")
    for (plabel in 1:length(new_labels)) {
      label_df$label[label_df$taxon == taxon][plabel] <- as.character(new_labels[plabel]) 
    }
  }
  
  ## Save:
  write.csv(label_df, file = here(paste0("sNMF_runs/samples/labels_", group, "_mi", miss_ind, "_a", a, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Add labels to qmatrix:
  qmatrix$label <- factor(qmatrix$cluster_assign, levels = label_df$cluster, labels = label_df$label)
  
  ## Color palette:
  palette <- get_palette(K)
  
  ## Finally, plot map:
  plot <- ggplot() +
        
    ## Adding baseline map:
    geom_sf(data = world, fill = "gray99", color = "gray30", size = 0.15) +
  
    ## Adding lat-longs for sampled individuals:
    geom_point(data = qmatrix, aes(x = LON, y = LAT, fill = cluster_assign), 
               shape = 21, alpha = 0.7,
               size = 2, color = "black", stroke = 0.2) +
                  
    ## Set map boundaries:
    coord_sf(xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")]) +
        
    ## Setting a bunch of aesthetic parameters:
    theme_void() +
    guides(fill = "none", size = "none", alpha = "none", shape = "none") +
    scale_fill_manual(values = rev(palette)) +
        
    ## Setting up facets by cluster:
    facet_wrap(~label, ncol = 2) +
        
    ## Setting some other elements:
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
          strip.text = element_text(size = 9, face = "italic", margin = margin(t = 0, r = 0, b = 2, l = 0)),
          strip.background = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank())
          
  ## Return:
  return(plot)
    
} ## End of function.

## PART 9: Function: Plot trees ----

## Function:
plot_tree <- function(a, group, miss_ind) {

  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5 ; a <- 500
  
  ## Print status:
  print(paste0("Now plotting the tree for ", group, "!"))
  
  ## Load sample info edited previously:
  assignments <- read.csv(file = here(paste0("sNMF_runs/assignments/", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
    
  ## List of samples in each taxon:
  t_list <- split(x = assignments, f = assignments$cluster)
  t_list <- purrr::map(t_list, .f = function(x) as.character(x$ID))
  
  ## Read tree:
  tree <- read.nexus(paste0(path, "RAxML/Ctenotus.nex"))
   
  ## Clade span:
  clade_samples <- as.character(unlist(t_list))
  tree_span <- clade_samples[clade_samples %in% tree$tip.label] 
  
  ## Keep only clade including focal taxa:
  node <- findMRCA(tree = tree, type = "node", tips = tree_span)
  tree <- extract.clade(tree, node = node)
    
  ## If needed, keep only samples used in sNMF analyses:
  remove_tips <- tree$tip.label[!tree$tip.label %in% clade_samples]
  tree <- drop.tip(tree, tree$tip.label[match(remove_tips, tree$tip.label)])
      
  ## Set clade color order:
  plot_order <- read.csv(file = here(paste0("sNMF_runs/samples/order_", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
  pure_ind <- assignments[assignments$coeff >= 0.95, ]
  order_df <- merge(plot_order, pure_ind, "ID")
  order_df <- arrange(order_df, plot_order)
  levels <- as.character(unique(order_df$cluster))
  t_list_o <- t_list[levels]
  
  ## Group tips by taxon to color tips in ggtree. We'll have to first rename clades to keep color order.
  get_names <- function(t) { paste0(letters[t], "_", names(t_list_o[t])) }
  new_names <- unlist(purrr::map(1:length(t_list_o), get_names))
  names(t_list_o) <- new_names
  tree <- groupOTU(tree, t_list_o)
  
  ## Color palette:
  K <- length(unique(assignments$cluster))
  palette <- get_palette(K)
  palette <- c("gray90", palette)
  
  ## Plot tree:
  plot <- ggtree(tree, size = 0.5, aes(color = group)) + ## Size = branch line thickness.
    
    ## Editing tree tips:
    scale_color_manual(values = palette) +
    
    ## Other edits:
    theme(legend.position = "none",
          plot.title = element_blank(),
          plot.caption = element_blank())
          
  ## Return:
  return(plot)

} ## End of function.

## PART 10: Function: Running PCA on the SNPs ----

## Function:
run_PCA <- function(a, group, miss_ind) {
  
  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5 ; a <- 500
  
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
  assignments <- read.csv(header = TRUE, file = here(paste0("sNMF_runs/assignments/", group, "_mi", miss_ind, "_a", a, ".csv")))
  
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

## PART 11: Function: Plotting PCA results ----

## Plotting PCA results:
plot_PCA <- function(a, group, miss_ind) {
  
  ## Testing:
  #group <- "inornatus_gr" ; miss_ind <- 0.5 ; a <- 500
  
  ## Run PCA:
  PCA_results <- run_PCA(a = a, group = group, miss_ind = miss_ind)
  
  ## Data we'll use:
  pcadata <- PCA_results$pcadata[!is.na(PCA_results$pcadata$cluster), ]
  pcadata <- pcadata[, c("ID", "cluster", "label", c(paste0("PC", rep(1:4)))) ]
  
  ## Ordering factors to keep order in plot:
  pcadata$cluster <- factor(pcadata$cluster, levels = unique(pcadata$cluster))
  
  ## Best fit-K:
  ## Read corresponding qmatrix:
  qmatrix <- read.csv(file = here(paste0("sNMF_runs/qmatrices/", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE)
  K <- length(unique(qmatrix$cluster_assign))
  
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
				 expand = unit(2, "mm"), size = 0.1, alpha = 0.3) +

      ## Configuring point size and shape:
      geom_point(aes(x = .data[[PCx]], y = .data[[PCy]], fill = cluster), size = 1.5, shape = 21, 
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
      theme(panel.border = element_rect(size = 0.5, colour = "gray30"),
            axis.text.x = element_text(size = 8, hjust = 0.5, margin = margin(t = 2, r = 0, b = 0, l = 0)),
            axis.text.y = element_text(size = 8, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
            axis.title.x = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 10, hjust = 0.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
            panel.grid = element_blank(), ## Removing background grid.
            axis.ticks = element_line(size = 0.25, color = "gray30"), ## Making axis thicks thicker.
            axis.ticks.length = unit(0.1, "cm")) ## Making axis ticks longer.
            
    ## Function will return:
    return(plot)
    
  } ## End of function.
  
  ## Plot for PC combinations while getting rid of redundant axis labels:
  plot <- plot_PCaxes(PCx = "PC1", PCy = "PC2")
  
  ## Remove pca-related files:
  unlink(here("*.pcaProject"))
  unlink(here("sNMF_runs/data/*.lfmm"))
  unlink(here("sNMF_runs/data/*.pca"), recursive = TRUE)
  
  ## Return:
  return(plot)
  
} ## End of function.

## PART 12: Function: Running it all for each species group ----

## Function:
run_group <- function(a, group, miss_ind, project) {

  ## Testing:
  #group <- "atlas_gr" ;  miss_ind <- 0.5 ;  a <- 200 ;  project <- "previous"
  
  ## Filter the data?
  #if (project == "new") { prep_data(group = group, miss_ind = miss_ind) }
        
  ## Run sNMF:
  #run_sNMF(a = a, cpu = cpu, group = group, maxK = maxK, minK = minK, miss_ind = miss_ind, project = project, rep = rep)
  
  ## Get assignments:
  #get_assigments(a = a, group = group, miss_ind = miss_ind)
  
  ## Plot ancestry coefficients:
  plot_b <- plot_sNMF_bars(a = a, group = group, minK = minK, miss_ind = miss_ind)
  
  ## Plot maps:
  plot_m <- plot_maps(a = a, group = group, minK = minK, miss_ind = miss_ind)
  
  ## Plot phylogeny:
  plot_d <- plot_tree(a = a, group = group, miss_ind = miss_ind)
  
  ## Plot genetic PCA:
  plot_p <- plot_PCA(a = a, group = group, miss_ind = miss_ind)
  
  ## Print status:
  print(paste0("Now combining plots for ", group, "!"))
  
  ## Combining plots:
  plot_mp <- ( plot_m / plot_p ) + plot_layout(heights = c(3, 1))
  plot_f <- ( plot_d + plot_b + plot_mp ) + plot_layout(widths = c(1.15, 1, 2)) +
  
  ## Plot annotation: 
  plot_annotation(tag_levels = "A",
                  title = paste0("Ctenotus ", gsub(x = group, pattern = "_gr", replacement = ""), " group")) &
  theme(plot.tag = element_text(size = 18),
        plot.tag.position = c(-0.005, 0.99),
        plot.title = element_text(size = 18, hjust = 0.5, face = "italic", margin = margin(t = 10, r = 0, b = 10, l = 0)))
  
  ## Saving plot:
  ggsave(plot = plot_f, width = 18, height = 19, units = "cm", limitsize = FALSE, device = "pdf",
         filename = here(paste0("sNMF_runs/plots/", group, "_mi", miss_ind, "_a", a, ".pdf")))
  
  ## Print status:
  print(paste0("All analyses for ", group, " completed!"))
  
} ## End of function.

## PART 13: Running it all for all species groups ----

## Parameters:
minK <- 1
maxK <- 20
cpu <- 32
rep <- 50
project <- "previous"
groups <- c("atlas_gr", "colletti_gr", "inornatus_gr", "leonhardii_gr", "essingtonii_gr", "pantherinus_gr", "schomburgkii_gr", "taeniolatus_gr")

## Estimating missing data:
#purrr::map(groups, calc_miss)

## Run:
for (a in c(400, 500)) {
  purrr::map(groups, run_group, a = a, miss_ind = 0.5, project = project)
}

## PART 14: Function: Combining assigments across groups ----

## Function:
join_assignments <- function(group, miss_ind) {
  
  ## Testing:
  #group <- "atlas_gr"
    
  ## alpha parameter:
  if (group == "atlas_gr") { a <- 500 }
  if (group == "colletti_gr") { a <- 400 }
  if (group == "essingtonii_gr") { a <- 400 }
  if (group == "inornatus_gr") { a <- 400 }
  if (group == "leonhardii_gr") { a <- 500 }
  if (group == "pantherinus_gr") { a <- 400 }
  if (group == "schomburgkii_gr") { a <- 500 }
  if (group == "taeniolatus_gr") { a <- 500 }
    
  ## Load assignments:
  assignments <- read.csv(header = TRUE, file = here(paste0("sNMF_runs/assignments/", group, "_mi", miss_ind, "_a", a, ".csv")), stringsAsFactors = FALSE)
  assignments$group <- group
  assignments <- dplyr::select(assignments, -coeff)
  assignments$cluster <- paste0(assignments$group, "_", assignments$cluster)
    
  ## Load taxon labels for final OTU scheme:
  label_df <- read.csv(file = here(paste0("sNMF_runs/samples/labels_", group, "_mi", miss_ind, "_a", a, ".csv")), header = TRUE, stringsAsFactors = FALSE)
  label_df$OTU <- paste0(group, "_", label_df$cluster)
  
  ## Save labels:
  write.csv(label_df, file = here(paste0("sNMF_runs/samples/labels_", group, "_mi", miss_ind, "_delimited_OTUs.csv")), row.names = FALSE, quote = FALSE)
    
  ## Return:
  return(as.data.frame(assignments))
    
} ## End of function.
  
## Apply:
assign_df <- purrr::map_dfr(.x = groups, .f = join_assignments, miss_ind = 0.5)

## Save assignments:
write.csv(assign_df, file = here("outputs/assignments_Ctenotus_OTUs.csv"), row.names = FALSE, quote = FALSE)

## End of script.
