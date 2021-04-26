####################
### The goals of this script are:
### To estimate expected heterozygosity to Ctenotus OTUs.
### To test for associations between range size, IBD slope, and heterozygosity.
### To make plots of these relationships.
### Checked for functionality on April 23rd 2021.

## PART 1: Getting ready ----

## Packages:
library(patchwork)
library(pegas)
library(tidyverse)
library(vcfR)

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## rhinellax.
runs_path <- "/home/ivan/sNMF_runs/" ## rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## PART 2: Function: Estimating heterozygosity ----

## OTU assignments:
assignments <- read.csv(file = here("outputs/assignments_Ctenotus_OTUs.csv"), header = TRUE, stringsAsFactors = FALSE)
assignments <- assignments[c("ID", "cluster")] 
assignments <- arrange(assignments, ID)

## Creating function to calculate H_ad:
calc_H_ad <- function(group){
  
  ## Testing:  
  #group <- "schomburgkii_gr"

  ## Status:
  print(paste0("Now estimating heterozigosity for ", group, "!"))
  
  ## Reading vcf using vcfR:
  gendata_vcfR <- read.vcfR(file = here("genetic_data_filtered_vcftools/", group, "/usnps.vcf"))
  
  ## Converting to genind:
  genind <- vcfR2genind(gendata_vcfR)
  
  ## What samples are present in both the VCF and delimitation analyses?
  overlapping <- intersect(indNames(genind), as.character(assignments$ID))
  
  ## Select samples:
  genind2 <- genind[indNames(genind) %in% overlapping]
  
  ## Are individuals (and their order) the same?
  table(indNames(genind2) == overlapping)
  
  ## Assigning individuals to populations:
  pop(genind2) <- assignments$cluster[assignments$ID %in% overlapping]
  
  ## Calculating expected heterozygosity using adegenet:
  Expected_heterozygosity <- Hs(genind2)
  gen_stats <- data.frame(Expected_heterozygosity = Expected_heterozygosity)
  gen_stats$OTU <- row.names(gen_stats) 
  
  ## Return:
  return(gen_stats)
    
} ## End of function.
  
## Groups:
groups <- c("atlas_gr", "colletti_gr", "inornatus_gr", "leonhardii_gr", "essingtonii_gr", "pantherinus_gr", "schomburgkii_gr", "taeniolatus_gr")

## Apply function to groups:  
H_df <- purrr::map_df(groups, calc_H_ad)
H_df <- H_df[c("OTU", "Expected_heterozygosity")]

## Save:
write.csv(H_df, file = here("outputs/expected_heterozygosity_Ctenotus_OTUs.csv"), row.names = FALSE)

## PART 3: Linear regressions ----

## Range size and IBD slope info:
slope_df <- read.csv(header = TRUE, file = here("outputs/IBD_slopes_Ctenotus_OTUs.csv"))
rs_nb_df <- read.csv(header = TRUE, file = here("outputs/range_size_climatic_niche_breadth_Ctenotus_OTUs.csv"))
names(rs_nb_df)[names(rs_nb_df) == "group"] <- "OTU"

## Merge:
estimates_df <- full_join(slope_df, H_df, by = "OTU")
estimates_df <- full_join(estimates_df, rs_nb_df, by = "OTU")

## Keep only OTUs with slope and range information:
estimates_df <- estimates_df[!is.na(estimates_df$slope) &
                             !is.na(estimates_df$log_rs), ]

## Linear regression: Effect of expected heterozygosity on IBD slope:
lm_1 <- lm(slope ~ Expected_heterozygosity, data = estimates_df) ; summary(lm_1)
plot(y = estimates_df$slope, x = estimates_df$Expected_heterozygosity, xlab = "Expected heterozygosity", ylab = "Slope of IBD")

## Linear regression: Effect of range size on expected heterozygosity:
lm_2 <- lm(Expected_heterozygosity ~ log_rs, data = estimates_df) ; summary(lm_2)
plot(x = estimates_df$log_rs, y = estimates_df$Expected_heterozygosity, ylab = "Expected heterozygosity", xlab = "Range size")

## PART 4: Plot linear regressions ----

## Column names:
names(estimates_df)[names(estimates_df) == "log_rs"] <- "Range size"
names(estimates_df)[names(estimates_df) == "slope"] <- "Slope of IBD"
names(estimates_df)[names(estimates_df) == "Expected_heterozygosity"] <- "Expected heterozygosity"
names(estimates_df)[names(estimates_df) == "Observed_heterozygosity"] <- "Observed heterozygosity"

## Function: Plot:
plot_H <- function(yvar, xvar) {

  ## Plot:
  plot <- ggplot(data = estimates_df, aes(y = .data[[yvar]], x = .data[[xvar]])) +
    
    ## Adding regression line:
    geom_smooth(method = "lm", color = "blue", size = 0.75, se = FALSE) +
    
    ## Scatter plot:    
    geom_point(size = 2.5, alpha = 0.5, color = "gray30") +
    
    ## Title and axes:
    scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    
    ## Changing overall theme:
    theme_bw() +
    
    ## Other params:
    theme(
      plot.margin = margin(t = 5, b = 5, l = 5, r = 5),
      axis.title.x = element_text(size = 14, margin = margin(t = 10, b = 0, l = 0, r = 0)),
      axis.title.y = element_text(size = 14, margin = margin(t = 0, b = 0, l = 0, r = 5)),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12), 
      legend.text = element_text(size = 10),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.border = element_rect(size = 1, colour = "gray20"), ## Changing box around plot.
      axis.ticks = element_line(size = 0.75, color = "gray20"), ## Making axis ticks thicker.
      panel.grid = element_blank()) ## Removing background grid.
  
  ## Return:
  return(plot)
  
} ## End of function.

## Run:
plot_A <- plot_H(yvar = "Expected heterozygosity", xvar = "Range size")
plot_B <- plot_H(yvar = "Slope of IBD", xvar = "Expected heterozygosity")

## Combine:
Fig_S6 <- (plot_A + plot_B) + 
  
  ## Annotation: 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22),
        plot.tag.position = c(0.024, 0.98))

## Save plot:
ggsave(plot = Fig_S6, height = 5, width = 12, filename = here("plots/FigS6.jpeg"))

## End of script.
