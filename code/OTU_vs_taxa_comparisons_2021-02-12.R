####################
### The goals of this script are:
### To compare OTU- vs. taxon-based estimates of range size, climatic niche breadth, and slope of IBD in Ctenotus.
### To plot these comparisons.
### Checked for functionality on April 23rd 2021.

## PART 1: Getting ready ----

## Loading packages:
library(patchwork)
library(tidyverse)

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## PART 2: Pairing OTU and taxon assigments ----

## Load sample info and assigment:
assignments <- read.csv(file = here("outputs/assignments_Ctenotus_OTUs.csv"), header = TRUE)
sample_info <- read.csv(file = here("sample_information/ddRAD_sample_information.csv"), header = TRUE)
sample_info <- sample_info[c("SAMPLE_ID", "UPDATED_SP")]
names(sample_info)[names(sample_info) == "SAMPLE_ID"] <- "ID" ## Change column name.

## Pairing OTU and taxon assignment:
paired_info <- merge(assignments, sample_info, by = "ID")
paired_info <- paired_info[c("UPDATED_SP", "cluster")]
names(paired_info) <- c("taxon", "OTU")

## Keep unique combinations:
paired_info <- paired_info %>% group_by(taxon, OTU) %>% sample_n(1)

## PART 3: Pairing OTU- and taxon-based estimates ----

## Load range size and climatic niche breadth estimates:
rs_nb_OTUs <- read.csv(header = TRUE, file = here("outputs/range_size_climatic_niche_breadth_Ctenotus_OTUs.csv"))
rs_nb_taxa <- read.csv(header = TRUE, file = here("outputs/range_size_climatic_niche_breadth_Ctenotus_taxa.csv"))
names(rs_nb_OTUs)[names(rs_nb_OTUs) == "group"] <- "OTU" ## Change column name.
names(rs_nb_taxa)[names(rs_nb_taxa) == "group"] <- "taxon" ## Change column name.

## Add OTU info to range size and niche breadth data:
rs_nb_OTUs <- full_join(rs_nb_OTUs, paired_info, by = "OTU")
rs_nb_OTUs <- rs_nb_OTUs[c("OTU", "taxon", "log_rs", "log_nb")] ## Log-transformed.
names(rs_nb_OTUs) <- c("OTU", "taxon", "rs_OTU", "nb_OTU")

## Add taxon info to range size and niche breadth data:
rs_nb_taxa <- full_join(rs_nb_taxa, paired_info, by = "taxon")
rs_nb_taxa <- rs_nb_taxa[c("taxon", "OTU", "log_rs", "log_nb")] ## Log-transformed.
names(rs_nb_taxa) <- c("taxon", "OTU", "rs_taxon", "nb_taxon")

## Merge:
rs_nb_df <- full_join(rs_nb_taxa, rs_nb_OTUs, by = c("OTU", "taxon"))

## Now add taxon info to OTU slope data:
slope_OTU <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_OTUs.csv")))
slope_OTU <- full_join(slope_OTU, paired_info, by = "OTU")
names(slope_OTU)[names(slope_OTU) == "slope"] <- "slope_OTU"

## Add OTU info to taxon slope data:
slope_taxa <- read.csv(header = TRUE, file = here(paste0("outputs/IBD_slopes_Ctenotus_taxa.csv")))
names(slope_taxa)[names(slope_taxa) == "OTU"] <- "taxon"
slope_taxa <- full_join(slope_taxa, paired_info, by = "taxon")
names(slope_taxa)[names(slope_taxa) == "slope"] <- "slope_taxon"

## Combine slope data:
slopes <- full_join(slope_OTU, slope_taxa, by = c("OTU", "taxon"))

## Merge:
estimates_df <- full_join(rs_nb_df, slopes, by = c("taxon", "OTU"))

## Keep only taxon and OTUs with no missing estimates:
estimates_df <- estimates_df[!is.na(estimates_df$rs_OTU) & !is.na(estimates_df$rs_taxon) &
                             !is.na(estimates_df$nb_OTU) & !is.na(estimates_df$nb_taxon) &
                             !is.na(estimates_df$slope_OTU) & !is.na(estimates_df$slope_taxon), ]

## Save:
write.csv(estimates_df, file = here("outputs/OTU_vs_taxon_based_estimates_Ctenotus.csv"), row.names = FALSE)

## PART 4: Function: Plot ranges and niches using OTU vs. taxa ----

## Function:
plot_OTU_vs_taxa <- function(var) {
  
  ## Testing:
  #var <- "Range size"
  
  ## Data to use:
  if (var == "Range size") { estimates_df <- estimates_df[c("rs_taxon", "rs_OTU")] }
  if (var == "Climatic niche breadth") { estimates_df <- estimates_df[c("nb_taxon", "nb_OTU")] }
  if (var == "Slope of IBD") { estimates_df <- estimates_df[c("slope_taxon", "slope_OTU")] }
  names(estimates_df) <- c("var_taxa", "var_OTUs")
  
  ## Plot:
  plot <- ggplot(data = estimates_df) +
    
    ## Adding 1:1 line:
    geom_abline(intercept = 0, slope = 1, color = "blue", size = 0.75) +  
    
    ## Scatter plot:    
    geom_point(aes(y = var_taxa, x = var_OTUs), size = 2.5, color = "gray50", alpha = 0.8) +
    
    ## Title and axes:
    xlab(label = paste0(var, " (OTUs)")) +
    ylab(label = paste0(var, " (taxa)")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    
    ## Changing overall theme:
    theme_bw() +
    
    ## Other params:
    theme(
      plot.margin = margin(t = 0, b = 0, l = 0, r = 20),
      axis.title.x = element_text(size = 18, margin = margin(t = 10, b = 0, l = 0, r = 0)),
      axis.title.y = element_text(size = 18, margin = margin(t = 0, b = 0, l = 0, r = 5)),
      axis.text = element_text(size = 16),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.border = element_rect(size = 1, colour = "gray20"), ## Changing box around plot.
      axis.ticks = element_line(size = 0.75, color = "gray20"), ## Making axis ticks thicker.
      panel.grid = element_blank()) ## Removing background grid.
  
  ## Return:
  return(plot)
  
} ## End of function.

## PART 5: Plot and combine plots ----

## Plot:
plot_a <- plot_OTU_vs_taxa(var = "Range size")
plot_b <- plot_OTU_vs_taxa(var = "Climatic niche breadth")
plot_c <- plot_OTU_vs_taxa(var = "Slope of IBD")

## Combine plots:
Fig3 <- plot_a + plot_b + plot_c +
  
  ## Annotation: 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22),
        plot.tag.position = c(0.025, 0.98))

## Save plot:
ggsave(plot = Fig3, height = 5, width = 16, dpi = 500, device = "jpeg", units = "in", filename = here("plots/Fig3.jpeg"))

## PART 6: Number of intra-taxon OTU vs. range size ----

## Import biome information:
biome_df <- read.csv(header = TRUE, file = here("outputs/biomes_Ctenotus_OTUs.csv"))
biome_df$biome <- gsub(x = biome_df$biome, pattern = "desert", replacement = "Desert")
biome_df$biome <- gsub(x = biome_df$biome, pattern = "mediter_woodland", replacement = "Mediterranean woodlands")
biome_df$biome <- gsub(x = biome_df$biome, pattern = "temp_forest", replacement = "Temperate forests")
biome_df$biome <- gsub(x = biome_df$biome, pattern = "trop_grassland", replacement = "Tropical grasslands")

## Add biome information to estimates of range size:
range_size_df <- full_join(estimates_df, biome_df, by = "OTU")

## OTU by taxa:
n_OTU <- data.frame(table(range_size_df$taxon))
names(n_OTU) <- c("taxon", "n_OTU")

## Subset:
range_size_df <- range_size_df[c("taxon", "rs_taxon", "biome", "biome_prop")]
range_size_df <- arrange(range_size_df, taxon)

## Assign taxa to biomes based on the highest biome proportion:
range_size_df <- range_size_df %>% group_by(taxon) %>% top_n(n = 1, wt = biome_prop)
range_size_df <- range_size_df[!duplicated(range_size_df), ] ## Removing duplicates.

## Merge:
range_size_df <- merge(range_size_df, n_OTU, by = "taxon")

## Remove missing:
range_size_df <- range_size_df[!is.na(range_size_df$rs_taxon), ]

## Log-transform::
range_size_df$rs_taxon <- log10(range_size_df$rs_taxon) 

## Basic plot:
plot(n_OTU ~ rs_taxon, data = range_size_df)

## Linear regression:
lm_res <- lm(n_OTU ~ rs_taxon, data = range_size_df) ; summary(lm_res)

## Number of OTUs per biome:
#ANOVA <- lm(n_OTU ~ biome, data = range_size_df) ; summary(ANOVA)

## PART 7: Plotting the number of intra-taxon OTUs vs. taxon range size ----

## Plot:
plot <- ggplot(data = range_size_df) +
  
  ## Scatter plot:    
  geom_point(aes(y = n_OTU, x = rs_taxon, fill = biome), size = 2.5, alpha = 0.5, color = "black", shape = 21) +
  
  ## Title and axes:
  guides(fill = guide_legend(title = "Most frequent biome for taxon:", )) +
  xlab(label = "Taxon range size (log-transformed)") +
  ylab(label = "Number of intra-taxon OTUs") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_fill_manual(values = c("red", "green", "yellow", "blue")) +
  
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

## Save plot:
ggsave(plot = plot, height = 5, width = 8, filename = here("plots/FigS5.jpeg"), device = "jpeg", units = "in", dpi = 200)

## End of script.
