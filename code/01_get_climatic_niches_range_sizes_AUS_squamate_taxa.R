####################
### The goals of this script are:
### To estimate climatic niche breadths for the Australian squamate taxa.
### To test for significant relationships between climatic niche breadth and geographic range size.
### To plot Fig. 1 of the manuscript based on the resulting estimates.
### Checked for functionality on April 20th 2021.
### Written by Ivan Prates (ivanprates.org).

## PART 1: Getting ready ----

## Loading packages:
require(crayon)
require(rgdal)
require(ggrepel)
require(hypervolume)
require(patchwork)
require(raster)
require(rnaturalearth)
require(rnaturalearthdata)
require(rgeos)
require(scico)
require(sf)
require(sp)
require(spThin)
require(tidyverse)

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## Creating folders to save results:
dir.create(path = here("plots"))
dir.create(path = here("climatic_data_Australian_taxa"))

## PART 2: Load taxon ranges and environmental data ----

## Read shape file with taxon distribution polygons:
shape_all <- read_sf(dsn = paste0(path, "GIS/Meiri_etal_ranges/GARD1.1_dissolved_ranges"), layer = "modeled_reptiles")
print("All reptile shapes loaded!")

## Read list of Australian taxa to consider:
AUS_taxa <- read.csv(here("outputs/Australian_squamate_taxa_list.csv"), header = TRUE)

## Keep only Australian taxa:
shape_AUS <- shape_all[shape_all$Binomial %in% AUS_taxa$taxon, ]

## Remove to free up memory:
rm(shape_all)

## Converting to spatial polygon class:
shape_used <- as(shape_AUS, "Spatial")
print("Shapes converted to 'spatial' class!")

## Remove to free up memory:
rm(shape_AUS)

## List environmental layers:
layers_clim <- list.files(path = paste0(path, "GIS/chelsa_bioclim"),  pattern = '*.tif', full.names = TRUE)

## If other bio variables are included in this folder, keep only bio1, bio4, bio12, and bio15:
#layers_clim <- layers_clim[c(1, 4, 12, 15)]

## Stacking variables:
layers_clim <- raster::stack(layers_clim)

## PART 3: Function: Extract environmental data from random points for each taxon ----

## Function:
get_clim_data <- function(taxon) {

  ## For testing purposes:
  #taxon = "Ctenotus robustus"
  
  ## Inform progress:
  cat(blue(paste0("\n", "Now extracting climate data for ", taxon, ":\n")))
  
  ## Extracting polygon for taxon:
  polygon <- shape_used[shape_used$Binomial == taxon, ] ## If using spatialpolygondataframe.
  
  ## Setting minimum range size:
  if (polygon$Area < 20) { extracted <- data.frame()
  } else {

    ## Extract random point from the spatial polygon:
    random_pts <- spsample(x = polygon, n = 2000, type = "random")
    
    ## Formatting points for spThin:
    random_pts_form <- as.data.frame(random_pts@coords)
    random_pts_form$taxon = taxon
    
    ## Thinning points to ensure a minimum distance in km:
    thinned_pts <- thin(loc.data = random_pts_form, 
                        lat.col = "y", long.col = "x", spec.col = "taxon", 
                        thin.par = 5, ## Distance between points in km.
                        reps = 1, ## Number of replicates. One is enough for our purposes.
                        locs.thinned.list.return = TRUE, write.files = FALSE, write.log.file = FALSE)
    
    ## List to dataframe:
    thinned_pts_df <- plyr::ldply(thinned_pts, data.frame)
    
    ## How many points?
    n_pts <- dim(thinned_pts_df)[1]
  
    ## Extract raster to points:
    extracted <- raster::extract(layers_clim, thinned_pts_df)
  
    ## Get rid of rows with missing data:
    extracted <- as.data.frame(na.omit(extracted))
    
    ## Taxon:
    extracted$taxon <- taxon
  
  } ## Close if statement.
  
  ## Save:
  write.csv(extracted, file = here(paste0("climatic_data_Australian_taxa/", taxon, ".csv")), row.names = FALSE)
  
} ## End of function.

## Get clim data:
purrr::map_df(AUS_taxa$taxon, get_clim_data)

## PART 4: Transform climatic data for niche estimation ----

## Combining data:
clim_data <- purrr::map_df(AUS_taxa$taxon, function(taxon) read.csv(file = here(paste0("climatic_data_Australian_taxa/", taxon, ".csv")), header = TRUE, stringsAsFactors = FALSE))

## Tranform the data using a z-scores transformation.
clim_data[1:4] <- scale(clim_data[1:4], center = FALSE, scale = TRUE)

## Save:
write.csv(clim_data, file = here("climatic_data_Australian_taxa/Australian_squamate_taxa_climatic_data_for_niche_estimation.csv"), row.names = FALSE)

## PART 5: Function: Estimate climatic niches ----

## Function: generate random points and extract environmental data:
get_hypervolume <- function(taxon){
  
  ## Testing:
  #taxon <- "Acanthophis antarcticus"
  
  ## Inform progress:
  cat(blue(paste0("\n", "Now estimating niche for ", taxon, ":\n")))
  
  ## Extracting data for taxon:
  taxon_clim_data <- clim_data[clim_data$taxon == taxon, ] ## If using spatialpolygondataframe.
  taxon_clim_data <- subset(taxon_clim_data, select = -c(taxon))
  
  ## Test if all values are the same for a variable, that is, if there's no environ variation.
  unique_vals <- function(cl){length(unique(cl)) > 1} ## Function: Is there more than one value?
  test_unique <- apply(taxon_clim_data, 2, unique_vals) ## Test by applying function.
  
  ## Run next steps only if all variables are not invariable.
  if(length(test_unique[test_unique == "FALSE"]) > 0) {
    volume <- NA
    } else {
  
      
    ## Estimating hypervolume based on the environmental data:
    hv <- hypervolume(taxon_clim_data, method = "box")

    ## Extracting volume value:
    volume <- get_volume(hv)
    
  } ## Close if statement.
  
  ## Create dataframe with outputs:
  hv_output <- data.frame(taxon = taxon, 
                          range_size = shape_used$Area[shape_used$Binomial == taxon], 
                          niche_breadth = volume,
                          n_pts = nrow(taxon_clim_data))
			   
  ## Append results to file out of R:
  write.table(hv_output, col.names = FALSE, row.names = FALSE, append = TRUE, eol = "\n", sep = ",",
            file = here("outputs/range_size_climatic_niche_breadth_Australian_squamates_.csv"))
    
  ## Inform progress:
  cat(red(paste0("\n", "Niche breadth for ", taxon, ": ", hv_output$niche_breadth, "\n")))
  
  ## Function will return:
  return(hv_output)
  
} ## End of funtion.

## Running:
purrr::map(unique(clim_data$taxon), get_hypervolume)

## PART 6: Testing significance ----

## Load range size and climatic niche breadth data:
AUS_df <- read.csv(file = here("outputs/range_size_climatic_niche_breadth_Australian_squamates.csv"), header = FALSE) ## If reloading data, may need header = TRUE.
names(AUS_df) <- c("taxon", "range_size", "niche_breadth", "n_pts") ## Add colnames.

## Load biome data and merge:
biome_df <- read.csv(header = TRUE, file = here("outputs/Australian_squamate_taxa_biomes.csv"))
AUS_df <- merge(AUS_df, biome_df, by = "taxon")

## Log-transform:
AUS_df$log_nb <- log10(AUS_df$niche_breadth)
AUS_df$log_rs <- log10(AUS_df$range_size)

## Histograms:
hist(AUS_df$log_nb)
hist(AUS_df$log_rs)

## Scatter plot:
plot(log_rs ~ log_nb, data = AUS_df)

## Linear regression:
lm_res <- lm(log_rs ~ log_nb, data = AUS_df) ; summary(lm_res)

## PART 7: Fig. 1A ----

## Ctenotus taxa to highlight:
highlight <- c("Ctenotus dux", ## Quadrant I.
               "Ctenotus pantherinus", ## Quadrant II.
               "Ctenotus storri", ## Quadrant III.
               "Ctenotus australis") ## Quadrant VI.
  
## Selecting four Ctenotus taxa to highlight in plots:
four_taxa <- AUS_df[AUS_df$taxon %in% highlight, ]
four_taxa$taxon <- factor(four_taxa$taxon, levels = highlight)
  
## Setting labels of highlighted taxa and their position in the plot:
label_df <- data.frame(taxon = highlight,
                       label = c("I", "II", "III", "IV"),
                       nudge_x = c(-1, -1, -1, 1.2), 
                       nudge_y = c(0.15, 0.15, 0.15, -0.15))
  
## Merge label and range size/niche breadth information:
highlight_df <- merge(four_taxa, label_df, by = "taxon")

## Color palette:
plot_palette <- scico(n = max(AUS_df$n_biomes, na.rm = TRUE), palette = "tokyo", direction = -1, begin = 0, end = 1)

## Plotting with ggplot:
Fig1A <- ggplot() +
    
  ## Plotting dots: 
  geom_point(data = AUS_df, aes(x = log_nb, y = log_rs, fill = as.factor(n_biomes)), size = 2.75, color = "black", stroke = 0.80, shape = 21, alpha = 0.90) +
  
  ## Adding regression line (based on all biomes):
  geom_smooth(data = AUS_df, aes(x = log_nb, y = log_rs), method = "lm", se = FALSE, size = 1, color = "gray40") +
  
  ## Add lines corresponding to mean values:
  geom_vline(data = AUS_df, aes(xintercept = mean(log_nb)), color = "gray30", size = 0.5, linetype = "dashed", show.legend = FALSE) +
  geom_hline(data = AUS_df, aes(yintercept = mean(log_rs)), color = "gray30", size = 0.5, linetype = "dashed", show.legend = FALSE) +
  
  ## Plotting colored dots highlighted taxa:
  geom_point(data = four_taxa, aes(x = log_nb, y = log_rs, fill = as.factor(n_biomes)), color = "black", size = 5, shape = 21, stroke = 1.5, show.legend = FALSE) +
  
  ## Add taxon labels:
  geom_text_repel(data = highlight_df, aes(x = log_nb, y = log_rs, label = label), size = 7, segment.size = 0.5, color = "black", direction = "y",
                  nudge_x = highlight_df$nudge_x, nudge_y = highlight_df$nudge_y, show.legend = FALSE) +
  
  ## Editing legends:
  guides(fill = guide_legend(title = "No. biomes")) +
  scale_fill_manual(guide = "legend", values = plot_palette) +

  ## Editing legends and axis labels:
  labs(x = "Climatic niche breadth (log-transformed volume)", 
       y = "Range size (log-transformed area)") +
  
  ## Also setting values on axes:
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  
  ## Some edits: 
  theme_bw() +
  theme(
    
    ## Margin around plot:
    plot.margin = margin(t = 0, r = 10, b = 30, l = 0),
    plot.title = element_blank(),
    axis.title = element_text(size = 18), #margin = margin(t = 5, l = 5, b = 5, r = 5)),
    axis.text = element_text(size = 16), #margin = margin(t = 5, l = 5, b = 5, r = 5)),
    axis.ticks = element_line(size = 0.75, color = "gray20"),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1, colour = "gray20"), 
    
    ## Legend:
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = c(0.9, 0.182),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), ## Getting rid of boxes around legend elements.
    legend.key.size = unit(0.5, 'lines'), ## Distance between legend elements.
    legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80")) ## Box around legend.

## Saving plot (with cowplot):
#ggsave(file = here("plots/Fig1A.jpeg"), plot = Fig1A, width = 8.5, height = 7.5, dpi = 100)

## PART 8: Fig. 1B ----

## Import biome data (using the WWF layer):
biomes <- read_sf(paste0(path, "GIS/terrestrial_ecoregions_WWF/wwf_terr_ecos.shp"))

## Selecting only biomes found in Australia:
biomes_AUS <- filter(biomes, biomes$BIOME %in% c(1, 2, 4, 7, 8, 10, 12, 13))

## Removing to free up memory:
rm(biomes)

## Variables of interest:
biomes_AUS <- biomes_AUS[c("BIOME", "geometry")]

## Replacing numbers by biome names:
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "10", replacement = "Montane grasslands and shrublands")
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "12", replacement = "Mediterranean forests and shrublands")
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "13", replacement = "Deserts and xeric shrublands")
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "1", replacement = "Tropical and subtropical forests")
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "2", replacement = "Tropical and subtropical forests") ## Actually Deciduous forest in Indonesia to the north.
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "4", replacement = "Temperate forests")
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "7", replacement = "Tropical grasslands and shrublands")
biomes_AUS$BIOME <- gsub(x = biomes_AUS$BIOME, pattern = "8", replacement = "Temperate grasslands and shrublands")

## Creating a vector to define order of colors in palette:
biome_info <- data.frame(biome = c("Tropical and subtropical forests", 
                                   "Montane grasslands and shrublands",
                                   "Mediterranean forests and shrublands", 
                                   "Temperate forests", 
                                   "Tropical grasslands and shrublands", 
                                   "Temperate grasslands and shrublands", 
                                   "Deserts and xeric shrublands"),
                          color = c("darkolivegreen","#402716","#883E3A","#CD524C","#E48551","#EDB555","#FFFFCC"))

## Force the order of elements in our plots (maps).
biomes_AUS$BIOME <- factor(biomes_AUS$BIOME, levels = biome_info$biome)

## Biome color palette:
biome_palette <- as.character(biome_info$color)

## Loading world map:
world_map <- ne_countries(scale = "medium", returnclass = "sf")

## Restricting geographic area to map:
bbox_a <- st_bbox(world_map)
bbox_a[["xmin"]] = 113
bbox_a[["xmax"]] = 154
bbox_a[["ymax"]] = -10
bbox_a[["ymin"]] = -45 ## Keeping Tasmania.

## Plot map:
Fig1B <- ggplot() +
  
  ## Adding biome map:
  geom_sf(data = biomes_AUS, aes(fill = BIOME), alpha = 0.75, color = "transparent", show.legend = TRUE, size = 0) +
  
  ## Adding baseline map:
  geom_sf(data = world_map, fill = "transparent", color = "gray30") +
  
  ## Set map boundaries:
  coord_sf(xlim = bbox_a[c("xmin", "xmax")], ylim = bbox_a[c("ymin", "ymax")]) +
  
  ## Labels:
  labs(x = "", y = "") +
  
  ## Colors, legends:
  guides(fill = guide_legend(title = "Biomes")) +
  scale_fill_manual(guide = "legend", values = biome_palette) +
  
  ## Setting a bunch of aesthetic parameters:
  theme_bw() +
  theme(
        ## Legend:
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = c(0.28, 0.15),
        legend.key = element_rect(fill = "transparent", colour = "transparent"), ## Getting rid of boxes around legend elements.
        legend.key.size = unit(1, 'lines'), ## Distance between legend elements.
        legend.background = element_rect(fill = "gray95", size = 0.75, linetype = "solid", colour = "gray80"),  ## Box around legend.
        
        ## Other parameters:
        plot.margin = margin(t = 0, r = 0, b = 20, l = 0),
        plot.title = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(size = 1, colour = "gray20")) 

## Save:
#ggsave(file = here("plots/Fig1B.jpeg"), plot = Fig1B, width = 8, height = 7.5, dpi = 100, unit = "in")

## PART 9: Fig. 1C ----

## Read shape file with taxon distribution polygons:
shape_all <- read_sf(dsn = paste0(path, "GIS/Meiri_etal_ranges/GARD1.1_dissolved_ranges"), layer = "modeled_reptiles")

## Extracting shapefile for the four selected Ctenotus taxa:
shape_highlight <- shape_all[shape_all$Binomial %in% highlight, ]

## Save:
st_write(obj = shape_highlight, dsn = paste0(path, "GIS/Meiri_etal_ranges/selected_taxa_Fig1.shp"))

## Remove shape_all to release memory:
rm(shape_all)
 
## Load:
shape_highlight <- read_sf(dsn = paste0(path, "GIS/Meiri_etal_ranges"), layer = "selected_taxa_Fig1")

## Setting levels of Binomial to fix order of maps in plots (otherwise, it would be alphabetic):
shape_highlight$Binomial <- factor(shape_highlight$Binomial, levels = highlight)

## Numbers:
site <- data.frame(lat = -14, lon = 117, ## Position of labels.
                   text_i = c("I", "II", "III", "IV"),
                   Binomial = highlight)

## Box around map:
bbox_b <- st_bbox(world_map)
bbox_b[["xmin"]] = 114
bbox_b[["xmax"]] = 153
bbox_b[["ymax"]] = -11
bbox_b[["ymin"]] = -39 ## Exclude Tasmania.

## Plot map with ggplot:
Fig1C <- ggplot() +
  
  ## Adding biome map:
  geom_sf(data = biomes_AUS, aes(fill = BIOME), color = "transparent", show.legend = FALSE, size = 0, alpha = 0.8) +
  scale_fill_grey(start = 0, end = 1) +
  
  ## Adding baseline map:
  geom_sf(data = world_map, fill = "transparent", color = "gray30", size = 0.5) +
  
  ## Adding taxon range shapes:
  geom_sf(data = shape_highlight, fill = "blue", color = "blue", size = 0.2, alpha = 0.3) +
  
  ## Set map boundaries to Australia:
  coord_sf(xlim = bbox_b[c("xmin", "xmax")], ylim = bbox_b[c("ymin", "ymax")]) +
  
  ## Numbers:
  geom_text(data = site, aes(y = lat, x = lon, label = text_i), size = 8) +
  
  ## Set facets by taxon:
  facet_wrap(~Binomial, ncol = 4, strip.position = "bottom") +
  
  ## Setting a bunch of aesthetic parameters:
  theme_void() +
  guides(fill = "none") +
  labs(x = "", y = "") +
  
  ## Setting some other elements:
  theme(plot.margin = margin(t = 20, r = 0, b = 0, l = 0),
        panel.background = element_rect(fill = "aliceblue"),
        plot.title = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray30", size = 1),
        strip.text = element_text(size = 18, face = "italic", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        strip.background = element_blank())

## Saving plot: 
#ggsave(plot = Fig1C, width = 15, height = 4, filename = here("plots/Fig1C.jpeg"), dpi = 300, unit = "in")

## PART 10: Combine plots in Fig. 1 ----

## Combining A and B:
AB <- ( Fig1A + Fig1B ) + plot_layout(ncol = 2, widths = c(1.1, 1))
ABC <- AB / Fig1C + plot_layout(nrow = 2, heights = c(1, 0.46))
Fig1 <- ABC + plot_annotation(tag_levels = "A") &
              theme(plot.tag = element_text(size = 24),
                    plot.tag.position = c(0, 0.97),
                    plot.margin = margin(l = 5, t = 5, r = 5, b = 5))

## Save:
ggsave(plot = Fig1, dpi = 100, width = 15, height = 10, units = "in", filename = here("plots/Fig1.jpeg"))
#ggsave(plot = Fig1, dpi = 300, width = 15, height = 10, units = "in", filename = here("plots/Fig1.pdf"))

## End of script.
