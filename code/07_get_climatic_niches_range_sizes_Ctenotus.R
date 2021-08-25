####################
### The goals of this script are:
### To estimate geographic range sizes and climatic niche breadths for the delimited Ctenotus OTUs and corresponding taxa.
### To plot Fig. 4 of the manuscript based on the resulting estimates.
### Checked for functionality on April 20th 2021.
### Written by Ivan Prates (ivanprates.org).

## PART 1: Getting ready ----

## Loading packages:
library(ggimage)
library(hypervolume)
library(pals)
library(patchwork)
library(raster)
library(reshape2)
library(scico)
library(sf)
library(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## Creating folders to save results:
dir.create(path = here("climatic_data_Ctenotus_OTUs"))

## PART 2: Preparing environmental data ----

## List environmental layers:
layers_clim <- list.files(path = paste0(path, "GIS/chelsa_bioclim"),  pattern = '*.tif', full.names = TRUE)

## If other bio variables are included in this folder, keep only bio1, bio4, bio12, and bio15:
#layers_clim <- layers_clim[c(1, 4, 12, 15)]

## Stacking variables:
layers_clim <- raster::stack(layers_clim)

## Importing OTU assignments:
assignments <- read.csv(file = here("outputs/assignments_Ctenotus_OTUs.csv"), header = TRUE)

## Importing lat-lon information for all samples:
sample_info <- read.csv(file = here("sample_information/ddRAD_sample_information.csv"), header = TRUE)

## Keep only samples assigned to an OTU:
sample_info <- sample_info[sample_info$SAMPLE_ID %in% assignments$ID, ]

## Change column names:
names(sample_info)[names(sample_info) == "SAMPLE_ID"] <-"ID"
names(sample_info)[names(sample_info) == "LAT"] <- "Latitude"
names(sample_info)[names(sample_info) == "LON"] <- "Longitude"
names(sample_info)[names(sample_info) == "UPDATED_SP"] <- "taxon"

## PART 3: Extracting environmental data from the sampled sites of individuals ----

## Extract raster to points:
## Important to specify package as there are many functions 'extract' out there!
environ_data <- as.data.frame(raster::extract(layers_clim, sample_info[c("Longitude", "Latitude")]))
  
## Add some data back to extracted data:
environ_data$ID <- sample_info$ID
environ_data$Latitude <- sample_info$Latitude
environ_data$Longitude <- sample_info$Longitude
environ_data$taxon <- sample_info$taxon
  
## Change column names for variables:
names(environ_data)[names(environ_data) == "CHELSA_bio10_01"] <- "Annual_mean_temperature"
names(environ_data)[names(environ_data) == "CHELSA_bio10_04"] <- "Temperature_seasonality"
names(environ_data)[names(environ_data) == "CHELSA_bio10_12"] <- "Annual_precipitation"
names(environ_data)[names(environ_data) == "CHELSA_bio10_15"] <- "Precipitation_seasonality"
  
## Keep only samples with environmental information:
environ_data <- environ_data[!is.na(environ_data[["Annual_mean_temperature"]]), ]
  
## Saving:
write.csv(x = environ_data, file = here("climatic_data_Ctenotus_OTUs/extracted_climatic_data.csv"), row.names = FALSE, quote = FALSE)

## PART 4: Function: Estimating climatic niches ----

## Function:
get_niches <- function(group_var) {
  
  ## Testing:
  #group_var <- "cluster"
  #group_var <- "taxon"
  
  ## Read environ data:
  environ_data <- read.csv(header = TRUE, file = here("climatic_data_Ctenotus_OTUs/extracted_climatic_data.csv"))
  
  ## Merge with sample info:
  environ_data <- merge(environ_data, assignments, by = "ID")
  
  ## Restrict to group_var:
  environ_data <- environ_data[c("ID", "group", group_var, 
                                 "Latitude", "Longitude",
                                 "Annual_mean_temperature", "Temperature_seasonality", 
                                 "Annual_precipitation", "Precipitation_seasonality")]
  
  ## Transform the data:
  for (var in c("Annual_mean_temperature", "Temperature_seasonality", "Annual_precipitation", "Precipitation_seasonality")) {
   environ_data[[var]] <- scale(environ_data[[var]], center = FALSE, scale = TRUE) ## z-scores transformation.
  }
  
  ## Function:
  get_hv <- function(group) {
  
    ## Testing:
    #group <- "atlas_gr_cluster_2"
    
    ## Print status:
    print(paste0("Now processing ", group, "!"))
    
    ## Select environmental data for corresponding group:
    climate_df <- environ_data[environ_data[[group_var]] == group, 
                                     c("Annual_mean_temperature", "Temperature_seasonality", 
                                       "Annual_precipitation", "Precipitation_seasonality")]
     
    ## Test if all values are the same for a variable, that is, if there's no environ variation.
    unique_vals <- function(cl){ length(unique(cl)) > 1 } ## Function: Is there more than one value?
    test_unique <- apply(climate_df, 2, unique_vals) ## Test by applying function.
    
    ## Run next steps only if all variables are not invariable.
    if(length(test_unique[test_unique == "FALSE"]) > 0) { volume <- data.frame() ## Same as adding NA.
      } else {
      
      ## Estimating niche hypervolumes based on the environmental data:
      #hv_out <- hypervolume(climate_df, method = "gaussian")
      hv_out <- hypervolume(climate_df, method = "box") ## Much faster.
      
      ## Extracting volume value:
      volume <- get_volume(hv_out)
  
      ## Return:
      return(data.frame(group = group, niche_breadth = volume))
      
    } ## Close if statement.
      
  } ## End of function.

  ## Apply to all groups:
  hv_df <- map_dfr(unique(environ_data[[group_var]]), get_hv)
  
  ## Save results:
  if (group_var == "cluster") { group_label <- "OTUs" }
  if (group_var == "taxon") { group_label <- "taxa" }
  write.csv(hv_df, file = here(paste0("climatic_data_Ctenotus_OTUs/climatic_niche_breadths_Ctenotus_", group_label, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Function will return:
  return(hv_df)
  
} ## End of function.

## PART 5: Function: Estimating range sizes ----

## Function:
get_ranges <- function(group_var) {
  
  ## Testing:
  #group_var <- "cluster"
  #group_var <- "taxon"
  
  ## Read environ data:
  environ_data <- read.csv(header = TRUE, file = here("climatic_data_Ctenotus_OTUs/extracted_climatic_data.csv"))
  
  ## Merge with assignments:
  environ_data <- merge(environ_data, assignments, by = "ID")
  
  ## Restrict to group_var:
  environ_data <- environ_data[c(group_var, "group",
                                 "Latitude", "Longitude",
                                 "Annual_mean_temperature", "Temperature_seasonality", 
                                 "Annual_precipitation", "Precipitation_seasonality")]
  
  ## Function:
  get_hulls <- function(group) {
    
    ## Testing:
    #group = "atlas_gr_cluster_2"
    
    ## Import shape of Australia to crop polygons, if needed:
    AUS_admin <- read_sf(dsn = paste0(path, "GIS/AUS_adm/", layer = "AUS_adm0.shp"))
  
    ## Print status:
    print(paste0("Now processing ", group, "!"))
    
    ## Keeping only records for group:
    group_info <- environ_data[(environ_data[[group_var]] == group), ]
    
    ## Defining convex polygon based on x and y variables.
    open_hull <- chull(x = group_info$Latitude, y = group_info$Longitude)
    
    ## Closing hull by repeting the first point. Also, keeping only a few columns of the resulting dataframe:
    closed_hull <- group_info[c(open_hull, open_hull[1]), c(group_var, "Longitude", "Latitude")]
    
    ## Creating spatial polygon from hull using the sp package.
    ## This involves a weird sequence of iteratively converting the products of functions into lists:
    spatial_polygon <- closed_hull[, c("Longitude", "Latitude")] %>% 
      Polygon() %>% list() %>% Polygons(ID = 1) %>% list() %>%
      SpatialPolygons(proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    ## Crop spatial polygon using map of Australia:
    if (length(closed_hull$subcluster) >= 3) { ## We can't do it for groups with less than three sites.
      spatial_polygon <- raster::crop(spatial_polygon, AUS_admin) 
    }
    
    ## Calculate area of polygon (in squared meters):
    polygon_area <- raster::area(x = spatial_polygon)
    
    ## Convert to squared kilometers and eliminate decimal digits:
    polygon_area <- round(x = (polygon_area/1000000), digits = 0)
    
    ## Return:
    return(data.frame(group = group, range_size = polygon_area))
    
  } ## End of function.
  
  ## Apply to all groups:
  rs_df <- map_dfr(unique(environ_data[[group_var]]), get_hulls)
  
  ## Replace area for groups == 0, that is, with too few records to estimate area.
  rs_df$range_size[rs_df$range_size == 0] <- NA ## Replace with NA.
  #rs_df$range_size[rs_df$range_size == 0] <- 10 ## Replace with a small number.
  
  ## Save results:
  if (group_var == "cluster") { group_label <- "OTUs" }
  if (group_var == "taxon") { group_label <- "taxa" }
  write.csv(rs_df, file = here(paste0("climatic_data_Ctenotus_OTUs/range_sizes_Ctenotus_", group_label, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Function will return:
  return(rs_df)
  
} ## End of function.

## PART 6: Function: Estimate range sizes and climatic niches and combine them ----

## Function:
get_rs_nb <- function(group_var) {

  ## Testing:
  #group_var <- "cluster"
  #group_var <- "taxon"
  
  ## Run functions:
  nb_df <- get_niches(group_var = group_var)
  rs_df <- get_ranges(group_var = group_var)
  
  ## Combine results:
  dist_df <- merge(nb_df, rs_df, by = "group")
  
  ## Estimate log10:
  dist_df$log_nb <- log10(dist_df$niche_breadth)
  dist_df$log_rs <- log10(dist_df$range_size)
  
  ## Save results:
  if (group_var == "cluster") { group_label <- "OTUs" }
  if (group_var == "taxon") { group_label <- "taxa" }
  write.csv(dist_df, file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_", group_label, ".csv")), row.names = FALSE, quote = FALSE)
  
  ## Return:
  return(dist_df)
  
} ## End of function.

## Use function to estimate range sizes and niche breadths for both Ctenotus OTUs and traditional taxa:
get_rs_nb(group_var = "cluster")
get_rs_nb(group_var = "taxon")

## PART 7: Counting how many OTUs and taxa are being used in analyses ----

## How many unique OTUs?
n_OTUs <- length(unique(assignments$cluster)) ; n_OTUs
n_samples <- length(unique(assignments$ID)) ; n_samples
n_groups <- length(unique(assignments$group)) ; n_groups

## Counting taxa used in OTU analyses:
sample_info <- read.csv(file = here("sample_information/ddRAD_sample_information.csv"), header = TRUE)
info_df <- merge(sample_info, assignments, by.x = "SAMPLE_ID", by.y = "ID")
n_taxa <- length(unique(info_df$UPDATED_SP)) ; n_taxa

## PART 8: Plotting range size vs. climatic niche breadth for Ctenotus OTUs ----

## Read results from previous run for AUS taxa:
AUS_df <- read.csv(file = here("outputs/range_size_climatic_niche_breadth_Australian_squamates.csv"), header = TRUE)

## Estimate log10:
AUS_df$log_nb <- log10(AUS_df$niche_breadth)
AUS_df$log_rs <- log10(AUS_df$range_size)
  
## Range and niche data data:
rs_nb_df <- read.csv(header = TRUE, file = here(paste0("outputs/range_size_climatic_niche_breadth_Ctenotus_OTUs.csv")))

## Biome data:
biome_df <- read.csv(header = TRUE, file = here(paste0("outputs/biomes_Ctenotus_OTUs.csv")))     

## Merge:
rs_nb_df <- merge(biome_df, rs_nb_df, by.x = "OTU", by.y = "group")

## Color palette:
n_colors <- max(rs_nb_df$n_biomes)+1
palette <- scico(n = n_colors, palette = "lajolla", direction = 1, begin = 0, end = 1)
palette <- palette[c(1, 3:n_colors)]

## Plot:
plot <- ggplot() +
    
  ## Plotting data for AUS taxa: 
  geom_point(data = AUS_df, aes(x = log_nb, y = log_rs), fill = "gray70", color = "gray70", size = 2.5, shape = 21, stroke = 0.75) +
  geom_smooth(data = AUS_df, aes(x = log_nb, y = log_rs), method = "lm", se = FALSE, size = 0.5, color = "gray50") +
    
  ## Adding data for OTUs:
  geom_smooth(data = rs_nb_df, aes(x = log_nb, y = log_rs), method = "lm", se = FALSE, size = 0.5, color = "blue") +
  geom_point(data = rs_nb_df, aes(x = log_nb, y = log_rs, fill = as.factor(n_biomes)), size = 2.75, color = "black", stroke = 0.5, shape = 21, alpha = 0.75) +
      
  ## Editing legends:
  guides(fill = guide_legend(title = "No. biomes")) +
  scale_fill_manual(guide = "legend", values = palette) +
  
  ## Editing legends and axis labels:
  labs(x = "Climatic niche breadth (log-transformed volume)", 
       y = "Range size (log-transformed area)") +
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
        legend.text = element_text(size = 10),
        legend.position = c(0.88, 0.15),
        legend.key = element_rect(fill = "transparent", colour = "transparent"), ## Getting rid of boxes around legend elements.
        legend.key.size = unit(0.5, 'lines'), ## Distance between legend elements.
        legend.background = element_rect(fill = "gray90", size = 0.75, linetype = "solid", colour = "gray80")) + ## Box around legend.
                  
  ## Add lizard image on top of plot:
  geom_image(aes(x = -4.5, y = 6.4), image = paste0(path, "figures/IMG_8004_Pascal_edited.jpg"), size = 0.4) +
    
  ## Add label on top of plot:
  geom_text(aes(x = -4, y = 6.05), label = "Ctenotus pantherinus", size = 3, fontface = "italic", color = "gray20")
      
## Density plot for range size:
plot_d_nb <- ggplot() + theme_classic() + #theme_nothing() +
  geom_density(data = AUS_df, aes(x = log_nb), color = "gray40", fill = "gray40", alpha = 0.5, trim = FALSE, size = 0.25) +
  geom_density(data = rs_nb_df, aes(x = log_nb), color = "blue", fill = "blue", alpha = 0.2, trim = FALSE, size = 0.25) +
  labs(y = "Density") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.25, colour = "gray80"),
        axis.line = element_line(size = 0.5, colour = "gray20"),
        axis.ticks = element_line(size = 0.75, color = "gray20"), ## Ticks thickness.
        axis.title = element_text(size = 14), # margin = margin(t = 5, l = 5, b = 5, r = 5)),
        axis.text = element_text(size = 12)) #, margin = margin(t = 5, l = 5, b = 5, r = 5)))
  
## Density plot for range size:
plot_d_rs <- ggplot() + theme_classic() + #theme_nothing() +
  geom_density(data = AUS_df, aes(y = log_rs), color = "gray40", fill = "gray40", alpha = 0.5, trim = FALSE, size = 0.25) +
  geom_density(data = rs_nb_df, aes(y = log_rs), color = "blue", fill = "blue", alpha = 0.2, trim = FALSE, size = 0.25) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Density") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.25, colour = "gray80"),
        axis.line = element_line(size = 0.5, colour = "gray20"),
        axis.ticks = element_line(size = 0.75, color = "gray20"), ## Ticks thickness.
        axis.title = element_text(size = 14), # margin = margin(t = 5, l = 5, b = 5, r = 5)),
        axis.text = element_text(size = 12)) #, margin = margin(t = 5, l = 5, b = 5, r = 5)))
  
## Combine plots:
plot_a <- ( plot_d_nb | plot_spacer() ) + plot_layout(widths = c(3.4, 1))
plot_b <- ( plot | plot_d_rs ) + plot_layout(widths = c(3.5, 1))
Fig <- ( plot_a / plot_b ) + plot_layout(heights = c(1, 3.5))

## Saving plot (with cowplot):
ggsave(file = here("plots/Fig4.jpeg"), plot = Fig, width = 8, height = 7, device = "jpeg", units = "in", dpi = 200)

## End of script.
