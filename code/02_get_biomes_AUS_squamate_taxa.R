####################
### The goals of this script are:
### To extract biomes occupied by Australian squamate taxa.
### Checked for functionality on April 20th 2021.
### Written by Ivan Prates (ivanprates.org).

## PART 1: Getting ready ----

## Packages:
require(crayon)
require(rgdal)
require(hypervolume)
library(magrittr)
require(raster)
require(sf)
require(sp)
require(spThin)
require(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## Creating folders to save results:
dir.create(path = here("outputs"))

## PART 2: Prepare taxon distribution shape files ----

## Read shape file with taxon distribution polygons:
shape_all <- read_sf(dsn = paste0(path, "GIS/Meiri_etal_ranges/GARD1.1_dissolved_ranges"), layer = "modeled_reptiles")
print("All reptile shapes loaded!")

## Read list of Australian taxa to consider:
AUS_taxa <- read.csv(here("outputs/Australian_squamate_taxa_list.csv"), header = TRUE)

## Keep shapes for Australian taxa only:
shape_AUS <- shape_all[shape_all$Binomial %in% AUS_taxa$taxon, ]

## Remove to free up memory:
rm(shape_all)

## Converting to spatial polygon class:
shape_used <- as(shape_AUS, "Spatial")
print("Taxon range shapes converted to 'spatial' class!")
    
## PART 3: Prepare biome data ----

## Read shape file with biomes:
biomes <- read_sf(paste0(path, "GIS/terrestrial_ecoregions_WWF/wwf_terr_ecos.shp"))
print("Biome shapes loaded!")

## Biomes are indicated by numbers, but we want them to be read as characters:
biomes$BIOME %<>% as.character()

## Replacing numbers by biome names:
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "10", replacement = "alp_grassland")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "11", replacement = "tundra")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "12", replacement = "mediter_woodland")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "13", replacement = "desert")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "14", replacement = "mangroves")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "1", replacement = "moist_forest")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "2", replacement = "decid_forest")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "3", replacement = "trop_subtrop_conifer_forest")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "4", replacement = "temp_forest")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "5", replacement = "temp_conifer_forest")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "6", replacement = "boreal_forest_taiga")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "7", replacement = "trop_grassland")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "8", replacement = "temp_grassland")
biomes$BIOME <- gsub(x = biomes$BIOME, pattern = "9", replacement = "flood_grassland")

## Converting to spatial polygon class:
biomes <- as(biomes, "Spatial")
print("Biome shapes converted to 'spatial' class!")

## PART 4: Function: Extract biomes for each taxon ----

## Function:
get_biome <- function(taxon){
  
  ## For testing purposes:
  #taxon = "Ctenotus robustus"

  ## Inform progress:
  cat(blue(paste0("\n", "Now processing ", taxon, ":\n")))
  
  ## Extracting polygon for taxon:
  polygon <- shape_used[shape_used$Binomial == taxon, ]
  
  ## Setting minimum range size (20 Km2):
  if (polygon$Area < 20) { ## Will exclude species from islets or known from the type locality only.
    n_biomes <- NA
    most_freq_biome <- NA
    biome_proportion <- NA
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
    
    ## Extract raster to points. Important to specify the raster package as there are many functions "extract" out there!
    extracted_biome <- raster::extract(biomes, thinned_pts_df)
    
    ## If all biome rows == NA:
    if (is.na(unique(extracted_biome$BIOME)[1])) {
      n_biomes <- NA
      most_freq_biome <- NA
      biome_proportion <- NA
    } else {
    
      ## Remove rows with NA:
      extracted_biome <- extracted_biome[is.na(extracted_biome$BIOME) == FALSE, ]
      
      ## How many biomes the taxon occupies?
      n_biomes <- length(unique(extracted_biome$BIOME))
      
      ## Keep the most frequent biome occupied by taxon:
      most_freq_biome <- names(which.max(table(extracted_biome$BIOME)))
      
      ## What proportion of total points correspond to the most frequent biome?
      biome_proportion <- table(extracted_biome$BIOME == most_freq_biome)/length(extracted_biome$BIOME)
      biome_proportion <- round(biome_proportion["TRUE"], digits = 2)
      
    }} ## Close if statements.
    
  ## Create dataframe with outputs:
  biome_output <- data.frame(taxon = taxon, 
                             n_biomes = n_biomes,
                             main_biome = most_freq_biome,
                             prop_main_biome = biome_proportion)
  
  ## Append results to file out of R:
  write.table(biome_output, col.names = FALSE, row.names = FALSE, append = TRUE, eol = "\n", sep = ",",
              file = here("outputs/Australian_squamate_taxa_biomes.csv"))
  
  ## Inform progress:
  cat(red(paste0("\n", "Biome for ", taxon, ": ", biome_output$biome, "\n")))
  
  ## Function will return:
  return(biome_output)
  
} ## End of function.
                    
## PART 5: Apply function to extract biomes ----

## Apply function across taxa:
biome_df <- map_dfr(.x = shape_used$Binomial, .f = function(x) get_biome(taxon = x))

## End of script.
