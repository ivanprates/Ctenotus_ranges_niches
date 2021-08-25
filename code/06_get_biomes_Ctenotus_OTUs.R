####################
### The goals of this script are:
### To extract biomes occupied by delimited Ctenotus OTUs.
### Checked for functionality on April 21st 2021.
### Written by Ivan Prates (ivanprates.org).

## PART 1: Getting ready ----

## Loading packages:
require(crayon)
require(raster)
require(sf)
require(tidyverse)

## If needed, clearing working space:
rm(list = ls())

## Working directory:
path <- "~/Dropbox/Science/MYPAPERS_ongoing/spheno_ddRAD/" ## Rhinellax.
setwd(paste0(path, "Dryad"))
detach("package:here", unload = TRUE)
library(here)

## PART 2: Prepare data ----

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

## OTU assigments and sample information:
assignments <- read.csv(file = here("outputs/assignments_Ctenotus_OTUs.csv"), header = TRUE)
sample_info <- read.csv(file = here("sample_information/ddRAD_sample_information.csv"), header = TRUE)
sample_info <- merge(sample_info, assignments, by.x = "SAMPLE_ID", by.y = "ID")

## PART 3: Function: Extract biomes for each OTU ----

## Function: generate random points and extract environmental data:
get_biome <- function(OTU) {
  
  ## Testing:
  #OTU <- "leonhardii_gr_cluster_2"
  
  ## Inform progress:
  cat(blue(paste0("\n", "Now processing ", OTU, ":\n")))
  
  ## Extracting lat-lon for OTU:
  OTU_df <- sample_info[sample_info$cluster == OTU, c("LON", "LAT")]
  
  ## Extract raster to points:
  ## Important to specify the raster package as there are many functions "extract" out there!
  extracted_biome <- raster::extract(biomes, OTU_df)
    
  ## How many biomes the taxon occupies?
  n_biomes <- length(unique(extracted_biome$BIOME))
      
  ## Keep the most frequent biome occupied by species:
  most_freq_biome <- names(which.max(table(extracted_biome$BIOME)))
      
  ## What proportion of total points occupy the most frequent biome?
  biome_proportion <- table(extracted_biome$BIOME == most_freq_biome)/length(extracted_biome$BIOME)
  biome_proportion <- round(biome_proportion["TRUE"], digits = 2)
      
  ## Create dataframe with outputs:
  biome_output <- data.frame(OTU = OTU, 
                             n_biomes = n_biomes,
                             biome = most_freq_biome,
                             biome_prop = biome_proportion)
  
  ## Inform progress:
  cat(red(paste0("\n", "Number of biomes for ", OTU, ": ", biome_output$n_biome, "\n")))
  
  ## Function will return:
  return(biome_output)
  
} ## End of function.
                    
## PART 4: Run for all OTUs ----

## Using purrr:
biome_df <- map_dfr(unique(sample_info$cluster), get_biome)

## Save:
write.csv(biome_df, file = here("outputs/biomes_Ctenotus_OTUs_.csv"), row.names = FALSE)

# End of script.
