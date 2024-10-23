
## Pre-process mammal range map rasters ##
rm(list=ls())

# These range maps come from the PHYLACINE database (you can download the "Data" folder from their online repository)
# The range maps are based on IUCN but were checked for the development of PHYLACINE. 
# The folder also has the phylogeny that we will be using for PART 2. 
# For simplicity, we will do our analysis on country-scale (.shp file). I will use the TDWG botanical country scale (level 3). 
# It mostly captures countries, but for very diverse countries, the country is divided into smaller regions that are individual botanical countries.


# Libraries ================
library(terra)
library(here)
library(sf)
library(dplyr)

# Set Folder Path ===========
here()
out_path <- here("Diversification", "data")




# ======================================================== #
# Country-scale map
# ======================================================== #

bot_v <- st_read(here("Diversification", "data_raw", "shp")) %>%
  vect()


# ======================================================== #
## Species ranges 
# !!!  Paste the Ranges folder from PHYLACINE into the data_raw folder for this to work !!!
# ======================================================== #


mammals_ranges_paths <- list(
    curr = here("Diversification", "data_raw" , "Data", "Ranges", "Current"),
    past = here("Diversification", "data_raw", "Data", "Ranges","Present_natural"))

# ======================================================== #
# Read ranges into a raster 
# ======================================================== #

for (i in seq_along(names(mammals_ranges_paths))){

tif_files_i <- list.files( 
    file.path(mammals_ranges_paths[[i]]), 
    pattern = "\\.tif$", 
    full.names = TRUE)

r_mammals_i <- rast(tif_files_i)
saveRDS(r_mammals_i, paste0(out_path, "/mammals_rasterstack_", names(mammals_ranges_paths)[i], ".rds"))


}
# ======================================================== #
# ## Read in .tif files with range maps to a raster stack
# ======================================================== #


r_mammals <- list(
  curr = readRDS(here(out_path, "mammals_rasterstack_curr.rds")),
  past = readRDS(here(out_path, "mammals_rasterstack_past.rds"))
)

## They are not the same, so we have to project the botanical countries to the same CRS as the raster
bot_v2 <- terra::project(bot_v, r_mammals[[1]])
crs(r_mammals[[2]]) == crs(bot_v2)
bot_v3 <- bot_v2 %>% st_as_sf()

# ======================================================== #
# Extract occurrence in botanical countries ==============
# ======================================================== #

for (i in seq_along(r_mammals)){
  r_mammals_i <- r_mammals[[i]]
  length(names(r_mammals_i))


  df_values_i <- data.frame()
  # Use lapply and bind the results to the dataframe
  df_values <- do.call(rbind, lapply(seq_along(names(r_mammals_i)), function(tax) {
    sp <- names(r_mammals_i)[tax] # species name
    print(sp)
    # Extract presence/absence from range maps on country scale
    v <- unique(terra::extract(r_mammals_i[[tax]], bot_v2, touches = TRUE))
    v$species <- names(v[2])
    names(v) <- c("LEVEL_3_ID", "presence", "species")
    head(v)
    # Convert 'v' to a dataframe and add the species name
    df <- as.data.frame(v)
    return(df)
     }
    )
)

saveRDS(df_values, paste0(out_path, "/extracted_values_mammals_df_", names(mammals_ranges_paths)[i], ".rds"))

}

# ======================================================== #
df_values <- list(
  curr =  readRDS(here(out_path, "extracted_values_mammals_df_curr.rds")) %>% filter(presence != 0),
  past =  readRDS(here(out_path, "extracted_values_mammals_df_past.rds")) %>% filter(presence != 0))


mammals_occu <- list(
  curr = left_join(bot_v3, df_values[[1]]),
  past = left_join(bot_v3, df_values[[2]])
)

saveRDS(mammals_occu, here(out_path, "mammals_occu.rds"))





## Sanity checks ##
mammals_occu <- readRDS(here(out_path, "mammals_occu.rds"))

setdiff(unique(mammals_occu$curr$species), unique(mammals_occu$past$species))

setdiff(unique(mammals_occu$past$species), unique(mammals_occu$curr$species))   %>% length() # 380 species went extinct


unique(mammals_occu$past$species) %>% length() #5815
unique(mammals_occu$curr$species) %>% length() #5435
