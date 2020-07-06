# TITLE: Recovery domain and crop acreage intersect
# AUTHOR: Braeden Van Deynze
# DATE: June, 2020
# INPUTS: Recovery domain files (subdomains-ver7.shp)
# OUTPUTS: Stats for land cover by recovery domain (and species) for most recent years (and last five years?) w/ map

# Prepare environment ====
# Clear environment
rm(list = ls())

# Load libraries
library(cdlTools)
library(raster)
library(tidyverse)
library(janitor)
library(sf)
library(stars)
library(ggspatial)
library(tmap)
library(units)
library(tabularaster)
library(osmdata)


# Build functions
# Splitting function based on discussion at https://stackoverflow.com/a/41618326
# (From above here) = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# The function spatially aggregates the original raster it turns each aggregated
# cell into a polygon then the extent of each polygon is used to crop the
# original raster. The function returns a list with all the pieces in case you
# want to keep them in the memory. it saves and plots each piece
# The arguments are:
# raster = raster to be chopped            (raster object)
# ppside = pieces per side                 (integer)
# save   = write raster                    (TRUE or FALSE)
# plot   = do you want to plot the output? (TRUE or FALSE)
split_raster <- function(raster,ppside=2,save=F,plot=F){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- raster::aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- extent(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  if(save==T){
    for(i in 1:length(r_list)){
      writeRaster(r_list[[i]],filename=paste("SplitRas",i,sep=""),
                  format="GTiff",datatype="FLT4S",overwrite=TRUE)  
    }
  }
  if(plot==T){
    par(mfrow=c(ppside,ppside))
    for(i in 1:length(r_list)){
      plot(r_list[[i]],axes=F,legend=F,bty="n",box=FALSE)  
    }
  }
  return(r_list)
}
# = = = = = = = = = = = = = = = = = = = = = = = = =

# Compute totals for polygons
count_cdl_pixels <-
  function(sf, cdl){
    for (i in 1:nrow(sf)) {
      # For testing
      # sf = sf_recoverydomain
      # cdl = raster_cdl_merge
      # i = 2
      # End testing block
      
      # Select polygon and reduce size for memory management
      border <- sf[i,]
      clip1 <- raster::crop(cdl, extent(border)) # Clip cdl to rectangle extents of the polygon
      clip2 <- mask(clip1, border) # Mask cdl to only what's within polygon
      clip2_split <- split_raster(clip2) # Final here is a list of four quadrants of the original polygon
      
      # Calculate pixel counts and percentages by class over each quadrant
      ext_long <- tibble()
      for (j in 1:4) {
        # j=2 # For testing
        cn <- cellnumbers(clip2_split[[j]], border) # Generates dummy tibble with cell ids for the split
        if(nrow(cn) > 0) { # Intended to skip quadrants with no cells
          ext <- # Extracts each individual cell from the raster then summarizes
            cn %>% 
            mutate(v = raster::extract(clip2_split[[j]], cell_)) %>% 
            group_by(object_, v) %>% 
            summarize(v_count = n()) %>% 
            mutate(v_pct = v_count / sum(v_count), v = updateNamesCDL(v)) %>%
            arrange(-v_count)
          ext_long <- # Merges it all together
            ext_long %>%
            bind_rows(ext)
        }
      }
      ext_summ <- # Summs across the four quadrants
        ext_long %>%
        group_by(v) %>%
        summarize(v_count = sum(v_count)) %>%
        mutate(v_pct = v_count / sum(v_count)) %>%
        arrange(-v_count)
      ext_wide <- # Projects wide for joining with sf
        ext_summ %>%
        mutate(v = make_clean_names(v)) %>%
        pivot_wider(names_from = v, values_from = c(v_count, v_pct))
      
      # Output to new sf
      border <-
        border %>%
        bind_cols(ext_wide)
      if(i>1) {
        new_sf <-
          new_sf %>%
          bind_rows(border)
      }
      if(i==1) {
        new_sf <-
          border
      }
      
      # For testing
      # border$area
      # plot(clip1)
      # plot(border$geometry, add = TRUE, col = "transparent", border = "red")
      # plot(clip2)
      # End testing block
    }
    
    return(new_sf)
  }


# Set working directory
# wd <- "D:/Braeden/OneDrive/Documents/My Work/NOAA/Agriculture/Analysis/"
# setwd(wd)

# Load recovery domain shape files
# setwd("./data")
sf_recoverydomain <- st_read("./data/recovery_subdomains/subdomains-ver7.shp") %>% clean_names()

# Get CDL data
# raster_cdl <-
#   getCDL(
#     c("CA", "WA", "OR", "ID"),
#     2017
#   )
# NAvalue(raster_cdl$CA2017) <- 0
# NAvalue(raster_cdl$OR2017) <- 0
# NAvalue(raster_cdl$WA2017) <- 0
# NAvalue(raster_cdl$ID2017) <- 0
# raster_cdl_merge <- raster::merge(raster_cdl$CA2017, raster_cdl$OR2017, raster_cdl$WA2017, raster_cdl$ID2017, overlap = TRUE)
# plot(raster_cdl_merge)
# writeRaster(raster_cdl_merge, "cdl_west.tif", format = "GTiff", overwrite = TRUE) # Saves data for downstream use
raster_cdl_merge <- raster("./data/cdl_west.tif")

# Prep maps ====
sf_recoverydomain <-
  sf_recoverydomain %>%
  distinct(.keep_all = TRUE) %>%
  filter(class == "Accessible") %>%
  select(-class, -sourcethm) %>%
  mutate(
    species2 = case_when(
      species == "Steelhead" ~ "Steelhead",
      str_detect(species, "Chinook") ~ "Chinook",
      str_detect(species, "Chum") ~ "Chum",
      str_detect(species, "Pink") ~ "Pink",
      TRUE ~ species
    ),
    species3 = case_when(
      species == "Steelhead" ~ "Trout",
      TRUE ~ "Salmon"
    ),
    domain = case_when(
      subdomain == "Washington Coast" ~ "Puget Sound",
      subdomain == "Snake River" ~ "Interior Columbia",
      subdomain == "Middle Columbia River" ~ "Interior Columbia",
      subdomain == "Upper Columbia River" ~ "Interior Columbia",
      subdomain == "Lower Columbia River" ~ "Willamette/Lower Columbia",
      subdomain == "Upper Willamette River" ~ "Willamette/Lower Columbia",
      subdomain == "California Central Valley" ~ "Central Valley",
      subdomain == "North-Central California Coast" ~ "North-Central California Coast",
      TRUE ~ subdomain
    ),
    area = st_area(geometry),
    area = set_units(area, km^2)
  )

df_rd_species <-
  sf_recoverydomain %>%
  st_drop_geometry() %>%
  mutate(
    area = set_units(area, km^2),
    area_endangered = area * as.numeric(I(status == "Endangered")),
    area_threatened = area * as.numeric(I(status == "Threatened")),
    area_concern = area * as.numeric(I(status == "Species of Concern"))
  ) %>%
  group_by(species) %>%
  summarize(
    n_domain = n_distinct(subdomain),
    n_endangered = sum(status == "Endangered"),
    n_threatened = sum(status == "Threatened"),
    n_concern = sum(status == "Species of Concern"),
    area_tot = sum(area),
    area_endangered = sum(area_endangered),
    area_threatened = sum(area_threatened),
    area_concern = sum(area_concern)
  )
df_rd_species

df_rd_species2 <-
  sf_recoverydomain %>%
  st_drop_geometry() %>%
  mutate(
    area = set_units(area, km^2),
    area_endangered = area * as.numeric(I(status == "Endangered")),
    area_threatened = area * as.numeric(I(status == "Threatened")),
    area_concern = area * as.numeric(I(status == "Species of Concern"))
  ) %>%
  group_by(species2) %>%
  summarize(
    n_domain = n_distinct(subdomain),
    n_endangered = sum(status == "Endangered"),
    n_threatened = sum(status == "Threatened"),
    n_concern = sum(status == "Species of Concern"),
    area_tot = sum(area),
    area_endangered = sum(area_endangered),
    area_threatened = sum(area_threatened),
    area_concern = sum(area_concern)
  )
df_rd_species2

# sf_recoverydomain_merge <-
#   sf_recoverydomain %>%
#   group_by(subdomain) %>%
#   summarize()
# 
# sf_recoverydomain_merge <- st_union_by(sf_recoverydomain$geometry, sf_recoverydomain$subdomain)

# Compute acreage summaries using CDL
# Set CRS to match
sf_recoverydomain <-
  sf_recoverydomain %>%
  st_transform(raster_cdl_merge@crs)

# Break up full count for memory management
# # Full count 1: 1 - 10
# df_crdcount_1 <-
#   count_cdl_pixels(
#     sf = sf_recoverydomain %>% arrange(area) %>% slice(1:10),
#     cdl = raster_cdl_merge
#   ) %>%
#   st_drop_geometry()
# # Full count 2: 11 - 20
# df_crdcount_2 <-
#   count_cdl_pixels(
#     sf = sf_recoverydomain %>% arrange(area) %>% slice(11:20),
#     cdl = raster_cdl_merge
#   ) %>%
#   st_drop_geometry()
# # Full count 3: 21 - 30
# df_crdcount_3 <-
#   count_cdl_pixels(
#     sf = sf_recoverydomain %>% arrange(area) %>% slice(21:30),
#     cdl = raster_cdl_merge
#   ) %>%
#   st_drop_geometry()
# # Full count 3: 31 - 40
# df_crdcount_4 <-
#   count_cdl_pixels(
#     sf = sf_recoverydomain %>% arrange(area) %>% slice(31:40),
#     cdl = raster_cdl_merge
#   ) %>%
#   st_drop_geometry()
# 
# # Combine
# df_crdcount <-
#   bind_rows(df_crdcount_1, df_crdcount_2, df_crdcount_3, df_crdcount_4)
# rm(df_crdcount_1, df_crdcount_2, df_crdcount_3, df_crdcount_4)
# 
# # Save
# write_csv(df_crdcount, "data_subdomains_crd.csv")
df_crdcount <- read_csv("data_subdomains_crd.csv")

# Add collapsed categories
df_crdcount_summ <-
  df_crdcount %>%
  rowwise() %>%
  mutate(
    v_summ_count_total = sum(c_across(contains("_count_")), na.rm = TRUE),
    v_summ_count_developed = sum(c_across(contains("developed") &
                                            contains("_count_")), na.rm = TRUE),
    v_summ_pct_developed = sum(c_across(contains("developed") &
                                          contains("_pct_")), na.rm = TRUE),
    v_summ_count_forest = sum(c_across(contains("forest") &
                                         contains("_count_")), na.rm = TRUE),
    v_summ_pct_forest = sum(c_across(contains("forest") &
                                       contains("_pct_")), na.rm = TRUE),
    v_summ_count_smgrains = sum(c_across((
      contains("rye") |
        contains("barley") |
        contains("winter_wheat") |
        contains("oats") |
        contains("sprint_wheat") |
        contains("triticale") |
        contains("sorghum") |
        contains("buckwheat") |
        contains("millet") |
        contains("durum_wheat")
    ) & contains("_count_")
    ), na.rm = TRUE),
    v_summ_pct_smgrains = sum(c_across((
      contains("rye") |
        contains("barley") |
        contains("winter_wheat") |
        contains("oats") |
        contains("sprint_wheat") |
        contains("triticale") |
        contains("sorghum") |
        contains("buckwheat") |
        contains("millet") |
        contains("durum_wheat")
    )  & contains("_pct_")
    ), na.rm = TRUE),
    v_summ_count_nuts = sum(c_across((
      contains("walnuts") |
        contains("almonds") |
        contains("pistachios")
    ) & contains("_count_")
    ), na.rm = TRUE),
    v_summ_pct_nuts = sum(c_across((
      contains("walnuts") |
        contains("almonds") |
        contains("pistachios")
    ) & contains("_pct_")
    ), na.rm = TRUE),
    v_summ_count_pasture = sum(c_across(contains("pasture") &
                                          contains("_count_")), na.rm = TRUE),
    v_summ_pct_pasture = sum(c_across(contains("pasture") &
                                        contains("_pct_")), na.rm = TRUE),
    v_summ_count_ag = sum(c_across((
      !contains("wetands") &
        !contains("forest") &
        !contains("clover_wildflowers") &
        !contains("aquaculture") &
        !contains("ice_snow") &
        !contains("background") &
        !contains("barren") &
        !contains("shrubland") &
        !contains("water") &
        !contains("summ")
    ) & contains("_count_")
    ), na.rm = TRUE),
    v_summ_pct_ag = sum(c_across((
      !contains("wetlands") &
        !contains("developed") &
        !contains("forest") &
        !contains("clover_wildflowers") &
        !contains("aquaculture") &
        !contains("ice_snow") &
        !contains("background") &
        !contains("barren") &
        !contains("shrubland") &
        !contains("water") &
        !contains("summ")
    ) & contains("_pct_")
    ), na.rm = TRUE),
    v_summ_count_ag_nopasture = v_summ_count_ag - v_summ_count_pasture,
    v_summ_pct_ag_nopasture = v_summ_pct_ag - v_summ_pct_pasture,
    v_summ_count_other = v_summ_count_total - v_summ_count_ag + v_summ_count_developed + v_summ_count_forest,
    v_summ_pct_other = 1 - v_summ_pct_ag - v_summ_pct_developed - v_summ_pct_forest
  )

# Transform for plotting
df_crdcount_long <-
  df_crdcount %>%
  select(-starts_with("v_pct")) %>%
  pivot_longer(starts_with("v_count"), names_to = "v", names_prefix = "v_count_", values_to = "v_count", values_drop_na = TRUE) %>%
  group_by(esu_dps, domain, subdomain) %>%
  mutate(v_pct = v_count / sum(v_count), v_tot = sum(v_count)) %>%
  ungroup()

df_crdcount_summ_long <-
  df_crdcount_summ %>%
  select(
    esu_dps, status, domain, subdomain, species, species2, species3,
    v_summ_pct_ag_nopasture, v_summ_pct_pasture, v_summ_pct_forest, v_summ_pct_developed, v_summ_pct_other
  ) %>%
  pivot_longer(starts_with("v_summ_pct"), names_to = "v", names_prefix = "v_summ_pct_", values_to = "v_summ_pct", values_drop_na = TRUE) %>%
  mutate(
    v = factor(v, levels = c("pasture", "ag_nopasture", "forest", "developed", "other")),
    esu_dps = ordered(esu_dps, df_crdcount_summ %>% arrange(-v_summ_pct_ag) %>% pull(esu_dps)),
    domain = factor(domain, levels = c( "Puget Sound", "Interior Columbia", "Willamette/Lower Columbia", "Oregon Coast", "Southern Oregon/Northern California Coast", "North-Central California Coast", "Central Valley", "South-Central/Southern California Coast")),
    species2 = factor(species2, levels = c( "Steelhead", "Chinook", "Coho", "Sockeye", "Chum", "Pink"))
  )
df_crdcount_summ_long %>% print(n=200)

# Broad Categories All Species
df_crdcount_summ_long %>%
  filter(status != "Not Warranted") %>%
  ggplot() +
  geom_col(aes(x = esu_dps, y = v_summ_pct, fill = v), position = "stack", width = 0.75) +
  geom_point(aes(x = esu_dps, y = -0.05, color = status), fill = NA, size = 5, shape = "square", position = "identity") +
  geom_text(aes(x = esu_dps, y = 1.25, label = domain), vjust = 0.5, size = 3.5, fontface = "plain") +
  scale_fill_manual(
    limits = c("pasture", "ag_nopasture", "forest", "developed", "other"),
    labels = c("Pasture", "Cropland", "Forest", "Developed", "Other"),
    values = c("darkgoldenrod", "goldenrod", "darkgreen", "darkgrey", "lightgrey"),
    name = "Land Use"
  ) +
  scale_color_manual(
    limits = c("Not Warranted", "Species of Concern", "Threatened", "Endangered"),
    values = c("green", "yellow", "darkorange", "red"),
    name = "Status"
  ) +
  scale_y_continuous(breaks = c(-0.05, 0.5, 1.25), labels = c("Status", "Land Use", "Domain"), limits = c(-0.05, 1.4), position = "right") +
  ggthemes::theme_clean() +
  theme(
    axis.text.y = element_text(hjust = 1, vjust = 0.5),
    panel.spacing = unit(0.1, "lines"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_blank(),
    # panel.background = element_rect(color = "black"),
    strip.background = element_rect(fill = "lightgrey"),
    strip.placement = "outside"
  ) +
  facet_grid(rows = "species2", switch = "y", scales = "free_y", space = "free_y") +
  coord_flip()


# df_crdcount_long %>%
#   ggplot(aes(x = "")) +
#   geom_bar(aes(y = v_pct, fill = v), stat = "identity", width = 1) +
#   facet_wrap(~ esu_dps) +
#   coord_polar("y", start = 0)

# Plot maps ====
# For testing, grab just snake river fall-run chinook
# sf_test <- sf_recoverydomain %>% filter(esu_dps == "Snake River Fall-run Chinook Salmon")
# raster_clip_test <- raster::crop(raster_cdl_merge, extent(sf_test)) %>% mask(sf_test) # Mask cdl to only what's within polygon
# raster_df_clip_test <- raster_clip_test %>% as_tibble(xy = TRUE) %>% mutate(cellvalue = updateNamesCDL(cellvalue) %>% factor)
# 
# 
# ggplot() +
#   geom_raster(data = raster_df_clip_test, aes(x = x, y = y)) +
#   coord_sf
# 
# ggplot() +
#   geom_bar(data = raster_cdl_merge, aes())
# 
# plot(raster_cdl_merge, bgc = NA)
# plot(st_geometry(sf_recoverydomain), col = NA, border = "red", bgc = NA, add = TRUE)
# 
# ggplot() +
#   geom_sf(data = sf_recoverydomain, fill = NA, mapping = aes(color = subdomain, linetype = species)) +
#   coord_sf() +
#   facet_wrap(~ species2)
# 
# ggplot() +
#   geom_raster(data = raster_cdl$CA2017) +
#   coord_sf()

states_core <- st_as_sf(map("state", regions = c("california", "oregon", "washington", "idaho"), plot = FALSE, fill = TRUE)) %>% st_transform(st_crs(sf_recoverydomain))
states_expand <- st_as_sf(map("state", regions = c("california", "oregon", "washington", "idaho", "montana", "wyoming", "arizona", "nevada", "utah", "colorado", "new mexico"), plot = FALSE, fill = TRUE)) %>% st_transform(st_crs(sf_recoverydomain))
land <- st_as_sf(map("world", regions = c("Canada", "Mexico"), plot = FALSE, fill = TRUE)) %>% st_transform(st_crs(sf_recoverydomain))


ggplot() +
  geom_sf(data = land, fill = "antiquewhite1") +
  geom_sf(data = states_expand, fill = "antiquewhite1") +
  geom_sf(data = sf_recoverydomain %>% slice(32), fill = NA, color = "red") +
  coord_sf(
    xlim = st_bbox(states_core)[c(1,3)],
    ylim = st_bbox(states_core)[c(2,4)]
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "aliceblue", color = NA)
  )
