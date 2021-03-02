# TITLE: ESU habitat use and crop land cover intersect
# AUTHOR: Braeden Van Deynze
# DATE: Aug., 2020
# INPUTS: Recovery domain files (huc6_esu-withadds-ver5-alb27.shp)
# OUTPUTS: Stats for land cover by species, habitat use

# Prepare environment ----
# Clear environment
rm(list = ls())

# Load libraries
library(CropScapeR)
library(raster)
library(sf)
library(tidyverse)
library(janitor)
library(nhdplusTools)
library(here)




# Load data ----
# Set working directory
# Load recovery domain shape files
(sf_recoverydomain <- st_read("./data/ESU_hucdata/huc6_esu-withadds-ver5-alb27.shp") %>% clean_names()) %>% names

# Key variables are, with * as species code
# r_*: ESU name
# sp_*: ESU spawn/rear
# re_*: ESU rear/migrate
# mi_*: ESU migrate only

# Species codes are (first code is for r_*, second is for sp_*, re_*, mi_*)
# ch/l_chin: chinook
# ch_sp/l_chsp: spring chinook
# ch_sp_su/l_chsps: spring/summer chinook
# ch_fa/l_chfa: fall chinook
# ch_wi/l_chwi: winter chinook
# coho/l_coho: coho
# sthd/l_sthd: steelhead
# chum/l_chup: chum
# chum_su/l_chum: summer chum
# ch_su_fa: summer/fall chinook (no use variables)
# sock/l_sock: sockeye
# pink_*/l_pink: pink (oy/ey: odd year/even year; use data shared)

# Load CDL key
data("linkdata")
linkdata

# Test maps ----

# Idea here is we start with coho, try to find a smaller region so the CDL work
# doesn't take so long

# Generates list of species (switch out "r_coho" in the tabyl call to select another species)
species <- sf_recoverydomain %>%
  tabyl(r_coho) %>% pull(1)
species
i = 2
species[i]

tibble("hab_use" = sf_recoverydomain %>%
         st_drop_geometry() %>%
         filter(r_coho == species[i]) %>%
         select(sp_l_coho:mi_l_coho) %>%
         as.matrix() %*% 1:3) %>% tabyl(hab_use)
# Note that several species do not have 

sf_recoverydomain %>%
  filter(r_coho == species[i]) %>%
  select(c(objectid:states), c(regime01b:r_ch03_c), geometry) %>%
  bind_cols(
    tibble("hab_use" = sf_recoverydomain %>%
      st_drop_geometry() %>%
      filter(r_coho == species[i]) %>%
      select(sp_l_coho:mi_l_coho) %>%
      as.matrix() %*% 1:3)
  ) %>%
  mutate(
    hab_use = case_when(hab_use == 0 ~ NA_real_, TRUE ~ hab_use),
    hab_use = factor(hab_use, levels = c(1, 2, 3), labels = c("spawn+rearing", "rearing+migration", "migration")),
  ) %>%
  ggplot() +
  geom_sf(aes(fill = hab_use)) +
  scale_fill_brewer(type = "seq", direction = -1, drop = FALSE) +
  ggtitle(species[i])

# Test CDL pull ----

# Test polygon
(
  sf_recoverydomain_species <-
    sf_recoverydomain %>%
    filter(r_coho == species[i]) %>%
    select(-everything()) %>%
    summarise()
) %>%
  ggplot() +
  geom_sf()
projection(sf_recoverydomain_species)

# Test raster pull (can be really slow!)
# Full box
(
  cdl_test <- GetCDLData(
    aoi =
      sf_recoverydomain %>%
      filter(r_coho == species[i]) %>%
      select(-everything()) %>%
      summarise(), 
    year = "2018", 
    type = "b",
    format = "table"
  )
) %>% plot()

cdl_test_mask %>% raster::resample()
  ggplot() + geom_raster(aes(x = x, y = y, fill = factor(value)))

projection(cdl_test)
# Different projections

# Mask box
(
  cdl_test_mask <-
    sf_recoverydomain_species %>%
    st_transform(., projection(cdl_test)) %>%
    mask(cdl_test, .)
) %>%
  ggplot() + geom_raster()

# Test post mask stats
(
  cdl_test_stats <-
    freq(cdl_test_mask) %>% 
    #--- matrix to data.frame ---#
    data.frame(.) %>% 
    #--- find share ---#
    mutate(share = count/sum(count)) %>%
    left_join(linkdata, by = c('value' = 'MasterCat'))
)
