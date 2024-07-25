#### Extinction Risk Assessment with IUCN Criterion B ConR Code ####

#### Required Packages ####

install.packages("devtools")
devtools::install_github("gdauby/ConR")
install.packages("lwgeom")
install.packages("tidyverse")

library(ConR)
library(lwgeom)
library(tidyr)
library(dplyr)

#### Data ####

occs <- read.csv("Vellozia_torquataCoords.csv")
spp <- unique(occs$species)

## Data must be in format: latitude, longitude, species
## with the names "ddlat" "ddlon" and "tax"

## Code to rearrange columns from species, longtiude, latitude to correct order

occs <- occs %>% relocate(species, .after = Longitude)
occs <- occs %>% relocate(Latitude, .before = Longitude)
colnames(occs) <- c("ddlat", "ddlon", "tax")

## Calculate EOO

EOO.hull <- EOO.computing(XY = occs, method.range = "convex.hull",
                          export_shp = TRUE, show_progress = FALSE)

## Calculate AOO

AOO <- AOO.computing(occs,
                     cell_size_AOO = 2, nbe.rep.rast.AOO = 30,
                     show_progress = FALSE)

## Radius, and therefore subpopulation and severe fragmentation calculations
## only work if there are at least 4 occurrence records. If subpopulation
## is not calculated, it defaults to the number of locations.

## Calculate radius for subpopulations

radius <- subpop.radius(XY = occs[, c(1:3)],
                        quant.max = 0.9)

## Calculate number of subpopulations

sub <- subpop.comp(XY = occs,
                   resol_sub_pop = radius[,c("tax", "radius")],
                   show_progress = FALSE)

## Determine if fragmentation is severe

sever.frag <- severe_frag(XY = occs,
                          resol_sub_pop = radius[,c("tax", "radius")],
                          dist_isolated = radius$radius,
                          show_progress = FALSE)

## Calculate number of locations

locs <- locations.comp(occs,
                           method = "fixed_grid",
                           nbe_rep = 30,
                           cell_size_locations = 10,
                           rel_cell_size = 0.05,
                           method_polygons = "no_more_than_one",
                           show_progress = FALSE)

## Determine if there is decline in habitat

rel.loss <- 1 # We set to decreasing for all species based on our SDMs
declineB <- ifelse(rel.loss >= 1, "Decreasing", "Not Decreasing")

#### Testing Criterion B ####

## If subpopulation and severe fragmentation could be calculated,
## use the first one. If it couldn't be calculated, use the second.

critB_bahiana <- criterion_B(x = occs,
                             AOO = AOO,
                             EOO = EOO.hull$results,
                             locations = locs$locations,
                             severe.frag = sever.frag,
                             subpops = sub, 
                             decline = declineB,
                             show_progress = FALSE)

critB_torquata <- criterion_B(x = occs,
                              AOO = AOO,
                              EOO = EOO.hull$results,
                              locations = locs$locations,
                              decline = declineB,
                              show_progress = FALSE)

write.table(critB_torquata, "V_torquata_IUCN.csv", 
            sep = ",",
            row.names = FALSE)

#### Compile the results into a single table ####

combined_species <- bind_rows(combined_species, critB_torquata,
                              .id = NULL)
write.table(combined_species, "Vellozia_IUCN.csv", 
            sep = ",",
            row.names = FALSE)
