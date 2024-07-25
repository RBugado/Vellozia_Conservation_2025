################## DOWNLOAD AND CLEAN DATA FROM GBIF ###########

install.packages("rgbif")

# fill in your gbif.org credentials 
user <- "" # your gbif.org username 
pwd <- "" # your gbif.org password
email <- "" # your email

library(dplyr)
library(readr)  
library(rgbif) # for occ_download

myspecies = c("Vellozia") # we applied the same for Nanuza plicata

# match the names 
gbif_taxon_keys <- name_backbone_checklist(name_data = myspecies, verbose = TRUE)

key <- gbif_taxon_keys$usageKey

# gbif_taxon_keys should be a long vector like this 
#c(2977832,2977901,2977966,2977835,2977863)
# !!very important here to use pred_in!!

occ_download(
  pred_in("taxonKey", key),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

d <- occ_download_get('0067610-240506114902167') %>%
  occ_download_import()

# follow the instructions

write.csv(d, file = "gbif_Vellozia.csv")

################# Cleaning data ####################

library(tidyverse)
library(rgbif)
library(sp)
library(countrycode)
library(CoordinateCleaner)
library(rnaturalearthdata)

d <- read.csv("gbif_Vellozia.csv") # if you start from the beginning d will already be an object

names(d)

## we selected the columns that we thought would be the most useful from the gbif columns

dat <- d %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, eventDate, recordNumber, recordedBy,
         identifiedBy, catalogNumber, occurrenceID, locality)%>% # you might find other ones useful depending on your downstream analyses
  mutate(countryCode = countrycode(d$countryCode, origin =  'iso2c', destination = 'iso3c'))

#Visualize the coordinates on a map

world.inp <- map_data("world")

ggplot() + geom_map(data = world.inp, map = world.inp, aes(x = long, y = lat, 
  map_id = region), 
  fill = "grey80") + xlim(min(dat$decimalLongitude, na.rm = T), 
   max(dat$decimalLongitude, na.rm = T)) + ylim(min(dat$decimalLatitude, na.rm = T),
   max(dat$decimalLatitude, na.rm = T)) + geom_point(data = dat, aes(x = decimalLongitude,                                                                                                                                                                                                            y = decimalLatitude), size = 1) + coord_fixed() + theme_bw() + 
  theme(axis.title = element_blank())

# remove records without coordinates
dat_cl <- dat %>% filter(!is.na(decimalLongitude)) %>% filter(!is.na(decimalLatitude))

# remove records with low coordinate precision
hist(dat_cl$coordinateUncertaintyInMeters/1000, breaks = 30)

dat_cl <- dat_cl %>% filter(coordinateUncertaintyInMeters/1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# remove unsuitable data sources, especially fossils
table(dat$basisOfRecord)
table(dat_cl$basisOfRecord)

dat_cl <- filter(dat_cl, basisOfRecord == "PRESERVED_SPECIMEN" | is.na(basisOfRecord))

# Individual count
table(dat_cl$individualCount)

dat_cl <- dat_cl %>% filter(individualCount > 0 | is.na(individualCount)) %>% 
  filter(individualCount < 99 | is.na(individualCount))  # high counts are not a problem

# Age of records
table(dat_cl$year)

dat_cl <- dat_cl %>% filter(year > 1945)  # remove records from before second world war

table(dat_cl$family)  #that looks good


table(dat_cl$taxonRank)  # We will only include records identified to species level
dat_cl <- dat_cl %>% filter(taxonRank == "SPECIES" | is.na(taxonRank))

# flag problems
dat_cl <- data.frame(dat_cl)
flags <- clean_coordinates(x = dat_cl, lon = "decimalLongitude", lat = "decimalLatitude", 
countries = "countryCode", species = "species", tests = c("capitals", "centroids", 
"equal", "gbif", "zeros", "countries", "seas"))  # most test are on by default

plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

dat_cl <- dat_cl[flags$.summary, ]

world.inp <- map_data("world")

ggplot() + geom_map(data = world.inp, map = world.inp, aes(x = long, y = lat, 
map_id = region), fill = "grey80") + xlim(min(dat$decimalLongitude, na.rm = T), 
max(dat$decimalLongitude, na.rm = T)) + ylim(min(dat$decimalLatitude, na.rm = T), 
max(dat$decimalLatitude, na.rm = T)) + geom_point(data = dat, aes(x = decimalLongitude, 
y = decimalLatitude), colour = "darkred", size = 1) + geom_point(data = dat_cl, 
aes(x = decimalLongitude, y = decimalLatitude), colour = "darkgreen", size = 1) + 
coord_fixed() + theme_bw() + theme(axis.title = element_blank())


ggplot() + geom_map(data = world.inp, map = world.inp, aes(x = long, y = lat, 
map_id = region), fill = "grey80") + xlim(min(dat$decimalLongitude, na.rm = T), 
max(dat$decimalLongitude, na.rm = T)) + ylim(min(dat$decimalLatitude, na.rm = T), 
max(dat$decimalLatitude, na.rm = T)) + geom_point(data = dat_cl, aes(x = decimalLongitude, 
y = decimalLatitude), size = 1) + coord_fixed() + theme_bw() + 
theme(axis.title = element_blank())

write_csv(dat_cl, "Vellozia.record.clean.csv")

dat_cl$species
table(dat_cl$species)

data.frame(table(dat_cl$species))->speciesNum
write_csv(speciesNum, "SpeciesbyNumberVellozia.csv")
