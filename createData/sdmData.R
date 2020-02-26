
#### Aim of prog: Prepare the presence-absence data for the 14 species
## Format the presence absence data:
#		- Subset for the 14 species, and for "recent" data (from 1997 to 2010)
#		- Calculate Presence/Absence = TRUE/FALSE
#		- Rename columns (for the random forest, otherwise problems in formula)
#		- Change the coordinates to match climate coordinates

#### Load packages and tool functions; clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

library(data.table)
library(stringi)
library(raster)

#### Tool functions
## To sort species, extract tsn, and create colnames data.table can handle
sortingSpecies = function(ls_species)
{
	tsn_dt = data.table(species = ls_species, codeName = character(length(ls_species)),
		tsn = character(length(ls_species))) # tsn = Taxonomic Serial Number

	tsn_dt[, codeName := stri_sub(str = ls_species,
		from = stri_locate_first(str = ls_species, regex = "-")[,1] + 1)]

	tsn_dt[, tsn := stri_sub(str = ls_species, from = 1,
		to = stri_locate_first(str = ls_species, regex = "-")[,1] - 1)]

	tsn_dt[, reversed := stri_replace_all(paste0(codeName, "_" , tsn), regex = "-", "_")]

	tsn_dt = tsn_dt[order(codeName),]
	return (tsn_dt)
}

## For trees that were both absent and present the past selected years. Set the value to last state
delDichotomy = function(booleanVec)
{
	n = length(booleanVec)
	print(booleanVec)
	if (length(unique(booleanVec)) > 1) # means there are TRUE and FALSE
		booleanVec = rep(TRUE, n)
	return (booleanVec)
}

#### Tree data
## Load data
treeData = readRDS("./treeData_cleaned.rds")
names(treeData)

## Keep only useful columns and alive trees
treeData = unique(treeData[is_dead == "f", .(longitude, latitude, year_measured, species_id)])

## Keep only the 14 species of interest, comment the next two commands if you want the 106 species
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

treeData = treeData[species_id %in% ls_14species]
freq = setorderv(treeData[, .N, by = year_measured], "year_measured")

100 - freq[year_measured < 1997, sum(N)]*100/freq[, sum(N)]

## Keep only the last 22 years of presence recording
treeData = treeData[year_measured > 1996] # Keep 73.30% of the dataset
treeData[, year_measured := NULL]
treeData = unique(treeData)

#### Calculate presence/absence within treeData, and clean when dichotomic choice
## Presence/absence
treeData = data.table::dcast(treeData, formula = longitude + latitude ~ species_id)
treeData[, c(ls_14species) := lapply(.SD, function(x) (!is.na(x))), .SDcols = ls_14species]

## Change the colnames of treeData because of problems (- sign and starting by a number)
# Switch tsn and codename
species_dt = sortingSpecies(ls_14species)
setnames(x = treeData, old = species_dt[, species], new = species_dt[, reversed])

#### Aligne treeData on climate data
## Read a reference raster (already resampled)
ref_raster = raster("../clim60sec/clim_2010.grd")

## Reproject tree data on the same crs
treeData = sf::st_as_sf(treeData, coords = c("longitude", "latitude"),
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

treeData = sf::st_transform(treeData, crs = crs(ref_raster))

## Extract the climate coordinates
coords = sf::st_coordinates(treeData)
ref_indices = extract(ref_raster, coords, cellnumbers = TRUE)[, "cells"]
ref_coords = coordinates(ref_raster)[ref_indices, c("x", "y")]

## Coerce back treeData to data table
treeData = sf::st_drop_geometry(treeData)
treeData[, longitude := coords[, "X"]]
treeData[, latitude := coords[, "Y"]]

## Change coordinates treeData
treeData[, c("tree_lon", "tree_lat") := .(longitude, latitude)]
treeData[, longitude := ref_coords[, "x"]]
treeData[, latitude := ref_coords[, "y"]]

## Check the distance is less than sqrt(2) km, the resolution of the data (60 arc sec ~ 2 km²)
max(abs(treeData[, longitude - tree_lon]))
max(abs(treeData[, latitude - tree_lat]))
# For EPSG:2163, the unit is in meters https://epsg.io/2163
# Hence, it is ok

## Delete useless columns
treeData[, c("tree_lon", "tree_lat") := NULL]

#### Clear dichotomy appearing due to approximated coordinates
# Longitude-latitude appearing more than once
treeData = unique(treeData)
treeData[, count := .N, by = c("longitude", "latitude")]

treeData[count > 1, (species_dt[, reversed]) := lapply(.SD, delDichotomy),
	.SDcols = species_dt[, reversed], by = c("longitude", "latitude")]
treeData = unique(treeData)

treeData[, count := .N, by = c("longitude", "latitude")]
range(treeData[, count])
treeData[, count := NULL]

#### Save data.table
saveRDS(treeData, "./treeData_presence_absence.rds")
