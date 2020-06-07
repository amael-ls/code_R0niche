
#### Aims of prog
## Describe the growth and mortality databases:
#	- How many data from Canada and USA respectively
#	- Climatic/geographic range
#	- Sample size
#	- Time frame
#
## Write latex table
#	- Sort by alphabetical order
#	- xtable package
#
## Figures
#
## Frequency of measurements
#	- Load growth data
#
## Remark:
# Due to projection error (non planar coordinates), some trees might be considered in both USA and Canada.
# Therefore they are counted twice, and I remove them from USA (apologise but there are only 538 trees!).

#### Load packages and tool functions; clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

library(data.table)
library(tikzDevice)
library(stringi)
library(xtable)
library(sf)

source("../toolFunctions.R")

#### Tool function
sortingSpecies = function(ls_species, decreasingOrder = FALSE)
{
	tsn_dt = data.table(species = character(length(ls_species)),
		tsn = character(length(ls_species))) # tsn = Taxonomic Serial Number

	tsn_dt[, species := stri_sub(str = ls_species,
		from = stri_locate_first(str = ls_species, regex = "-")[,1] + 1)]

	tsn_dt[, tsn := stri_sub(str = ls_species, from = 1,
		to = stri_locate_first(str = ls_species, regex = "-")[,1] - 1)]

	if (!decreasingOrder)
		setorderv(tsn_dt, "species")
	if (decreasingOrder)
		setorderv(tsn_dt, "species", -1)
	return (tsn_dt)
}

#### Load data
treeData = readRDS("treeClim_sStar.rds")

## Number of plots
nbPlots = length(unique(treeData[, plot_id]))
print(paste0("There are ", nbPlots, " plots in the database"))

## Keep exclusively columns of interest
treeData = unique(treeData[, .(tree_id, year_measured, species_id,
	longitude, latitude,
	annual_mean_temperature, annual_precipitation,
	min_temperature_of_coldest_month, precipitation_of_driest_quarter)])

#### Provenance of the data (USA/Canada)
## Load data
canada = readRDS("~/database/shapefiles/North_America_Shapefile/canada_union.rds")
usa = readRDS("~/database/shapefiles/North_America_Shapefile/usa_union.rds")

canada = st_geometry(canada)
usa = st_geometry(usa)

## Coerce to simple feature and reproject to EPSG:2163
treeData_sf = unique(treeData[, .(longitude, latitude)])
treeData_sf = st_as_sf(treeData_sf, coords = c("longitude", "latitude"),
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

treeData_sf = st_transform(treeData_sf, crs = st_crs(usa))

## Which points are in Canada, and which are in the USA. Those on the borderline will be in Canada
in_can = st_contains(canada, treeData_sf, sparse = FALSE)
in_usa = st_contains(usa, treeData_sf, sparse = FALSE)

saveRDS(in_can, "./in_can.rds")
saveRDS(in_usa, "./in_usa.rds")

countedTwice = sum(in_can) + sum(in_usa) - nrow(treeData_sf) # If negative, then there are forgotten rather than counted twice

print(paste0("There are ", sum(in_can), " in Canada"))
print(paste0("There are ", sum(in_usa) - countedTwice, " in USA"))

print(paste0("There are ", round(sum(in_can)*100/nrow(treeData_sf), 2), "% in Canada"))
print(paste0("There are ", round((sum(in_usa) - countedTwice)*100/nrow(treeData_sf), 2), "% in USA"))

#### Climatic/geographic range for each species
## Keep only 14 species
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

treeData = treeData[species_id %in% ls_14species]

## Coerce to simple feature and reproject to EPSG:2163
treeData_sf = unique(treeData[, .(longitude, latitude)])
treeData_sf = st_as_sf(treeData_sf, coords = c("longitude", "latitude"),
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

treeData_sf = st_transform(treeData_sf, crs = st_crs(usa))

## Which points are in Canada, and which are in the USA. Those on the borderline will be in Canada
in_can = st_contains(canada, treeData_sf, sparse = FALSE)
in_usa = st_contains(usa, treeData_sf, sparse = FALSE)

saveRDS(in_can, "./in_can2.rds")
saveRDS(in_usa, "./in_usa2.rds")

countedTwice = sum(in_can) + sum(in_usa) - nrow(treeData_sf) # If negative, then there are forgotten rather than counted twice

print(paste0("There are ", sum(in_can), " in Canada"))
print(paste0("There are ", sum(in_usa) - countedTwice, " in USA"))

print(paste0("There are ", round(sum(in_can)*100/nrow(treeData_sf), 2), "% in Canada"))
print(paste0("There are ", round((sum(in_usa) - countedTwice)*100/nrow(treeData_sf), 2), "% in USA"))

## Climatic range
climVar = c("annual_mean_temperature", "min_temperature_of_coldest_month",
	"annual_precipitation", "precipitation_of_driest_quarter")

treeData[, paste0(climVar, "_min") := lapply(.SD, min), .SDcol = climVar, by = species_id]
treeData[, paste0(climVar, "_max") := lapply(.SD, max), .SDcol = climVar, by = species_id]
treeData[, paste0(climVar, "_mean") := lapply(.SD, mean), .SDcol = climVar, by = species_id]

## Geographic range
space = c("longitude", "latitude")
treeData[, paste0(space, "_min") := lapply(.SD, min), .SDcol = space, by = species_id]
treeData[, paste0(space, "_max") := lapply(.SD, max), .SDcol = space, by = species_id]

#### Number of individuals per species
treeData[, "nbIndiv" := .N, by = species_id]

#### Time frame
treeData[, "year_min" := lapply(.SD, min), .SDcol = "year_measured", by = species_id]
treeData[, "year_max" := lapply(.SD, max), .SDcol = "year_measured", by = species_id]

#### Remove all the useless columns, and do unique
treeData[, c("tree_id", "year_measured", climVar, space) := NULL]
treeData = unique(treeData)

treeData[, .N] == length(ls_14species)

#### Write latex tables
## Sort by alphabetical order, remove tsn
alphaSpecies = sortingSpecies(ls_14species)
alphaSpecies[, species_id := paste(tsn, species, sep = "-")]

treeData = treeData[alphaSpecies, on = "species_id"]

## Save in rds format
saveRDS(treeData, "./speciesSpecificInfo.rds")

## Create the xtable
describe_xt = xtable(treeData[, .(species_id, nbIndiv,
	longitude_min, longitude_max, latitude_min, latitude_max)])
print(describe_xt, file = "./describe.tex", include.rownames = FALSE, booktabs = TRUE)

## Calculate min, max, total
max_sp = treeData[nbIndiv == max(nbIndiv), species_id]
min_sp = treeData[nbIndiv == min(nbIndiv), species_id]
print(paste0("Most measured species: ", max_sp, " with ", treeData[species_id == max_sp, nbIndiv], " individuals"))
print(paste0("Least measured species: ", min_sp, " with ", treeData[species_id == min_sp, nbIndiv], " individuals"))
print(paste0("Size of database: ", sum(treeData[, nbIndiv]), " individual measurements"))

#### Climate range
## Label x axis
xLabels = stri_replace_all(climVar, replacement =  " ", regex = "_")

## Decreasing species order of treeData for plot
setorderv(treeData, "species", -1)

## Plots
for (i in seq(1, length(climVar), by = 2))
{
	# Bounds of plot for climVar
	min_bound = min(treeData[, paste0(climVar[i], "_min"), with = FALSE])
	max_bound = max(treeData[, paste0(climVar[i], "_max"), with = FALSE])
	av_global = mean(unlist(treeData[, paste0(climVar[i], "_mean"), with = FALSE]))

	tikz(paste0("./", climVar[i], "-",climVar[i + 1], "_sp_range.tex"),
		width = 3.1, height = 3.1)
	op = par(mar = c(2.5, 5, 0.1, 0.4), mgp = c(1.5, 0.2, 0),
		oma=c(0,0,0.9,0), tck = -0.015, las = 1)
	plot(x = NULL, y = NULL, xlim = c(min_bound, max_bound),
		ylim = c(1, length(ls_14species)), axes = FALSE, xlab = xLabels[i],
		ylab = "")

	if (stri_detect(str = climVar[i], regex = "precipitation"))
		pos = round(seq(min_bound, max_bound, length.out = 5), 0)
	if (stri_detect(str = climVar[i], regex = "temperature"))
		pos = round(seq(min_bound, max_bound, length.out = 5), 1)
	axis(side = 1, at = pos)
	axis(side = 2, at = seq(1, length(ls_14species), length.out = length(ls_14species)),
		labels = treeData[, species], tck = 0, lwd = "")
	for (j in 1:length(ls_14species))
	{
		m = unlist(treeData[j, paste0(climVar[i], "_min"), with = FALSE])
		M = unlist(treeData[j, paste0(climVar[i], "_max"), with = FALSE])
		a = unlist(treeData[j, paste0(climVar[i], "_mean"), with = FALSE])

		segments(x0 = m, y0 = j, x1 = M, y1 = j, lwd = 2)
		points(x = a, y = j, pch = 19)
	}
	abline(v = av_global, lwd = 2, lty = 5)

	i = i + 1
	# Bounds of plot for climVar
	min_bound = min(treeData[, paste0(climVar[i], "_min"), with = FALSE])
	max_bound = max(treeData[, paste0(climVar[i], "_max"), with = FALSE])
	av_global = mean(unlist(treeData[, paste0(climVar[i], "_mean"), with = FALSE]))

	op = par(mar = c(2.5, 5, 0.1, 0.4), mgp = c(1.5, 0.2, 0),
		oma=c(0,0,0.9,0), tck = -0.015, las = 1)
	plot(x = NULL, y = NULL, xlim = c(min_bound, max_bound),
		ylim = c(1, length(ls_14species)), axes = FALSE, xlab = xLabels[i],
		ylab = "")

	if (stri_detect(str = climVar[i], regex = "precipitation"))
		pos = round(seq(min_bound, max_bound, length.out = 5), 0)
	if (stri_detect(str = climVar[i], regex = "temperature"))
		pos = round(seq(min_bound, max_bound, length.out = 5), 1)
	axis(side = 1, at = pos)
	for (j in 1:length(ls_14species))
	{
		m = unlist(treeData[j, paste0(climVar[i], "_min"), with = FALSE])
		M = unlist(treeData[j, paste0(climVar[i], "_max"), with = FALSE])
		a = unlist(treeData[j, paste0(climVar[i], "_mean"), with = FALSE])

		segments(x0 = m, y0 = j, x1 = M, y1 = j, lwd = 2)
		points(x = a, y = j, pch = 19)
	}
	abline(v = av_global, lwd = 2, lty = 5)
	dev.off()
}

#### Frequency of measurements
## Growth
growth_dt = readRDS("./growth_dt.rds")
sum(growth_dt[(3 <= deltaYear) & (deltaYear <= 10), table(deltaYear)]*100/growth_dt[, .N])
sum(growth_dt[(3 <= deltaYear) & (deltaYear <= 7), table(deltaYear)]*100/growth_dt[, .N])
sum(growth_dt[(3 <= deltaYear) & (deltaYear <= 15), table(deltaYear)]*100/growth_dt[, .N])
sum(growth_dt[deltaYear == 5, table(deltaYear)]*100/growth_dt[, .N])
growth_dt[, table(deltaYear)]

## Mortality
mortality_dt = readRDS("./mortality_dt.rds")
sum(mortality_dt[(3 <= deltaYear) & (deltaYear <= 10), table(deltaYear)]*100/mortality_dt[, .N])
sum(mortality_dt[(3 <= deltaYear) & (deltaYear <= 7), table(deltaYear)]*100/mortality_dt[, .N])
sum(mortality_dt[(5 <= deltaYear) & (deltaYear <= 11), table(deltaYear)]*100/mortality_dt[, .N])
sum(mortality_dt[deltaYear == 5, table(deltaYear)]*100/mortality_dt[, .N])
mortality_dt[, table(deltaYear)]

#### linear model s* with latitude
## Keep variable of interest of growth_dt
growth_dt = unique(growth_dt[, .(plot_id, year_measured, latitude, longitude, s_star)])
growth_dt[, mean(s_star)]
setDF(growth_dt)
s_star_lat_lm = lm(formula = s_star ~ latitude, data = growth_dt)
summary(s_star_lat_lm)
