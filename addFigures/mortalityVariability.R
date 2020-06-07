
#### Aim of prog: Show the variability of tree mortality using histograms
## Description:
# There will be three histogram per species:
#	- northern part
#	- middle part
#	- southern part
#
## Remark:
# I did a quick and dirty function using global access variable

#### Load package and clear memory
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool function
## Get centroids (centro, north, and south) of data
computeCentroLat = function(lat)
{
	lat = unique(lat)

	min_lat = min(lat)
	max_lat = max(lat)
	av_lat = mean(lat)

	return(list(centro_lat = av_lat, north_lat = mean(c(av_lat, max_lat)), south_lat = mean(c(av_lat, min_lat))))
}

## Set region
region = function(lat, centro, north, south)
{
	results = rep("middle", length(lat))
	results[lat <= south] = "south"
	results[lat >= north] = "north"

	return (list(region = results))
}

## Function to plot (quick and dirty, unsing global access variable)
varMortality_fct = function(nameFig, rangeFig = 1:3, stand = FALSE)
{
	tikz(nameFig, width = 6, height = 5, standAlone = stand)
	op = par(mar = c(0, 3, 4, 0), mgp = c(1.5, 0.75, 0),
		oma = c(0,2,0,0), tck = -0.015, las = 1)

	tikzAnnotate(paste0("\\definecolor{shadeN}{RGB}{", paste(shadeN, collapse = ","), "}"))
	tikzAnnotate(paste0("\\definecolor{shadeM}{RGB}{", paste(shadeM, collapse = ","), "}"))
	tikzAnnotate(paste0("\\definecolor{shadeS}{RGB}{", paste(shadeS, collapse = ","), "}"))

	plot(x = NULL, y = NULL, xlim = c(0, 3*speciesSpace + 2*interSpecies + 0.1),
		ylim = c(0, maxM), axes = FALSE, xlab = "",
		ylab = "")

	count = 0
	for (i in rangeFig)
	{
		## Species line and name
		sp = ls_species[i]
		species_coords = halfSpeciesSpace + count*(speciesSpace + interSpecies)
		tikzCoord(species_coords - halfSpeciesSpace, maxM + 0.01, paste0("pt_start_", i))
		tikzCoord(species_coords + halfSpeciesSpace, maxM + 0.01, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(species_coords, maxM + 0.01, 'species_coords')

		tikzAnnotate(paste0("\\node[above] (spPos) at (species_coords) {", sp, "};"))

		## Add poxplots manually
		# North
		north_prop = mortality_dt[(species == sp) & (region == "north"), unique(prop)]
		north_right = species_coords - midBox - interBox
		north_left = north_right - widthBox

		tikzCoord(north_left, 0, paste0("lower_north_", i)) # Lower left corner
		tikzCoord(north_right, north_prop, paste0("upper_north_", i)) # Upper right corner
		tikzAnnotate(paste0("\\draw[color=shadeN] (lower_north_", i, ") rectangle (upper_north_", i, ");")) # Rectangle

		# Middle
		middle_prop = mortality_dt[(species == sp) & (region == "middle"), unique(prop)]
		middle_right = species_coords + midBox
		middle_left = species_coords - midBox

		tikzCoord(middle_left, 0, paste0("lower_middle_", i)) # Lower left corner
		tikzCoord(middle_right, middle_prop, paste0("upper_middle_", i)) # Upper right corner
		tikzAnnotate(paste0("\\draw[color=shadeM] (lower_middle_", i, ") rectangle (upper_middle_", i, ");")) # Rectangle

		# South
		south_prop = mortality_dt[(species == sp) & (region == "south"), unique(prop)]
		south_right = species_coords + midBox + interBox
		south_left = south_right + widthBox

		tikzCoord(south_left, 0, paste0("lower_south_", i)) # Lower left corner
		tikzCoord(south_right, south_prop, paste0("upper_south_", i)) # Upper right corner
		tikzAnnotate(paste0("\\draw[color=shadeS] (lower_south_", i, ") rectangle (upper_south_", i, ");")) # Rectangle

		count = count + 1
	}

	axis(side = 2, at = seq(0, maxM, by = 0.05))
	mtext(text = "Relative mortality proportion", side = 2, outer = TRUE, las = 0)

	# Legend
	tikzAnnotate("
		\\matrix [below right] at (current bounding box.north west) {
			\\node [shape = rectangle, fill = shadeN, label = right:North] {}; &
			\\node [shape = rectangle, fill = shadeM, label = right:Middle] {}; &
			\\node [shape = rectangle, fill = shadeS, label = right:South] {}; \\\\
		};
	")
	dev.off()
}

## Round number to the closest .05 number
round05 = function(x, delta = 0.05)
{
	inf = floor(x)
	sup = ceiling(x)
	vec = seq(inf, sup, by = delta)

	return(vec[which.min(abs(x - vec))]) # Will give the first solution if non-uniqueness
}

#### Load data and compute data centroids
## Mortality data
mortality_dt = readRDS("../createData/mortality_dt.rds")[, .(species_id, latitude, deltaState)]

tsn = readRDS("../growth/tsn.rds")[, .(species_id, species, tolLevel)]
ls_species = sort(tsn[, species])
n = length(ls_species)

## Species-specific data centroids
# Latitude centroids
mortality_dt[, c("centroid", "north", "south") := computeCentroLat(latitude), by = species_id]
mortality_dt[, region := region(latitude, centroid, north, south), by = species_id]
mortality_dt[, table(region), by = species_id]

mortality_dt = mortality_dt[tsn, on = "species_id"]

## Proportion of death
mortality_dt[, prop := sum(deltaState)/.N, by = c("species", "region")]
maxM = round05(mortality_dt[, max(prop)])

#### Plot parameters
## Colours for Northern, Middle and Southern regions
shadeN = col2rgb("#010120")[,1]
shadeM = col2rgb("#0E51FF")[,1]
shadeS = col2rgb("#00F9FF")[,1]

## Space and width
widthBox = 0.75
interBox = 0.25
interSpecies = 1

midBox = widthBox/2
speciesSpace = 3*widthBox + 2*interBox
halfSpeciesSpace = speciesSpace/2

#### Plot
varMortality_fct("mortalityVar1-3.tex", rangeFig = 1:3) # , stand = TRUE)
varMortality_fct("mortalityVar4-6.tex", rangeFig = 4:6) # , stand = TRUE)
varMortality_fct("mortalityVar7-9.tex", rangeFig = 7:9) # , stand = TRUE)
varMortality_fct("mortalityVar10-12.tex", rangeFig = 10:12) # , stand = TRUE)
varMortality_fct("mortalityVar13-14.tex", rangeFig = 13:14) # , stand = TRUE)
