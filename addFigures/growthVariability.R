
#### Aim of prog: Show the variability of tree growth using boxplots
## Description:
# There will be three boxplots per species:
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
varGrowth_fct = function(nameFig, rangeFig = 1:3, stand = FALSE)
{
	tikz(nameFig, width = 6, height = 5, standAlone = stand)
	op = par(mar = c(0, 2.5, 2, 0), mgp = c(1.5, 0.75, 0),
		oma = c(0,2,0,0), tck = -0.015, las = 1)

	tikzAnnotate(paste0("\\definecolor{shadeN}{RGB}{", paste(shadeN, collapse = ","), "}"))
	tikzAnnotate(paste0("\\definecolor{shadeM}{RGB}{", paste(shadeM, collapse = ","), "}"))
	tikzAnnotate(paste0("\\definecolor{shadeS}{RGB}{", paste(shadeS, collapse = ","), "}"))

	plot(x = NULL, y = NULL, xlim = c(0, 3*speciesSpace + 2*interSpecies + 0.1),
		ylim = c(0, maxG + 1), axes = FALSE, xlab = "",
		ylab = "")

	count = 0
	for (i in rangeFig)
	{
		## Species line and name
		sp = ls_species[i]
		species_coords = halfSpeciesSpace + count*(speciesSpace + interSpecies)
		tikzCoord(species_coords - halfSpeciesSpace, maxG + 0.5, paste0("pt_start_", i))
		tikzCoord(species_coords + halfSpeciesSpace, maxG + 0.5, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(species_coords, maxG + 0.5, 'species_coords')

		tikzAnnotate(paste0("\\node[above] (spPos) at (species_coords) {", sp, "};"))

		## Add boxplots manually
		# North
		north_qt = growth_dt[(species == sp) & (region == "north"), quantile(growth, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))]
		outliers_above = growth_dt[(species == sp) & (region == "north") & (growth > north_qt["97.5%"]), unique(growth)]
		outliers_below = growth_dt[(species == sp) & (region == "north") & (growth < north_qt["2.5%"]), unique(growth)]
		north_right = species_coords - midBox - interBox
		north_left = north_right - widthBox

		segments(x0 = north_left, y0 = north_qt["2.5%"], x1 = north_right, y1 = north_qt["2.5%"], col = "#010120") # Lower horizontal
		segments(x0 = species_coords - interBox - widthBox, y0 = north_qt["2.5%"],
			x1 = species_coords - interBox - widthBox, y1 = north_qt["25%"], col = "#010120") # Lower vertical

		tikzCoord(north_left, north_qt["25%"], paste0("lower_north_", i)) # Lower left corner
		tikzCoord(north_right, north_qt["75%"], paste0("upper_north_", i)) # Upper right corner
		tikzAnnotate(paste0("\\draw[color=shadeN] (lower_north_", i, ") rectangle (upper_north_", i, ");")) # Rectangle

		segments(x0 = north_left, y0 = north_qt["50%"], x1 = north_right, y1 = north_qt["50%"], col = "#010120") # Lower horizontal

		segments(x0 = species_coords - interBox - widthBox, y0 = north_qt["75%"],
			x1 = species_coords - interBox - widthBox, y1 = north_qt["97.5%"], col = "#010120") # Upper vertical
		segments(x0 = north_left, y0 = north_qt["97.5%"], x1 = north_right, y1 = north_qt["97.5%"], col = "#010120") # Upper horizontal

		points(x = rep(species_coords - interBox - widthBox, length(outliers_above)),
			y = outliers_above, col = "#010120", pch = 20)

		points(x = rep(species_coords - interBox - widthBox, length(outliers_below)),
			y = outliers_below, col = "#010120", pch = 20)

		# Middle
		middle_qt = growth_dt[(species == sp) & (region == "middle"), quantile(growth, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))]
		outliers_above = growth_dt[(species == sp) & (region == "middle") & (growth > middle_qt["97.5%"]), unique(growth)]
		outliers_below = growth_dt[(species == sp) & (region == "middle") & (growth < middle_qt["2.5%"]), unique(growth)]
		middle_right = species_coords + midBox
		middle_left = species_coords - midBox

		segments(x0 = middle_left, y0 = middle_qt["2.5%"], x1 = middle_right, y1 = middle_qt["2.5%"], col = "#0E51FF") # Lower horizontal
		segments(x0 = species_coords, y0 = middle_qt["2.5%"], x1 = species_coords, y1 = middle_qt["25%"], col = "#0E51FF") # Lower vertical

		tikzCoord(middle_left, middle_qt["25%"], paste0("lower_middle_", i)) # Lower left corner
		tikzCoord(middle_right, middle_qt["75%"], paste0("upper_middle_", i)) # Upper right corner
		tikzAnnotate(paste0("\\draw[color=shadeM] (lower_middle_", i, ") rectangle (upper_middle_", i, ");")) # Rectangle

		segments(x0 = middle_left, y0 = middle_qt["50%"], x1 = middle_right, y1 = middle_qt["50%"], col = "#0E51FF") # Lower horizontal

		segments(x0 = species_coords, y0 = middle_qt["75%"], x1 = species_coords, y1 = middle_qt["97.5%"], col = "#0E51FF") # Upper vertical
		segments(x0 = middle_left, y0 = middle_qt["97.5%"], x1 = middle_right, y1 = middle_qt["97.5%"], col = "#0E51FF") # Upper horizontal

		points(x = rep(species_coords, length(outliers_above)),
			y = outliers_above, col = "#0E51FF", pch = 20)

		points(x = rep(species_coords, length(outliers_below)),
			y = outliers_below, col = "#0E51FF", pch = 20)

		# South
		south_qt = growth_dt[(species == sp) & (region == "south"), quantile(growth, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))]
		outliers_above = growth_dt[(species == sp) & (region == "south") & (growth > south_qt["97.5%"]), unique(growth)]
		outliers_below = growth_dt[(species == sp) & (region == "south") & (growth < south_qt["2.5%"]), unique(growth)]
		south_right = species_coords + midBox + interBox
		south_left = south_right + widthBox

		segments(x0 = south_left, y0 = south_qt["2.5%"], x1 = south_right, y1 = south_qt["2.5%"], col = "#00F9FF") # Lower horizontal
		segments(x0 = species_coords + interBox + widthBox, y0 = south_qt["2.5%"],
			x1 = species_coords + interBox + widthBox, y1 = south_qt["25%"], col = "#00F9FF") # Lower vertical

		tikzCoord(south_left, south_qt["25%"], paste0("lower_south_", i)) # Lower left corner
		tikzCoord(south_right, south_qt["75%"], paste0("upper_south_", i)) # Upper right corner
		tikzAnnotate(paste0("\\draw[color=shadeS] (lower_south_", i, ") rectangle (upper_south_", i, ");")) # Rectangle

		segments(x0 = south_left, y0 = south_qt["50%"], x1 = south_right, y1 = south_qt["50%"], col = "#00F9FF") # Lower horizontal

		segments(x0 = species_coords + interBox + widthBox, y0 = south_qt["75%"],
			x1 = species_coords + interBox + widthBox, y1 = south_qt["97.5%"], col = "#00F9FF") # Upper vertical
		segments(x0 = south_left, y0 = south_qt["97.5%"], x1 = south_right, y1 = south_qt["97.5%"], col = "#00F9FF") # Upper horizontal

		points(x = rep(species_coords + interBox + widthBox, length(outliers_above)),
			y = outliers_above, col = "#00F9FF", pch = 20)

		points(x = rep(species_coords + interBox + widthBox, length(outliers_below)),
			y = outliers_below, col = "#00F9FF", pch = 20)

		count = count + 1
	}

	axis(side = 2, at = seq(0, maxG, by = 5))
	mtext(text = "Growth data (in mm/yr)", side = 2, outer = TRUE, las = 0)

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

#### Load data and compute data centroids
## Growth data
growth_dt = readRDS("../createData/growth_dt.rds")[, .(species_id, latitude, growth)]

tsn = readRDS("../growth/tsn.rds")[, .(species_id, species, tolLevel)]
ls_species = sort(tsn[, species])
n = length(ls_species)

## Species-specific data centroids
# Latitude centroids
growth_dt[, c("centroid", "north", "south") := computeCentroLat(latitude), by = species_id]
growth_dt[, region := region(latitude, centroid, north, south), by = species_id]
growth_dt[, table(region), by = species_id]

growth_dt = growth_dt[tsn, on = "species_id"]

maxG = growth_dt[, max(growth)]

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
varGrowth_fct("growthVar1-3.tex", rangeFig = 1:3) # , stand = TRUE)
varGrowth_fct("growthVar4-6.tex", rangeFig = 4:6) # , stand = TRUE)
varGrowth_fct("growthVar7-9.tex", rangeFig = 7:9) # , stand = TRUE)
varGrowth_fct("growthVar10-12.tex", rangeFig = 10:12) # , stand = TRUE)
varGrowth_fct("growthVar13-14.tex", rangeFig = 13:14) # , stand = TRUE)
