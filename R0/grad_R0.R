
#### Aim of prog: Compute the gradient for each cell, and analyse the results for north and south
## Compute the gradient for each cell:
#	- Run the function sa2xy (inspired from the package rasterVis)
#
## Compute the average gradient within north and south:
#	- Compute the centroid of the cropped raster (important it is the cropped)
#	- North is the zone with latitudes above the centroid, and south the remaining part
#	- Compute the average of the arrows within each zone
#
## Remarks:
# The slope is the gradient, while the aspect is the gradient projected in the xy plan.
# See Paul Ritter (1987): A Vector-Based Slope and Aspect Generation Algorithm
# https://www.asprs.org/wp-content/uploads/pers/1987journal/aug/1987_aug_1109-1111.pdf

library(data.table)
library(raster)

options(max.print = 500)
rm(list = ls())

#### Tool function
## Compute x, y and, dx, dy
sa2xy = function(sa, scaleSlope, aspX, aspY)
{
    slope = subset(sa, 1)
    aspect = subset(sa, 2)

    # center = FALSE to get only positive values of slope
    if (scaleSlope)
        slope = scale(slope, center = FALSE)

    # sin due to the angular definition of aspect
    dx = slope * sin(aspect)
    dy = slope * cos(aspect)

    ## Returns a data.frame for panel.arrows; if reverse, put '-' sign
    dx = getValues(dx) * aspX
    dy = getValues(dy) * aspY
    x = getValues(init(sa, v = 'x'))
    y = getValues(init(sa, v = 'y'))
    data.frame(x, y, dx, dy)
}

## Qualitative direction arrow
dirArrow = function(dx, dy, reverse = FALSE)
{
	dirLat = ""
	dirLon = ""

	if (reverse)
	{
		dx = -dx
		dy = -dy
	}

	if (dy > 0)
		dirLat = "North"
	if (dy < 0)
		dirLat = "South"
	if (dx > 0)
		dirLon = "East"
	if (dx < 0)
		dirLon = "West"

	globalDir = paste0(dirLat, "-", dirLon)
	return (globalDir)
}

#### Common variables
## List folders
ls_folder = list.files(path = "./results/", pattern = "^[0-9]{4,}")
n = length(ls_folder)

## Data table averaged gradient north and south zones
average_gradient = data.table(species_id = ls_folder, dx_north_10 = numeric(n),
	dy_north_10 = numeric(n), dx_south_10 = numeric(n), dy_south_10 = numeric(n),
	dir_north_10 = character(n), dir_south_10 = character(n),
	dx_north_0 = numeric(n), dy_north_0 = numeric(n),
	dx_south_0 = numeric(n), dy_south_0 = numeric(n),
	dir_north_0 = character(n), dir_south_0 = character(n))

## Load bounding box determined by the data. Will be used to crop the results
data_bbox = readRDS("../createData/bbox_treeData.rds")
data_extent = extent(matrix(data_bbox, nrow = 2, ncol = 2))
north_extent = data_extent
south_extent = data_extent

#### Run for all the species
for (i in 1:n)
{
	## Load species-specific data
	species = ls_folder[i]
	folder = ls_folder[i]

	# R0 and climate (do not reorganise the rows before join)
	R0_10m = fread(paste0("./results/", folder, "/R0_10m.csv"))
	R0_0m = fread(paste0("./results/", folder, "/R0_0m.csv"))
	clim_2010 = readRDS(paste0("./Matlab_data/", folder, "/", species, "_clim2010.rds"))

	clim_2010 = clim_2010[, .(longitude, latitude)]
	clim_2010[, R0_10 := R0_10m]
	clim_2010[, R0_0 := R0_0m]

	## Set clim_2010 as raster
	clim_2010 = rasterFromXYZ(clim_2010,
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

	if (!identical(names(clim_2010), c("R0_10", "R0_0")))
	{
		print("For whaterver reason, R changed the names of the raster")
		names(clim_2010) = c("R0_10", "R0_0")
	}

	clim_2010 = crop(clim_2010, data_extent)

	# Compute centroid of cropped raster, get all coordinates
	centroid = colMeans(xyFromCell(clim_2010, which(getValues(is.na(clim_2010[["R0_10"]])) == 0))) # Using R0_0 shall not change anything

	## Northern zone
	# North
	north_extent@ymin = centroid["y"]
	north = crop(clim_2010, north_extent)

	# Compute slope and aspect in northern zone
	sa_north_10 = sa2xy(terrain(north[["R0_10"]], opt = c("slope", "aspect"), neighbors = 4),
    	scaleSlope = FALSE, aspX = 1, aspY = 1)

	sa_north_0 = sa2xy(terrain(north[["R0_0"]], opt = c("slope", "aspect"), neighbors = 4),
		scaleSlope = FALSE, aspX = 1, aspY = 1)

	setDT(sa_north_10)
	setDT(sa_north_0)
	sa_north_10 = sa_north_10[!is.na(dx)]
	sa_north_0 = sa_north_0[!is.na(dx)]

	print(paste0("Number of arrows: ", sa_north_10[, .N]))

	# Compute centroid and average
	centroid_north = colMeans(xyFromCell(north[["R0_0"]], which(getValues(is.na(north[["R0_0"]])) == 0))) # Using R0_0 shall not change anything
	north_average_10 = sa_north_10[, lapply(.SD, mean), .SDcols = c("dx", "dy")]
	north_average_0 = sa_north_0[, lapply(.SD, mean), .SDcols = c("dx", "dy")]

	## Southern zone
	# South
	south_extent@ymax = centroid["y"]
	south = crop(clim_2010, south_extent)

	# Compute slope and aspect in northern zone
	sa_south_10 = sa2xy(terrain(south[["R0_10"]], opt = c("slope", "aspect"), neighbors = 4),
    	scaleSlope = FALSE, aspX = 1, aspY = 1)

	sa_south_0 = sa2xy(terrain(south[["R0_0"]], opt = c("slope", "aspect"), neighbors = 4),
		scaleSlope = FALSE, aspX = 1, aspY = 1)

	setDT(sa_south_10)
	setDT(sa_south_0)
	sa_south_10 = sa_south_10[!is.na(dx)]
	sa_south_0 = sa_south_0[!is.na(dx)]

	print(paste0("Number of arrows: ", sa_south_10[, .N]))

	# Compute centroid and average
	centroid_south = colMeans(xyFromCell(south[["R0_0"]], which(getValues(is.na(south[["R0_0"]])) == 0))) # Using R0_0 shall not change anything
	south_average_10 = sa_south_10[, lapply(.SD, mean), .SDcols = c("dx", "dy")]
	south_average_0 = sa_south_0[, lapply(.SD, mean), .SDcols = c("dx", "dy")]

	## Combine and save results
	average_gradient[species_id == species, c("dx_north_10", "dy_north_10") := north_average_10]
	average_gradient[species_id == species, c("dx_north_0", "dy_north_0") := north_average_0]
	average_gradient[species_id == species, c("dx_south_10", "dy_south_10") := south_average_10]
	average_gradient[species_id == species, c("dx_south_0", "dy_south_0") := south_average_0]

	average_gradient[species_id == species, dir_north_10 := dirArrow(dx_north_10, dy_north_10, reverse = TRUE)]
	average_gradient[species_id == species, dir_north_0 := dirArrow(dx_north_0, dy_north_0, reverse = TRUE)]
	average_gradient[species_id == species, dir_south_10 := dirArrow(dx_south_10, dy_south_10, reverse = TRUE)]
	average_gradient[species_id == species, dir_south_0 := dirArrow(dx_south_0, dy_south_0, reverse = TRUE)]

	saveRDS(centroid, paste0("./results/", species, "/centroid.rds"))
	saveRDS(centroid_north, paste0("./results/", species, "/centroid_north.rds"))
	saveRDS(centroid_south, paste0("./results/", species, "/centroid_south.rds"))

	print(paste0("species: ", species, " done"))
}

saveRDS(average_gradient, "./average_gradient.rds")
