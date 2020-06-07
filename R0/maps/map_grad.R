
#### Aim of prog: Make maps from matlab results
library(data.table)
library(stringi)
library(scales)
library(raster)
library(fields)
library(sf)

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

## Get the neighbours values of a cell from a raster
getNeighboursVal = function(val_raster, id, ncell_r, ncol_r)
{
	# Get indices of the neighbours
	ind = c(id - ncol_r - 1, id - ncol_r, id - ncol_r + 1,
		id - 1, id + 1,
		id + ncol_r - 1, id + ncol_r, id + ncol_r + 1)

	# Remove indices outside of raster (happen when working on the edge)
	ind = ind[(ind > 0) & (ind < ncell_r)]
	return (val_raster[ind])
}

## Fill the gaps in raster, using neighbourhood
fillGap = function(r, p)
{
	vals = getValues(r)
	vals_orig = vals # Needs it
	ind_na_orig = sort(which(is.na(vals)))

	# Get values in polygon p
	ind_p = cellFromPolygon(r, p) # list of indices for each polygon (if p is a multipolygon)
	ind_p = sort(unlist(ind_p))
	ind_p = unique(ind_p)

	vals = vals[ind_p]

	# Update ind_na_orig, and create a new ind_na within polygon
	ind_na_orig = ind_p[ind_p %in% ind_na_orig]
	ind_na = sort(which(is.na(vals)))
	ind_na_orig_withinPolygon = ind_na

	ncell_r = ncell(r)
	ncol_r = ncol(r)

	unstable = TRUE
	count = 0

	ind_dt = data.table(ind_val = ind_na, ind_orig = ind_na_orig)

	while (length(ind_na) != 0 & unstable)
	{
		prev_stage = vals
		print(paste0("Number of NAs: ", length(ind_na)))
		for (i in ind_na)
		{
			neighVals = getNeighboursVal(vals_orig, ind_dt[ind_val == i, ind_orig], ncell_r, ncol_r)
			if (sum(is.na(neighVals)) == 8)
				next;
			vals[i] = mean(neighVals, na.rm = TRUE)
		}

		vals_orig[ind_dt[ind_val %in% ind_na, ind_orig]] = vals[ind_na]

		ind_na = sort(which(is.na(vals)))

		ind_dt = ind_dt[ind_val %in% ind_na]

		if(identical(prev_stage, vals))
		{
			if (length (ind_na) == 1)
			{
				unstable = FALSE
			} else {
				count = count + 1
				ind_na = sample(ind_na)
				if (count == 10)
					unstable = FALSE
			}
		}
	}

	r[ind_na_orig] = vals[ind_na_orig_withinPolygon]
	return (r)
}

#### Common variables
## List folders
ls_folder = list.files(path = "../results/", pattern = "^[0-9]{4,}")

height_canopy = "10m"

## Data table centroid and gradient
average_gradient = readRDS("../average_gradient.rds")

## Shapefiles Northern America
canada = readRDS("~/database/shapefiles/North_America_Shapefile/canada_union.rds")
usa = readRDS("~/database/shapefiles/North_America_Shapefile/usa_union.rds")

erie = st_read("~/database/shapefiles/North_America_Shapefile/lakes/erie")
huron = st_read("~/database/shapefiles/North_America_Shapefile/lakes/huron")
michigan = st_read("~/database/shapefiles/North_America_Shapefile/lakes/michigan")
ontario = st_read("~/database/shapefiles/North_America_Shapefile/lakes/ontario")
stClair = st_read("~/database/shapefiles/North_America_Shapefile/lakes/stClair")
superior = st_read("~/database/shapefiles/North_America_Shapefile/lakes/superior")

canada = st_transform(canada,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
usa = st_transform(usa,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
erie = st_transform(erie,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
huron = st_transform(huron,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
michigan = st_transform(michigan,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
ontario = st_transform(ontario,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
stClair = st_transform(stClair,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
superior = st_transform(superior,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

canada = st_geometry(canada)
usa = st_geometry(usa)
erie = st_geometry(erie)
huron = st_geometry(huron)
michigan = st_geometry(michigan)
ontario = st_geometry(ontario)
stClair = st_geometry(stClair)
superior = st_geometry(superior)

## Load bounding box determined by the data. Will be used to crop the results
data_bbox = readRDS("../../createData/bbox_treeData.rds")
data_extent = extent(matrix(data_bbox, nrow = 2, ncol = 2))

#### Run for all the species
## Plot maps
for (i in 1:length(ls_folder))
{
	## Load species-specific data
	species = ls_folder[i]
	folder = ls_folder[i]

	# R0 and climate (do not reorganise the rows before join)
	R0_vec = fread(paste0("../results/", folder, "/R0_", height_canopy, ".csv"))
	clim_2010 = readRDS(paste0("../Matlab_data/", folder, "/", species, "_clim2010.rds"))

	clim_2010 = clim_2010[, .(longitude, latitude)]
	clim_2010[, R0 := R0_vec]

	# Load the Little's 1971 shapefile, crop to data bbox
	ls_little = list.files(path = paste0("../../little1971/", folder), pattern = ".shp$")
	little1971 = st_read(paste0("../../little1971/", folder, "/", ls_little))
	little1971 = st_transform(little1971,
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
	little1971 = st_crop(little1971, data_bbox)
	little1971 = st_geometry(little1971)
	little1971_sp = as(little1971, "Spatial")

	# Centroids and gradients
	centroid_north = readRDS(paste0("../results/", folder, "/centroid_north.rds"))
	centroid_south = readRDS(paste0("../results/", folder, "/centroid_south.rds"))
	centroid = readRDS(paste0("../results/", folder, "/centroid.rds"))

	if (height_canopy == "10m")
	{
		dx_north = average_gradient[species_id == species, dx_north_10]
		dy_north = average_gradient[species_id == species, dy_north_10]
		dx_south = average_gradient[species_id == species, dx_south_10]
		dy_south = average_gradient[species_id == species, dy_south_10]
	}

	if (height_canopy == "0m")
	{
		dx_north = average_gradient[species_id == species, dx_north_0]
		dy_north = average_gradient[species_id == species, dy_north_0]
		dx_south = average_gradient[species_id == species, dx_south_0]
		dy_south = average_gradient[species_id == species, dy_south_0]
	}

	## Set clim_2010 as raster
	clim_2010 = rasterFromXYZ(clim_2010,
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

	if (names(clim_2010) != "R0")
	{
		print("For whaterver reason, R changed the name")
		names(clim_2010) = "R0"
	}

	clim_2010 = crop(clim_2010, data_extent)

	## Fill the NA within little1971
	clim_2010 = fillGap(clim_2010, little1971_sp)

	## Compute slope
	# Sample regularly (otherwise too many arrows)
	set.seed(1969-08-18) # Woodstock seed
	dat = sampleRegular(clim_2010, size = 300, asRaster = TRUE)

	sa = sa2xy(terrain(dat, opt = c("slope", "aspect"), neighbors = 4),
		scaleSlope = FALSE, aspX = 1, aspY = 1)

	setDT(sa)
	print(paste0("Number of arrows: ", sa[!is.na(dx), .N]))
	sa = sa[!is.na(dx)]

	## Calculate rho_0
	minR0 = cellStats(x = clim_2010, stat = "min", na.rm = TRUE)
	maxR0 = cellStats(x = clim_2010, stat = "max", na.rm = TRUE)
	clim_2010$R0_scaled = (clim_2010$R0 - minR0)/(maxR0 - minR0)

	# Save lon-lat, R0, rho_0 sp-specific data
	path_lonLat_rs = paste0("../results/", folder, "/", height_canopy, "/")
	if (!dir.exists(path_lonLat_rs))
		dir.create(path_lonLat_rs)

	writeRaster(clim_2010, filename = paste0(path_lonLat_rs, "lonLatR0cropped.grd"), bandorder='BIL', overwrite = TRUE)

	## Bounding box of plot
	lonMin = unname(data_bbox["xmin"]) - 100
	lonMax = unname(data_bbox["xmax"]) + 30000 # To get Nova Scotia completely
	latMin = unname(data_bbox["ymin"]) - 60000 # To get Florida completely
	latMax = unname(data_bbox["ymax"]) + 100

	## Spatial plot of R0(h* = 10m)
	# Format = jpeg. arrow.plot is reverse to get direction of increase (-sa[, .(dx, dy)])
	jpeg(paste0("./", species, "_R0_h=", height_canopy, "_scaled.jpg"), width = 1000, height = 1000, quality = 100)
	plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
		xlab = "", ylab = "", bg = "transparent")
	plot(clim_2010[["R0_scaled"]], cex = 0.5, pch = 19, add = TRUE, legend = FALSE, # , "#4269E2"
		col = c("#000000", "#2A3344", "#1122AA", "#2058DC", "#5A82EA", "#FDECBE", "#FDDB5B", "#FAB935", "#FD9859", "#B6343A"))
	plot(little1971, col = NA, border = "#000000", add = TRUE, lwd = 2) # "#DD925C"
	plot(canada, col = NA, add = TRUE, lwd = 2)
	plot(usa, col = NA, add = TRUE, lwd = 2)
	plot(erie, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(huron, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(michigan, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(ontario, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(stClair, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(superior, col = "#FFFFFF", border = "#000000", add = TRUE)
	arrow.plot(a1 = sa[, x], a2 = sa[, y], u = -sa[, dx], v = -sa[, dy],
		true.angle = TRUE, arrowfun = sfsmisc::p.arrows, fill = "white", lty = "blank")
	arrow.plot(a1 = centroid_north["x"], a2 = centroid_north["y"], u = -dx_north, v = -dy_north,
		true.angle = TRUE, arrowfun = sfsmisc::p.arrows, fill = "red", size = 2, lty = "blank")
	arrow.plot(a1 = centroid_south["x"], a2 = centroid_south["y"], u = -dx_south, v = -dy_south,
		true.angle = TRUE, arrowfun = sfsmisc::p.arrows, fill = "red", size = 2, lty = "blank")
	points(x = centroid["x"], y = centroid["y"], pch = 15, cex = 3, col = "#00ff00")
	# text(x = centroid["x"], y = centroid["y"], labels = "C", pos = 4, cex = 2.5, offset = 0.8, family = "CM Roman")
	dev.off()

	## Spatial plot of local s_inf
	# Format = jpeg
	if (height_canopy == "10m")
	{
		local_s_inf = readRDS(paste0("../Matlab_data/", folder, "/", species, "_clim2010.rds"))[, .(longitude, latitude)]
		local_s_inf[, s_inf := fread(paste0("../results/", folder, "/local_s_inf.csv"))]
		local_s_inf = rasterFromXYZ(local_s_inf,
			crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
		local_s_inf = crop(local_s_inf, data_extent)

		jpeg(paste0("./", species, "_local_s_inf.jpg"), width = 1000, height = 1000, quality = 100)
		plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
			xlab = "", ylab = "", bg = "transparent")
		plot(local_s_inf[["s_inf"]], cex = 0.5, pch = 19, add = TRUE, legend = FALSE, # , "#4269E2"
			col = c("#000000", "#2A3344", "#1122AA", "#2058DC", "#5A82EA", "#FDECBE", "#FDDB5B", "#FAB935", "#FD9859", "#B6343A"))
		plot(little1971, col = NA, border = "#000000", add = TRUE, lwd = 2)
		plot(canada, col = NA, add = TRUE, lwd = 2)
		plot(usa, col = NA, add = TRUE, lwd = 2)
		plot(erie, col = "#FFFFFF", border = "#000000", add = TRUE)
		plot(huron, col = "#FFFFFF", border = "#000000", add = TRUE)
		plot(michigan, col = "#FFFFFF", border = "#000000", add = TRUE)
		plot(ontario, col = "#FFFFFF", border = "#000000", add = TRUE)
		plot(stClair, col = "#FFFFFF", border = "#000000", add = TRUE)
		plot(superior, col = "#FFFFFF", border = "#000000", add = TRUE)
		dev.off()
	}

	print(paste0("species: ", species, " done"))
}
