
#### Aim of prog: Make maps from matlab results
library(RColorBrewer)
library(data.table)
library(stringi)
library(scales)
library(sf)

options(max.print = 500)
rm(list = ls())

#### Common variables
## List folders
ls_folder = list.files(path = "../results/", pattern = "^[0-9]{4,}")

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

## Load bounding box determined by the data. Willbe used to crop the results
data_bbox = readRDS("../../createData/bbox_treeData.rds")
st_crs(data_bbox) == st_crs(canada)
st_crs(data_bbox) == st_crs(usa)

#### Run for all the species
for (i in 1:length(ls_folder))
{
	## Load species-specific data
	species = ls_folder[i]
	folder = ls_folder[i]

	# R0 and climate (do not reorganise the rows before join)
	R0_10m = fread(paste0("../results/", folder, "/R0_10m.csv"))
	clim_2010 = readRDS(paste0("../Matlab_data/", folder, "/", species, "_clim2010.rds"))

	clim_2010 = clim_2010[, .(longitude, latitude)]
	clim_2010[, R0 := R0_10m]

	# Load the Little's 1971 shapefile, crop to data bbox
	ls_little = list.files(path = paste0("../../little1971/", folder), pattern = ".shp$")
	little1971 = st_read(paste0("../../little1971/", folder, "/", ls_little))
	little1971 = st_transform(little1971,
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
	little1971 = st_crop(little1971, data_bbox)
	little1971 = st_geometry(little1971)

	## Set clim_2010 as sf
	clim_2010 = st_as_sf(x = clim_2010, coords = c("longitude", "latitude"),
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

	# Crop clim_2010 to the data bounding box
	clim_2010 = st_crop(clim_2010, data_bbox)

	# Calculate rho_0
	clim_2010$R0_scaled = (clim_2010$R0 - min(clim_2010$R0))/(max(clim_2010$R0) - min(clim_2010$R0))

	# Save lon-lat, R0, rho_0 sp-specific data
	saveRDS(clim_2010, paste0("../results/", folder, "/lonLatR0cropped.rds"))

	## Bounding box of plot
	lonMin = unname(data_bbox["xmin"]) - 100
	lonMax = unname(data_bbox["xmax"]) + 100
	latMin = unname(data_bbox["ymin"]) - 100
	latMax = unname(data_bbox["ymax"]) + 100

	## Spatial plot of R0(h* = 10m)
	# Format = jpeg
	jpeg(paste0("./", species, "_R0_h=10m_scaled.jpg"), width = 1000, height = 1000, quality = 100)
	plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
		xlab = "", ylab = "", bg = "transparent")
	plot(clim_2010["R0_scaled"], cex = 0.5, pch = 19, add = TRUE, # , "#4269E2"
		pal = c("#000000", "#2A3344", "#1122AA", "#2058DC", "#5A82EA", "#FDECBE", "#FDDB5B", "#FAB935", "#FD9859", "#B6343A"))
	plot(little1971, col = NA, border = "#000000", add = TRUE, lwd = 2) # "#DD925C"
	plot(canada, col = NA, add = TRUE, lwd = 2)
	plot(usa, col = NA, add = TRUE, lwd = 2)
	plot(erie, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(huron, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(michigan, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(ontario, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(stClair, col = "#FFFFFF", border = "#000000", add = TRUE)
	plot(superior, col = "#FFFFFF", border = "#000000", add = TRUE)
	dev.off()

	print(paste0("species: ", species, " done"))
}
