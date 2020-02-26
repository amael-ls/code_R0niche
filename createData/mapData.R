
#### Aim of prog: Map the data with EPSG:4269
library(data.table)
library(sf)

options(max.print = 500)
rm(list = ls())

#### Load data
## Shapefiles Northern America
canada = readRDS("~/database/shapefiles/North_America_Shapefile/canada_union.rds")
usa = readRDS("~/database/shapefiles/North_America_Shapefile/usa_union.rds")

erie = st_read("~/database/shapefiles/North_America_Shapefile/lakes/erie")
huron = st_read("~/database/shapefiles/North_America_Shapefile/lakes/huron")
michigan = st_read("~/database/shapefiles/North_America_Shapefile/lakes/michigan")
ontario = st_read("~/database/shapefiles/North_America_Shapefile/lakes/ontario")
stClair = st_read("~/database/shapefiles/North_America_Shapefile/lakes/stClair")
superior = st_read("~/database/shapefiles/North_America_Shapefile/lakes/superior")

canada = st_geometry(canada)
usa = st_geometry(usa)
erie = st_geometry(erie)
huron = st_geometry(huron)
michigan = st_geometry(michigan)
ontario = st_geometry(ontario)
stClair = st_geometry(stClair)
superior = st_geometry(superior)

## Tree data
treeData = readRDS("./treeData_cleaned.rds")
ls_14species = readRDS("./ls_species.rds")

treeData = unique(treeData[species_id %in% ls_14species, .(longitude, latitude)])

## Subset to Eastern America
treeData = treeData[longitude >= -100, ]

#### Transform the data
## Coerce to an sf object
treeData = st_as_sf(treeData, coords = c("longitude", "latitude"),
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

## Reproject
treeData = st_transform(treeData, crs = st_crs(usa))

## Bounding box
data_bbox = st_bbox(treeData)

xmin = min(st_bbox(usa)["xmin"], st_bbox(canada)["xmin"])
xmax = max(st_bbox(usa)["xmax"], st_bbox(canada)["xmax"])
ymin = min(st_bbox(usa)["ymin"], st_bbox(canada)["ymin"])
ymax = max(st_bbox(usa)["ymax"], st_bbox(canada)["ymax"])

coords = st_coordinates(treeData)

#### Plot
## pdf
pdf("./mapDatabase.pdf", height = 8, width = 8)
	plot(0, pch = "", xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE,
		xlab = "", ylab = "", bg = "transparent")
	points(coords, pch = 20, cex = 0.5, col = "#434343")
	plot(canada, col = NA, add = TRUE)
	plot(usa, col = NA, add = TRUE)
	plot(erie, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(huron, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(michigan, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(ontario, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(stClair, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(superior, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(st_as_sfc(data_bbox), lwd = 6, col = NA, border = "#ff9933", add = TRUE)
dev.off()

## jpeg
jpeg("./mapDatabase.jpg", height = 1000, width = 1000, quality = 100)
	plot(0, pch = "", xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE,
		xlab = "", ylab = "", bg = "transparent")
	points(coords, pch = 20, cex = 0.5, col = "#434343")
	plot(canada, col = NA, add = TRUE)
	plot(usa, col = NA, add = TRUE)
	plot(erie, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(huron, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(michigan, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(ontario, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(stClair, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(superior, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(st_as_sfc(data_bbox), lwd = 6, col = NA, border = "#ff9933", add = TRUE)
dev.off()
