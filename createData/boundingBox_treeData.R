
#### Aim of prog: Bounding box of the data, using projection EPSG:2163
library(data.table)
library(sf)

options(max.print = 500)
rm(list = ls())

#### Load tree data
treeData = readRDS("./treeData_cleaned.rds")
ls_14species = readRDS("./ls_species.rds")

treeData = unique(treeData[species_id %in% ls_14species, .(longitude, latitude)])

## Subset to Eastern America
treeData = treeData[longitude >= -100, ]

#### Transform the data
## Coerce to an sf object
treeData = st_as_sf(treeData, coords = c("longitude", "latitude"),
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

## Reproject to EPSG:2163
treeData = st_transform(treeData,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

## Bounding box
data_bbox = st_bbox(treeData)
saveRDS(data_bbox, "bbox_treeData.rds")
