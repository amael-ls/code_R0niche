
#### Aim of prog: Determine the species-specific centroid of Little's 1971
#

#### Load libraries and clear memory
library(stringi)
library(sf)

rm(list = ls())
graphics.off()

#### Compute centroid for each species
ls_folders = dir(path = "./", pattern = "^[0-9]{4}")
for(folder in ls_folders)
{
	# Create centroid folder in folder
	if (!dir.exists(paste0(folder, "/centroid")))
		dir.create(paste0(folder, "/centroid"))

	# Centroid
	shapefile = st_read(dsn = folder)
	shapefile = st_geometry(st_union(shapefile))
	centro = st_centroid(shapefile)
	# centro = st_transform(x = centro,
		# crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # EPSG 4326
	st_write(obj = centro, dsn = paste0(folder, "/centroid"), layer = paste0(folder, "_centroid.shp"),
		driver = "ESRI Shapefile")
	centro_mat = st_coordinates(centro)
	saveRDS(centro_mat, paste0(folder, "/centroid/", folder, "_centroid.rds"))
}
