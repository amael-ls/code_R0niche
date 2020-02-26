
#### Plot species distribution from Little 1971
# All maps (Little's and IV) are in Albers equal area projection as follows
# (this would be useful only to ArcInfo users):
#
# Projection: ALBERS
# Units: METERS
# Spheroid: CLARKE1866
# Parameters:
# 1st standard parallel: 38 0 0.000
# 2nd standard parallel: 42 0 0.000
# central meridian: -82 0 0.000
# latitude of projection's origin: 40 0 0.000
# false easting (meters): 0.00000
# false northing (meters): 0.00000
#
# From this information, I do a homemade proj4string
# "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD27 +units=m +no_defs";

library(stringi)
library(sf)

#### Load data
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

## Keep geometry only
bounding_canada = st_bbox(canada)
canada = st_geometry(canada)

bounding_usa = st_bbox(usa)
usa = st_geometry(usa)

erie = st_geometry(erie)
huron = st_geometry(huron)
michigan = st_geometry(michigan)
ontario = st_geometry(ontario)
stClair = st_geometry(stClair)
superior = st_geometry(superior)

lonMin = min(bounding_canada["xmin"], bounding_usa["xmin"])
lonMax = max(bounding_canada["xmax"], bounding_usa["xmax"])

latMin = min(bounding_canada["ymin"], bounding_usa["ymin"])
latMax = max(bounding_canada["ymax"], bounding_usa["ymax"])

#### Plot Little's 1971 species range
ls_folders = list.dirs(path = ".", recursive = FALSE)

for (folder in ls_folders)
{
	if (sum(stri_detect(str = folder, regex = c("ACE-RUB", "ACE-SAC", "BET-PAP"))) != 0)
	{
		shapefile = st_read(dsn = folder,
			crs = "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD27 +units=m +no_defs")

		species = stri_sub(str = folder, from = stri_locate_first(str = folder, regex = "./")[2] + 1)
		st_write(obj = shapefile, dsn = folder, layer = paste0(species, ".shp"), delete_dsn = TRUE, driver = "ESRI Shapefile")
	} else {
		shapefile = st_read(dsn = folder)
	}
	shapefile = st_transform(x = shapefile,
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
	shapefile = st_geometry(shapefile)
	species = stri_sub(str = folder, from = stri_locate_first(str = folder, regex = "./")[2] + 1)

	jpeg(paste0(folder, "/", species, ".jpg"), width = 1000, height = 1000, quality = 100)
	plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
		xlab = "", ylab = "")
	plot(shapefile, col = rgb(43/255, 161/255, 72/255, 0.2), add = TRUE)
	plot(canada, col = NA, add = TRUE)
	plot(usa, col = NA, add = TRUE)
	plot(erie, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(huron, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(michigan, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(ontario, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(stClair, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	plot(superior, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
	dev.off()
}

# plot(0, pch = "", xlim = st_bbox(shapefiles)[c(1,3)], ylim = st_bbox(shapefiles)[c(2,4)], axes = FALSE,
# 		xlab = "", ylab = "")
