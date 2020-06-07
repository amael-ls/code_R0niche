
library(sf)

ls_species = dir(path = "../little1971/", pattern = "^[0-9]{1,}")

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

data_bbox = readRDS("../createData/bbox_treeData.rds")
st_crs(data_bbox) == st_crs(canada)
st_crs(data_bbox) == st_crs(usa)

jpeg("plotCentroid.jpg", height = 1080, width = 1080, quality = 100)
plot(canada, col = NA, lwd = 2)
plot(usa, col = NA, add = TRUE, lwd = 2)

for (sp in ls_species)
{
	centre = st_read(dsn = paste0("../little1971/", sp, "/centroid"))
	centre = st_transform(centre,
		crs = st_crs(usa))
	plot(centre, lwd = 4, cex = 4, col = "red", add = TRUE)
	coords = st_coordinates(centre)
	text(x = coords[1, "X"], y = coords[1, "Y"], labels = sp)
}
dev.off()
