
#### Aim of prog: Compute the difference of presence and absence in function of competition
library(data.table)
library(stringi)
library(raster)
library(sf)

options(max.print = 500)
rm(list = ls())

#### Tool functions
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
ls_folder = list.files(path = "../R0/results/", pattern = "^[0-9]{4,}")

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
data_bbox = readRDS("../createData/bbox_treeData.rds")
data_extent = extent(matrix(data_bbox, nrow = 2, ncol = 2))

#### Run for all the species
## Plot maps
for (i in 1:length(ls_folder))
{
	## Load species-specific data
	species = ls_folder[i]
	folder = ls_folder[i]

	# R0 and climate (do not reorganise the rows before join)
	R0_vec_10m = fread(paste0("../R0/results/", folder, "/R0_10m.csv"))
	R0_vec_0m = fread(paste0("../R0/results/", folder, "/R0_0m.csv"))
	clim_2010 = readRDS(paste0("../R0/Matlab_data/", folder, "/", species, "_clim2010.rds"))

	clim_2010 = clim_2010[, .(longitude, latitude)]
	clim_2010[, R0_10m := R0_vec_10m]
	clim_2010[, R0_0m := R0_vec_0m]

	## Compute presence/absence
	clim_2010[, presence_10m := ifelse(R0_vec_10m >= 1, 1, 0)] # 1 = presence, 0 = absence
	clim_2010[, presence_0m := ifelse(R0_vec_0m >= 1, 1, 0)] # 1 = presence, 0 = absence

	## Compute change when adding competition (color plot):
	# -1 = change from present to absent when adding competition (#2A3344)
	# 0 = no change in absence (#2058DC)
	# 1 = no change in presence (#FAB935)
	# 2 = change from absent to present when adding competition (#B6343A)
	clim_2010[, diffPresAbs_addedCompetition := 2*presence_10m - presence_0m]

	# Load the Little's 1971 shapefile, crop to data bbox
	ls_little = list.files(path = paste0("../little1971/", folder), pattern = ".shp$")
	little1971 = st_read(paste0("../little1971/", folder, "/", ls_little))
	little1971 = st_transform(little1971,
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
	little1971 = st_crop(little1971, data_bbox)
	little1971 = st_geometry(little1971)
	little1971_sp = as(little1971, "Spatial")

	## Set clim_2010 as raster
	clim_2010 = rasterFromXYZ(clim_2010[, .(longitude, latitude, diffPresAbs_addedCompetition)],
		crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

	if (names(clim_2010) != "diffPresAbs_addedCompetition")
	{
		print("For whaterver reason, R changed the name")
		names(clim_2010) = c("diffPresAbs_addedCompetition")
	}

	clim_2010 = crop(clim_2010, data_extent)

	## Fill the NA within little1971
	clim_2010 = fillGap(clim_2010, little1971_sp)

	## Detect borders
	borderlines = boundaries(clim_2010, classes = TRUE, asNA = TRUE, filename = paste0("./diff_presAbs/", species, ".grd"), overwrite = TRUE)
	print(paste0("species: ", species, " done"))
}

# Note that for the article I ran again the map_grad.R function from R0/maps/ and added manually the raster borderlines on top
