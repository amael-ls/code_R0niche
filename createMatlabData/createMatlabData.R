
#### Aim of prog: Create the data matlab use to calculate R0 within Little's 1971 distribution
## Polygons from Little 1971: The data were downloaded from two websites (actually the second website is enough):
#		- https://databasin.org/datasets/
#		- https://www.fs.fed.us/nrs/atlas/littlefia/species_table.html
# I found the second website later, otherwise I would have used only that one.
#
## Prepare the data for Matlab:
#		- Climate in 2010 within the Little's 1971 distribution
#		- 4 data tables, for growth and mortality (above and below canopy for each demography)

library(data.table)
library(doParallel)
library(stringi)
library(raster)
library(velox)
library(sf)
library(sp)

rm(list = ls())
graphics.off()

#### Create the cluster
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

######## PART 1: Climate within the range of the species
#### Load climate data (raster stack version)
clim_2010 = raster::stack("../clim60sec/clim_2010.gri")
crs_clim = crs(clim_2010, asText = TRUE)

#### Create folders and load polygon (from Little 1971, )
## Get species and polygon folder
ls_folder = list.dirs("../little1971", recursive = FALSE)
ls_species = stri_sub(ls_folder, from = stri_locate_last(ls_folder, regex = "/")[,2] + 1)
(species = ls_species[array_id])

nb_species = length(ls_folder)
if (nb_species != length(ls_species))
	print("Number of species and number of folders mismatch")

## Matlab folder
matlabPath = "../R0/Matlab_data/"
if (!dir.exists("../R0"))
	dir.create("../R0")

if (!dir.exists(matlabPath))
	dir.create(matlabPath)

(savePath = paste0(matlabPath, species, "/"))
if (!dir.exists(savePath))
	dir.create(savePath)

## Polygon
littleEnvelop = st_read(dsn = ls_folder[array_id])
littleEnvelop = st_transform(x = littleEnvelop, crs = crs_clim)
littleEnvelop = st_geometry(littleEnvelop)

#### Create raster with ids as cell values (to keep track of which cell to extract later)
rid = setValues(clim_2010[[1]], 1:ncell(clim_2010[[1]]))

#### Extractions using velox
## Create velox rasters
v = velox(clim_2010)
vid = velox(rid)
print("Veloxes created")

## Extractions (small = TRUE, cf doc velox)
vals = v$extract(littleEnvelop, small = TRUE)
print("Extraction val done")

ids = vid$extract(littleEnvelop, small = TRUE)
print("Extraction Id done")

## Extract all cell centroids
xy = v$getCoordinates()
print(paste0("Number of polygons: ", length(vals)))

#### Associate xy points with values
xyvals = lapply(seq_along(vals), function(i)
	{
		print(paste0("i: ", i))
		print(paste0("length: ", nrow(vals[[i]])))
		if (!is.null(dim(vals[[i]])))
		{
			x = cbind(xy[ids[[i]][,1], , drop=FALSE], vals[[i]])
			x = as.data.table(x)
			setnames(x,c("x", "y", names(clim_2010)))
		}
	}
)

#### Save as a data table and raster
## Data table
database = rbindlist(xyvals)

# Reorder and rename columns
setcolorder(x = database, neworder = sort(names(database)))
climVar = readRDS("../createData/climaticVariables.rds")
setnames(database, old = names(database), new = c(climVar, "longitude", "latitude"))

print(dim(database))
spdf = data.table::copy(database)
database = na.omit(database)
print(dim(database))

saveRDS(database, paste0(savePath, species, "_clim2010.rds"))

## Raster
# Create spatial points data frame
coordinates(spdf) = ~ longitude + latitude

# Coerce to SpatialPixelsDataFrame
gridded(spdf) = TRUE

# Coerce to raster
rasterDF = stack(spdf)
writeRaster(rasterDF, filename = paste0(savePath, species, ".grd"), bandorder = "BIL", overwrite = TRUE)

######## PART 2: growth parameters
#### Tool function
matlabGrowthParams = function(parameters, temp, precip)
{
	beta0 = unname(parameters[, Intercept] +
		unlist(parameters[, "annual_mean_temperature"])*temp +
		unlist(parameters[, "I(annual_mean_temperature^2)"])*temp^2 +
		unlist(parameters[, "annual_precipitation"])*precip +
		unlist(parameters[, "I(annual_precipitation^2)"])*precip^2)

	beta1 = unname(parameters[, dbh] +
		unlist(parameters[, "dbh:annual_mean_temperature"])*temp +
		unlist(parameters[, "dbh:I(annual_mean_temperature^2)"])*temp^2 +
		unlist(parameters[, "dbh:annual_precipitation"])*precip +
		unlist(parameters[, "dbh:I(annual_precipitation^2)"])*precip^2)

	beta2 = unname(parameters[, dbh2] +
		unlist(parameters[, "I(dbh^2):annual_mean_temperature"])*temp +
		unlist(parameters[, "I(dbh^2):I(annual_mean_temperature^2)"])*temp^2 +
		unlist(parameters[, "I(dbh^2):annual_precipitation"])*precip +
		unlist(parameters[, "I(dbh^2):I(annual_precipitation^2)"])*precip^2)

	return (list(beta0 = beta0, beta1 = beta1, beta2 = beta2))
}

#### Read the climatic scaling
## Temperature
scaling_T = fread("./growthTempScaling.csv")

## Precipitations
scaling_P = fread("./growthPrecipScaling.csv")

#### Add growth parameters (for each location)
## Subset the climatic variables used in growth
matlabGrowth_above_dt = database[, .(longitude, latitude, annual_mean_temperature, annual_precipitation)]
matlabGrowth_below_dt = database[, .(longitude, latitude, annual_mean_temperature, annual_precipitation)]

## Species-specific scaling
matlabGrowth_above_dt[, annual_mean_temperature :=
	(annual_mean_temperature - scaling_T[species_id == species, mu])/scaling_T[species_id == species, sd]]
matlabGrowth_above_dt[, annual_precipitation :=
	(annual_precipitation - scaling_P[species_id == species, mu])/scaling_P[species_id == species, sd]]

matlabGrowth_below_dt[, annual_mean_temperature :=
	(annual_mean_temperature - scaling_T[species_id == species, mu])/scaling_T[species_id == species, sd]]
matlabGrowth_below_dt[, annual_precipitation :=
	(annual_precipitation - scaling_P[species_id == species, mu])/scaling_P[species_id == species, sd]]

## Read the species-specific parameters
parameters_above = readRDS("./parameters_above_growth.rds")
parameters_below = readRDS("./parameters_below_growth.rds")

## Merge the slopes (beta0, 1, 2) to climate
matlabGrowth_above_dt[, paste0("beta", 0:2) := matlabGrowthParams(parameters_above[species_id == species],
	annual_mean_temperature, annual_precipitation)]

matlabGrowth_below_dt[, paste0("beta", 0:2) := matlabGrowthParams(parameters_below[species_id == species],
	annual_mean_temperature, annual_precipitation)]

#### Save the species-specific climate/slopes data table
fwrite(matlabGrowth_above_dt, paste0(savePath, "matlabGrowth_above.csv"))
fwrite(matlabGrowth_below_dt, paste0(savePath, "matlabGrowth_below.csv"))

######## PART 3: mortality parameters
#### Tool function
matlabMortalityParams = function(parameters, temp, precip)
{
	beta0 = unname(parameters[, Intercept] +
		unlist(parameters[, "min_temperature_of_coldest_month"])*temp +
		unlist(parameters[, "I(min_temperature_of_coldest_month^2)"])*temp^2 +
		unlist(parameters[, "precipitation_of_driest_quarter"])*precip +
		unlist(parameters[, "I(precipitation_of_driest_quarter^2)"])*precip^2)

	beta1 = parameters[, dbh]

	beta2 = parameters[, dbh2]

	return (list(beta0 = beta0, beta1 = beta1, beta2 = beta2))
}

#### Read the climatic scaling
## Temperature
scaling_T = fread("./mortalityTempScaling.csv")

## Precipitations
scaling_P = fread("./mortalityPrecipScaling.csv")

#### Add mortality parameters (for each location), climate 2010 already in memory
## Subset the climatic variables used in mortality
matlabMortality_above_dt =
	database[, .(longitude, latitude, min_temperature_of_coldest_month, precipitation_of_driest_quarter)]
matlabMortality_below_dt =
	database[, .(longitude, latitude, min_temperature_of_coldest_month, precipitation_of_driest_quarter)]

## Species-specific scaling
matlabMortality_above_dt[, min_temperature_of_coldest_month :=
	(min_temperature_of_coldest_month - scaling_T[species_id == species, mu])/scaling_T[species_id == species, sd]]
matlabMortality_above_dt[, precipitation_of_driest_quarter :=
	(precipitation_of_driest_quarter - scaling_P[species_id == species, mu])/scaling_P[species_id == species, sd]]

matlabMortality_below_dt[, min_temperature_of_coldest_month :=
	(min_temperature_of_coldest_month - scaling_T[species_id == species, mu])/scaling_T[species_id == species, sd]]
matlabMortality_below_dt[, precipitation_of_driest_quarter :=
	(precipitation_of_driest_quarter - scaling_P[species_id == species, mu])/scaling_P[species_id == species, sd]]

## Read the species-specific parameters
parameters_above = readRDS("./parameters_above_mortality.rds")
parameters_below = readRDS("./parameters_below_mortality.rds")

## Merge the slopes (beta0, 1, 2) to climate
matlabMortality_above_dt[, paste0("beta", 0:2) := matlabMortalityParams(parameters_above[species_id == species],
	min_temperature_of_coldest_month, precipitation_of_driest_quarter)]

matlabMortality_below_dt[, paste0("beta", 0:2) := matlabMortalityParams(parameters_below[species_id == species],
	min_temperature_of_coldest_month, precipitation_of_driest_quarter)]

#### Save the species-specific climate/slopes data table
fwrite(matlabMortality_above_dt, paste0(savePath, "matlabMortality_above.csv"))
fwrite(matlabMortality_below_dt, paste0(savePath, "matlabMortality_below.csv"))
