
#### Aim of prog: Extract climate for the whole Northern America, year = 2010
# This climate database will be used by the random forest.
# I chose to averaged the climatic data over 5 years, i.e., 2006 - 2010

#### Load library and clear memory
library(data.table)
library(stringi)
library(raster)
library(sp)

options(max.print = 500)
rm(list = ls())

#### Load data and average them
## Common variables
ls_years = 2006:2010
final_raster = stack()
count = 1

## Define variables that must be divided by 10 (cf 1st remark at the begening of the program)
varToDevideBy10 = paste0("bio60_", c("01", "02", "05", "06", "07", "08", "09", "10", "11"))
ls_17biovar = paste0("bio60_", c(paste0("0", 1:9), 10:12, 15:19))

## Work on the 17 biovar
for (j in 1:17)
{
	raster_stack_bio = stack()
	# Read for the 5 years 2006 - 2010
	for (year in ls_years)
	{
		loadPath_bio = paste0("climateData/bioclim/", year, "/")

		# List the bioclim files, remove suspicious bioclim BIO13 and BIO14 (cf remark 2 in subsetClimForTrees.R)
		raster_files_bio = paste0(loadPath_bio, list.files(path = loadPath_bio, pattern = ".asc"))
		raster_files_bio = raster_files_bio[!stri_detect_regex(str = raster_files_bio, pattern = "13|14")]

		rs = raster(raster_files_bio[j])
		raster_stack_bio = stack(raster_stack_bio, rs)
	}
	# Average over the 5 years
	rasterMean = calc(raster_stack_bio, fun = mean)
	names(rasterMean) = ls_17biovar[j]

	# Correction of factor 10 for some variables (cf subsetClimForTrees.R)
	if (sum(stri_detect(raster_files_bio[j], regex = varToDevideBy10)) != 0)
	{
		rasterMean = calc(rasterMean, function(x) {x/10})
		names(rasterMean) = varToDevideBy10[count]
		count = count + 1
	}
	final_raster = stack(final_raster, rasterMean)
	print(paste0("biovar: ", j, " done"))
}

crs_bioclim = crs(final_raster)

#### Stack the required rasters of precipitation, to create bioclim 13 & 14
bio60_13 = stack()
bio60_14 = stack()
for (year in ls_years)
{
	raster_stack_pcp = stack()
	loadPath_pcp = paste0("climateData/pcp/", year, "/")
	raster_files_pcp = paste0(loadPath_pcp, list.files(path = loadPath_pcp, pattern = ".asc"))

	# Load climate data (12 months)
	for (j in 1:length(raster_files_pcp))
		raster_stack_pcp = stack(raster_stack_pcp, raster(raster_files_pcp[j]))

	bio60_13_yr = calc(raster_stack_pcp, function(x) {max(x)})
	bio60_13 = stack(bio60_13, bio60_13_yr)

	bio60_14_yr = calc(raster_stack_pcp, function(x) {min(x)})
	bio60_14 = stack(bio60_14, bio60_14_yr)
}

bio60_13 = calc(bio60_13, fun = mean)
bio60_14 = calc(bio60_14, fun = mean)

names(bio60_13) = "bio60_13"
names(bio60_14) = "bio60_14"

final_raster = stack(final_raster, bio60_13, bio60_14)
final_raster = projectRaster(final_raster,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

#### Save the stacks in different format
## Raster stack format (really slow, so only if the file does not exist)
if (!file.exists("clim_2010.grd"))
	writeRaster(final_raster, filename = "clim_2010.grd", bandorder = "BIL", overwrite = TRUE)

## Data table format
# Convert to sp data frame
final_raster = rasterToPoints(x = final_raster, spatial = TRUE)
if (!file.exists("clim2010sp/bio60.dbf"))
	rgdal::writeOGR(obj = final_raster, dsn = "clim2010sp", "bio60", driver = "ESRI Shapefile")

# Convert to data table
final_raster = as.data.frame(final_raster)
setDT(final_raster)

# Rename and reorder
setnames(x = final_raster, old = c("x", "y"), new = c("longitude", "latitude"))

biovar = sort(colnames(final_raster)[!(colnames(final_raster) %in% c("latitude", "longitude"))])
setcolorder(x = final_raster, neworder = c("latitude", "longitude", biovar))

final_raster = na.omit(final_raster) # Loss of 107342 values

saveRDS(final_raster, "../createData/clim_2010.rds")
