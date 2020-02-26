
#### Aim of prog: Create climate data for trees
## Read data (trees)
#		- Tree data
#		- Extract year with array_id (I checked in advance how many I needed)
#
## Read data (climate)
#		- The climatic data are loaded over a time span defined by lengthAverage
#
## Average climate over 5 years
#		- Create a single data table
#		- Correction of some variables (dividing by 10, cf remarks)
#		- Apply mean to each biovar col
#
#### Remarks
## Remark 1
# There are 19 bioclimatic variables, named bio60_01, ..., bio60_19.
# Further information about them can be found at:
#		0/ https://pubs.usgs.gov/ds/691/ds691.pdf
#
#		1/ https://fennerschool.anu.edu.au/files/anuclim61.pdf
#
# *** Section 6.1, page 61
# The  descriptions  below  assume that BIOCLIM is using the weekly time step (the default).
# If the monthly time step is selected, monthly values are used when calculating these parameters.
# The quarterly parameters are not aligned to any calendar quarters. BIOCLIM's definition of a
# quarter is any 13 consecutive weeks, (or any consecutive 3 months if running with a monthly
# time step). For example, the driest quarter will be the 13 consecutive weeks that are drier than
# any other set of 13 consecutive weeks.
#
# *** /!\ Section 1.4, page 10 /!\
# When outputting results from ESOCLIM, BIOCLIMand GROCLIMof AN UCLIMVersion 5.1
# to Arc/InfoUNGENERATE files (point data), Arc/InfoASCIIGRID or IDRISI ASCII image files,
# the output values were multiplied by 10 or 100 and then rounded to the nearest integer.
# This was done to preserve appropriate precision in the output data with smaller output files.
# However this process led to some confusion and inconvenience, particularly when these outputs
# were compared with data from other sources.
#
#
#		2/ http://worldclim.org/bioclim
#
# BIO1  = Annual Mean Temperature
# BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3  = Isothermality (BIO2/BIO7)
# BIO4  = Temperature Seasonality (Coefficient of Variation)
# BIO5  = Max Temperature of Warmest Month
# BIO6  = Min Temperature of Coldest Month
# BIO7  = Temperature Annual Range (BIO5-BIO6)
# BIO8  = Mean Temperature of Wettest Quarter
# BIO9  = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
#
## Remark 2
# BIO13 is suspicious (at least for the year 1974).
# For instance, there is BIO12 = 1048mm, and BIO13 = 46mm for the same coordinate.
# This does not make sense, because if we assume all the months are the same, then
# I should get BIO12 < 12 x BIO13 (with equality if and only if my assumption is true).
# Same, if I use BIO16, I should get:
#	BIO12 < 4 x BIO16, which is the case: sum(climateData$bio60_12 > 4*climateData$bio60_16) == 0
#	BIO16 < 4 x BIO13, which is NOT the case...
# Hence, I decided to rebuild BIO13 from the raw data
# Similarly I rebuild BIO14.
# Actually, it seems that BIO13 is BIO14 from the raw data, maybe they made a mistake
# between min and max. However, what they have done from BIO14 is a mistery...

#### Load library and clear memory
library(matrixStats)
library(data.table)
library(doParallel)
library(stringi)
library(raster)

options(max.print = 500)
rm(list = ls())

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

#### Common variables
## Load data
treeData = readRDS("../createData/treeData_cleaned.rds")
lengthAverage = 5
ls_years = sort(treeData[, unique(year_measured)])

## Define variables that must be divided by 10 (cf 1st remark at the begening of the program)
varToDevideBy10 = paste0("bio60_", c("01", "02", "05", "06", "07", "08", "09", "10", "11"))

## Get current year and the associated coordinates
if (array_id > length(ls_years))
	print(paste0("*** ERROR: array_id must be smaller than ", length(ls_years)))

(year = ls_years[array_id])

coords = unique(treeData[year_measured == year, .(longitude, latitude)])
print(paste0("Number of coordinates to extract: ", coords[, .N]))

## List of climatic data
ls_climData = vector(mode = "list", length = lengthAverage)
count = 1

#### Stack the required rasters ((19 - 2) x lengthAverage), but bioclim 13 & 14
for (yr in (year - lengthAverage + 1):year)
{
	raster_stack = stack()

	loadPath = paste0("./climateData/bioclim/", yr, "/")
	raster_files = paste0(loadPath, list.files(path = loadPath, pattern = ".asc"))

	# Remove suspicious bioclim BIO13 and BIO14 (cf remark 2)
	raster_files = raster_files[!stri_detect_regex(str = raster_files, pattern = "13|14")]

	# Load climate data (19 biovar)
	for (j in 1:length(raster_files))
		raster_stack = stack(raster_stack, raster(raster_files[j]))

	ls_climData[[count]] = setDT(extract(x = raster_stack, y = coords,
		method = "simple", df = TRUE))
	count = count + 1
	print(paste0("current year ", yr, " loaded for the year ", year))
}

#### Average the data of all the bioclim but 13 & 14
## List of biovar
biovar = stri_sub(raster_files,
	from = stri_locate_last(str = raster_files, regex = "/")[1] + 1,
	to = stri_locate_last(str = raster_files, regex = ".asc")[1] - 1)

## Merge all the data tables, correction of factor 10
climateData = rbindlist(ls_climData)
climateData[, c(varToDevideBy10) := lapply(.SD, function (x) {return (x/10)}), .SDcols = varToDevideBy10]

## Averaging of lengthAverage years
averageClim = climateData[, lapply(.SD, mean), by = ID, .SDcols = c(biovar)]
averageClim[, ID := NULL]

#### Stack the required rasters of bioclim 13 & 14 (2 x lengthAverage)
ls_climData_pcp = vector(mode = "list", length = lengthAverage)
count = 1
for (yr in (year - lengthAverage + 1):year)
{
	raster_stack_pcp = stack()

	loadPath = paste0("./climateData/pcp/", yr, "/")
	raster_files_pcp = paste0(loadPath, list.files(path = loadPath, pattern = ".asc"))

	# Load climate data (12 months)
	for (j in 1:length(raster_files_pcp))
		raster_stack_pcp = stack(raster_stack_pcp, raster(raster_files_pcp[j]))

	ls_climData_pcp[[count]] = setDT(extract(x = raster_stack_pcp, y = coords,
		method = "simple", df = TRUE))

	## Calculate bioclim 13 & 14
	ls_climData_pcp[[count]][, bio60_13 := rowMaxs(as.matrix(.SD), na.rm = TRUE), .SDcols = !"ID"]
	ls_climData_pcp[[count]][, bio60_14 := rowMins(as.matrix(.SD), na.rm = TRUE), .SDcols = !"ID"]

	count = count + 1
	print(paste0("current year ", yr, " loaded for the year ", year))
}

#### Average the data of bioclim 13 & 14
## List of pcp/biovar
pcpvar = stri_sub(raster_files_pcp,
	from = stri_locate_last(str = raster_files_pcp, regex = "/")[1] + 1,
	to = stri_locate_last(str = raster_files_pcp, regex = ".asc")[1] - 1)

biovar13_14 = c("bio60_13", "bio60_14")

## Merge all the data tables, correction of factor 10
climateData_pcp = rbindlist(ls_climData_pcp)
climateData_pcp[, c(pcpvar) := NULL]

## Averaging of lengthAverage years
averageClim_pcp = climateData_pcp[, lapply(.SD, mean), by = ID, .SDcols = c(biovar13_14)]
averageClim_pcp[, ID := NULL]

#### Saving the results
## Add longitude, latitude, and year
averageClim[, c("longitude", "latitude") := coords]
averageClim[, year_measured := year]

## Add bioclim 13 & 14
averageClim[, c(biovar13_14) := averageClim_pcp]

biovar = sort(colnames(averageClim)[!(colnames(averageClim) %in% c("latitude", "longitude", "year_measured"))])
setcolorder(x = averageClim,
	neworder = c("latitude", "longitude", "year_measured", biovar))

## Save in ./climateData/averageClim/
if (!dir.exists("./climateData/averageClim/"))
	dir.create("./climateData/averageClim/")

saveRDS(averageClim, paste0("./climateData/averageClim/", year, ".rds"))
