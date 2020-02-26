
#### Aim of prog: Download climatic data from ANUSPLIN
## Create cluster
#		- One year per array
#
## Download data
#		Define url, year, resolution
#		Create saving folder
#
## Remark:
# It might be necessary to re-run this program several times to get all the data.
# Indeed, the server sometimes does not respond (maybe due to the parallelisation)
# The program detect automatically which years are missing hence, when you re-run it,
# the program only tries the missing years.
#
# I also got a problem with unzip (due to a lack of memory), hence I ceck everything
# is unzipped
#
## Metadata 60 arc-sec
# ID DATA:  northam.52j.latest.FLs.52j.tar
# ORIGINATOR: LANDSCAPE ANALYSIS AND APPLICATION SECTION (LAAS)
#   GREAT LAKES FORESTRY CENTRE (GLFC), CANADIAN FOREST SERVICE (CFS),
#   NATURAL RESOURCES CANADA (NRCAN)
# ID PARAMETER:  MAXT, MINT, PCP, CMI, BIO, SGRO
# TIME INTERVAL: MONTHLY
# PERIOD:        1900-2015
#
# SPATIAL FORMAT: GRID
# DATA FORMAT: ARC/Info ASCIIGRID
# GRIDS DATUM: GEOGRAPHIC NAD83
# SPATIAL RESOLUTION: 60 arc-second (~2 km)
# GRIDS DOMAIN: 168.0 to 52.0W, 25.0 to 85.0N
# ESRI ASCII RASTER FORMAT:
# NCOLS   6960
# NROWS   3600
# XLLCORNER  -168.0
# YLLCORNER    25.0
# CELLSIZE 0.016666667536
# NODATA VALUE -99
#
# UNITS:
# MAXT/MINT (deg C)
# PCP (mm)
# CMI (cm)
# BIO, SGRO: The bioclim/seedgro parameters list and definitions are at: http://cfs.nrcan.gc.ca/projects/3/6


#### Load packages and clear memory
rm(list = ls())

#### Create the cluster
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

#### Get the data
## Common variables
allYear = 1956:2013 # Can be from 1900 to 2015
(year = allYear[array_id])

# url address Dan McKenney (ftp://ftp.nrcan.gc.ca/pub/outgoing/NAM_grids)
basurl = "ftp://ftp.nrcan.gc.ca/pub/outgoing/NAM_grids/zipfiles60/"

# Variable to download, among "bio", "cmi", "mint", "maxt", "pcp", "sg"
info = "bio"

# Resolution either "_300arcsec.zip" or "_60arcsec.zip", 300 = 10km², 60 = 2km²
end = "_60arcsec.zip"

## Download data
# File name & create folder
(zout = paste0("./climateData/", info, year, ".zip"))

if (!dir.exists("./climateData"))
	dir.create("./climateData")

if (!dir.exists("./climateData/bioclim"))
	dir.create("./climateData/bioclim")

if (!dir.exists(paste0("./climateData/bioclim/", year)))
{
	download.file(paste0(basurl, info, year, end),
		destfile = zout, method = "wget")

	## Unzip in the folder climateData
	unzip(zout, exdir = "./climateData/bioclim/")
	ls_asc = list.files(path = paste0("./climateData/bioclim/", year), pattern = ".asc")
	if (length(ls_asc) != 19)
		print(paste0("WARNING: unzip function had a problem, year = ", year))
} else {
	ls_asc = list.files(path = paste0("./climateData/bioclim/", year), pattern = ".asc")
	if (length(ls_asc) != 19)
	{
		## Unzip in the folder climateData
		unzip(zout, exdir = "./climateData/bioclim/")
		ls_asc = list.files(path = paste0("./climateData/bioclim/", year), pattern = ".asc")
		if (length(ls_asc) != 19)
			print(paste0("WARNING: unzip function had a problem, year = ", year))
	}
}

#### Quick and dirty
# for (year in 1956:2013)
# {
# 	ls_asc = list.files(path = paste0("./climateData/bioclim/", year), pattern = ".asc")
# 	if (length(ls_asc) != 19)
# 	{
# 		print(paste0("WARNING: unzip function had a problem, year = ", year))
# 		unzip(zout, exdir = "./climateData/bioclim/")
# 		if (length(ls_asc) != 19)
# 			print(paste0("WARNING: unzip function had a problem, year = ", year))
#
# 	}
# }
