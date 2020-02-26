
#### Aim of prog: Create a single climatic database for the trees
# This database will later be merged to tree data

#### Load library and clear memory
library(data.table)
library(stringi)

rm(list = ls())

#### Load data
## List files
(loadPath = "./climateData/averageClim/")
ls_files = list.files(loadPath)

## Create list
ls_clim = vector(mode = "list", length = length(ls_files))

for (i in 1:length(ls_files))
	ls_clim[[i]] = readRDS(paste0(loadPath, ls_files[i]))

#### Create database and save in a single file
## Database
climForTrees = rbindlist(ls_clim)

## Set the new names (Informations can be found at: http://cfs.nrcan.gc.ca/projects/3/8)
# For more informations, check program subsetClimForTrees
newNames = c(
	"annual_mean_temperature",
	"mean_diurnal_range",
	"isothermality",
	"temperature_seasonality",
	"max_temperature_of_warmest_month",
	"min_temperature_of_coldest_month",
	"temperature_annual_range",
	"mean_temperature_of_wettest_quarter",
	"mean_temperature_of_driest_quarter",
	"mean_temperature_of_warmest_quarter",
	"mean_temperature_of_coldest_quarter",
	"annual_precipitation",
	"precipitation_of_wettest_month",
	"precipitation_of_driest_month",
	"precipitation_seasonality",
	"precipitation_of_wettest_quarter",
	"precipitation_of_driest_quarter",
	"precipitation_of_warmest_quarter",
	"precipitation_of_coldest_quarter")

currentNames = names(climForTrees)[stri_detect(str = names(climForTrees), regex = "bio60_")]

setnames(climForTrees, old = currentNames, new = newNames)

## Remove the NAs
clim[, lapply(.SD, function (x) {sum(is.na(x))})]
sum(is.na(climForTrees))
climForTrees = na.omit(climForTrees)
sum(is.na(climForTrees))

## Save
saveRDS(climForTrees, "../createData/averagedClim_5yearsAllVar.rds")
saveRDS(newNames, "../createData/climaticVariables.rds")
