
#### Aims of prog: Create the data for Matlab to calculate R0 on the plots
# This data will be used by Matlab to calculate R0 on all the plots I have a sample
# Then, the random forest (already done) that have been trained on these data will
# make prediction of presence/absence for each species
# Finally a correlation between R0 and the random forest prediction will be done
#

#### Load packages
library(data.table)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Create folder
if (!dir.exists("./Matlab_data"))
	dir.create("./Matlab_data")

######## PART 1: climate data of the plots
#### Load data
## Coordinates of the plot
coords = readRDS("../createData/treeData_presence_absence.rds")[, .(longitude, latitude)]

## Climate database average of data over 5 years, from 2006 to 2010
database = readRDS("../createData/clim_2010.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")
setnames(database, old = names(database[, !c("longitude", "latitude")]), new = climaticVariables)

## Species list
ls_species = readRDS("../createData/ls_species.rds")
n = length(ls_species)

#### Subset climate to the longitude-latitude that have samples, lose of 1609 data
database = database[coords, on = c("latitude", "longitude"), nomatch = 0]

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
scaling_T = fread("../createMatlabData/growthTempScaling.csv")

## Precipitations
scaling_P = fread("../createMatlabData/growthPrecipScaling.csv")

for (i in 1:n)
{
	#### Create folder
	species = ls_species[i]
	savePath = paste0("./Matlab_data/", species, "/")
	if (!dir.exists(savePath))
		dir.create(savePath)

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
	parameters_above = readRDS("../createMatlabData/parameters_above_growth.rds")
	parameters_below = readRDS("../createMatlabData/parameters_below_growth.rds")

	## Merge the slopes (beta0, 1, 2) to climate
	matlabGrowth_above_dt[, paste0("beta", 0:2) := matlabGrowthParams(parameters_above[species_id == species],
		annual_mean_temperature, annual_precipitation)]

	matlabGrowth_below_dt[, paste0("beta", 0:2) := matlabGrowthParams(parameters_below[species_id == species],
		annual_mean_temperature, annual_precipitation)]

	#### Save the species-specific climate/slopes data table
	fwrite(matlabGrowth_above_dt, paste0(savePath, "matlabGrowth_above.csv"))
	fwrite(matlabGrowth_below_dt, paste0(savePath, "matlabGrowth_below.csv"))
}

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
scaling_T = fread("../createMatlabData/mortalityTempScaling.csv")

## Precipitations
scaling_P = fread("../createMatlabData/mortalityPrecipScaling.csv")

for (i in 1:n)
{
	#### Create folder
	species = ls_species[i]
	savePath = paste0("./Matlab_data/", species, "/")
	if (!dir.exists(savePath))
		print(paste0("*** ERROR: dir ", savePath, " does not exist"))

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
	parameters_above = readRDS("../createMatlabData/parameters_above_mortality.rds")
	parameters_below = readRDS("../createMatlabData/parameters_below_mortality.rds")

	## Merge the slopes (beta0, 1, 2) to climate
	matlabMortality_above_dt[, paste0("beta", 0:2) := matlabMortalityParams(parameters_above[species_id == species],
		min_temperature_of_coldest_month, precipitation_of_driest_quarter)]

	matlabMortality_below_dt[, paste0("beta", 0:2) := matlabMortalityParams(parameters_below[species_id == species],
		min_temperature_of_coldest_month, precipitation_of_driest_quarter)]

	#### Save the species-specific climate/slopes data table
	fwrite(matlabMortality_above_dt, paste0(savePath, "matlabMortality_above.csv"))
	fwrite(matlabMortality_below_dt, paste0(savePath, "matlabMortality_below.csv"))
}
