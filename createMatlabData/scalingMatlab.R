
#### Aim of prog: Create the csv files containing the scaling information for Matlab
## Growth part:
#		- annual mean temperature
#		- annual precipitation
#
## Mortality part:
#		- min temperature of coldest month
#		- precipitation of driest quarter
#

#### Load packages and clear memory
library(data.table)
library(stringi)

rm(list = ls())
graphics.off()

#### Growth
ls_folder = list.files(path = "../growth", pattern = "array_[0-9]")
scaling_G = vector(mode = "list", length = length(ls_folder))
scaling_DBH = vector(mode = "list", length = length(ls_folder))
scaling_T = vector(mode = "list", length = length(ls_folder))
scaling_P = vector(mode = "list", length = length(ls_folder))
sp_vector = character(length = length(ls_folder))

for (i in 1:length(ls_folder))
{
	(loadPath = paste0("../growth/", ls_folder[i], "/"))
	scaling = readRDS(paste0(loadPath, "normalisation_growth_data_log_final.rds"))
	scaling_G[[i]] = scaling[var == "growth",]
	scaling_DBH[[i]] = scaling[var == "dbh",]
	scaling_T[[i]] = scaling[var == "annual_mean_temperature",]
	scaling_P[[i]] = scaling[var == "annual_precipitation",]
	species = list.files(path = loadPath, pattern = "^[0-9]{4,}")
	species = species[stri_detect(str = species, regex = ".txt")]
	sp_vector[i] = stri_sub(str = species, to = stri_locate_last(str = species, regex = ".txt")[1] - 1)
}

scaling_G = rbindlist(scaling_G)
scaling_G[, species_id := sp_vector]

scaling_DBH = rbindlist(scaling_DBH)
scaling_DBH[, species_id := sp_vector]

scaling_T = rbindlist(scaling_T)
scaling_T[, species_id := sp_vector]

scaling_P = rbindlist(scaling_P)
scaling_P[, species_id := sp_vector]

fwrite(scaling_G, "./growthScaling.csv")
fwrite(scaling_DBH, "./growthDbhScaling.csv")
fwrite(scaling_T, "./growthTempScaling.csv")
fwrite(scaling_P, "./growthPrecipScaling.csv")

#### Mortality
ls_folder = list.files(path = "../mortality", pattern = "array_[0-9]")
scaling_DBH = vector(mode = "list", length = length(ls_folder))
scaling_T = vector(mode = "list", length = length(ls_folder))
scaling_P = vector(mode = "list", length = length(ls_folder))
sp_vector = character(length = length(ls_folder))

for (i in 1:length(ls_folder))
{
	(loadPath = paste0("../mortality/", ls_folder[i], "/"))
	scaling = readRDS(paste0(loadPath, "normalisation_mortality_data.rds"))
	scaling_DBH[[i]] = scaling[var == "dbh",]
	scaling_T[[i]] = scaling[var == "min_temperature_of_coldest_month",]
	scaling_P[[i]] = scaling[var == "precipitation_of_driest_quarter",]
	species = list.files(path = loadPath, pattern = "^[0-9]{4,}")
	species = species[stri_detect(str = species, regex = ".txt")]
	sp_vector[i] = stri_sub(str = species, to = stri_locate_last(str = species, regex = ".txt")[1] - 1)
}

scaling_DBH = rbindlist(scaling_DBH)
scaling_DBH[, species_id := sp_vector]

scaling_T = rbindlist(scaling_T)
scaling_T[, species_id := sp_vector]

scaling_P = rbindlist(scaling_P)
scaling_P[, species_id := sp_vector]

fwrite(scaling_DBH, "./mortalityDbhScaling.csv")
fwrite(scaling_T, "./mortalityTempScaling.csv")
fwrite(scaling_P, "./mortalityPrecipScaling.csv")
