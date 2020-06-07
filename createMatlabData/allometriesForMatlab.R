
#### Aim of prog: Make the allometries accessible to matlab
#		- A data table for dbh_star and dbh_inf, the lower and upper bounds of the integral
#		- A data table for the allometries (crown--dbh, etc...)
#
#### Units are very messy:
## - Adams2007 (from where the formulae are from):
#		* Growth in cm/yr
# 		* dbh in cm
#		* Height allometry (named H in his paper) in m
## - Purves2008 (from where the values are from):
#		* Growth in cm/yr
# 		* dbh in cm
# 		* Height allometry (named α in his paper) in m
## - On my own
#		* Growth in mm
#		* dbh in mm
# It implies I had to convert growth or dbh in some analytical formulae for the allometries

#### Load packages and clear memory
library(data.table)
library(xtable)

rm(list = ls())
graphics.off()

options(max.print = 500)

#### Source files allometries
source("../toolFunctions.R")
source("../createData/parametersAllometries.R")

#### Tool function
sortingSpecies = function(ls_species)
{
	tsn_dt = data.table(species = character(length(ls_species)),
		tsn = character(length(ls_species))) # tsn = Taxonomic Serial Number

	tsn_dt[, species := stri_sub(str = ls_species,
		from = stri_locate_first(str = ls_species, regex = "-")[,1] + 1)]

	tsn_dt[, tsn := stri_sub(str = ls_species, from = 1,
		to = stri_locate_first(str = ls_species, regex = "-")[,1] - 1)]
	tsn_dt = tsn_dt[order(species),]
	return (tsn_dt)
}

#### Common variables
h_star10 = 10 # h* = 10m
# h_inf = 45 # the infinit height is 45 m (for the upper bound of the integral)
ls_species = dir("../R0/Matlab_data")
nbSpecies = length(ls_species)

## Asymptotic height and dbh in m and mm respectively (from Burns 1990, Silvics 1 and 2)
# h_inf is the max height in general, and h_inf_max is the max height ever recorded/in optimal condition
# dbh_inf and dbh_inf_max follow the same rules
# For Fagus Grandifolia, it seems there is a mistake in the units. I use the table p. 660
# For Acer rubrum, dbh_inf_max is the bole circumference => does not reflect dbh

h_inf_dt = data.table(species_id = ls_species,
	h_inf = as.integer(c(18, 23, 30, 20, 20, 46, 30, 26.8, 30.5, 21, 25, 27, 37, 15)),
	h_inf_max = as.integer(c(27, 35, 55, 27, 30, 48.2, 49, 49, 34.7, 30, 26, 38.1, 23.8, 30)),
	dbh_inf = as.integer(c(460, 610, 900, 230, 250, 1020, 1020, 510, 760, 300, 300, 760, 910, 600)),
	dbh_inf_max = as.integer(c(750, 1330, 1200, 460, 640, 1680, 2130, 1350, 1440, 760, 1140, 4950, 2090, 1800)),
	age_max = as.integer(c(200, 400, 275, 250, 80, 200, 400, 366, 300, 150, 200, 125, 350, 400)),
	page_info = as.integer(c(33, 499, 410, 451, 566, 982, 1246, 660, 302, 348, 1093, 170, 202, 1197)))

dbh_params = data.table(species_id = ls_species, dbh_star10 = numeric(nbSpecies),
	dbh_infinity = numeric(nbSpecies))

#### Data table for the integral's bound
for (i in 1:nbSpecies)
{
	currentSpecies = ls_species[i]

	# Canopy height h* = 10m, associated species-specific dbh* (in mm)
	dbh_star = heightToDbh(height = h_star10,
		a = purves2007_allometries[species == currentSpecies, a],
		b = purves2007_allometries[species == currentSpecies, b])*10 # Times 10 to get it in mm

	# Infinity dbh (from Burns 1990, silvics)
	h_inf = h_inf_dt[species_id == currentSpecies, h_inf]
	dbh_inf = heightToDbh(height = h_inf,
		a = purves2007_allometries[species == currentSpecies, a],
		b = purves2007_allometries[species == currentSpecies, b])*10 # Times 10 to get it in mm

	# Fill dbh_params
	dbh_params[i, c("dbh_star10", "dbh_infinity") := .(dbh_star, dbh_inf)]
}

dbh_params = dbh_params[h_inf_dt, on = "species_id"]

#### Allometries data table
## Species-specific integral bounds
write.csv(dbh_params, "./dbh_params.csv", row.names = FALSE)

## Species-specific allometries
write.csv(C0_C1, "./C0_C1.csv", row.names = FALSE)
write.csv(purves2007_allometries[species %in% ls_species],
	"./purves2007_allometries.csv", row.names = FALSE)

## Species list
write.csv(ls_species, "./ls_species.csv", row.names = FALSE)

## Tex files
ls_species_dt = sortingSpecies(ls_species)
ls_species_dt[, species_id := paste0(tsn, "-", species)]
dbh_params = dbh_params[ls_species_dt, on = "species_id"]
dbh_params_xt = xtable(dbh_params[, .(species, dbh_inf, age_max, page_info)])
print(dbh_params_xt, file = "./dbh_age_max.tex", include.rownames = FALSE, booktabs = TRUE)
