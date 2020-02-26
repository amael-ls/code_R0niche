
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

rm(list = ls())
graphics.off()

options(max.print = 500)

#### Source files allometries
source("../toolFunctions.R")
source("../createData/parametersAllometries.R")

#### Common variables
h_star10 = 10 # h* = 10m
h_inf = 45 # the infinit height is 45 m (for the upper bound of the integral)
ls_species = dir("../R0/Matlab_data")
nbSpecies = length(ls_species)

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

	# Infinity dbh (here 45 meters)
	dbh_inf = heightToDbh(height = h_inf,
		a = purves2007_allometries[species == currentSpecies, a],
		b = purves2007_allometries[species == currentSpecies, b])*10 # Times 10 to get it in mm

	# Fill dbh_params
	dbh_params[i, c("dbh_star10", "dbh_infinity") := .(dbh_star, dbh_inf)]
}

#### Allometries data table
## Species-specific integral bounds
write.csv(dbh_params, "./dbh_params.csv", row.names = FALSE)

## Species-specific allometries
write.csv(C0_C1, "./C0_C1.csv", row.names = FALSE)
write.csv(purves2007_allometries[species %in% ls_species],
	"./purves2007_allometries.csv", row.names = FALSE)

## Species list
write.csv(ls_species, "./ls_species.csv", row.names = FALSE)
