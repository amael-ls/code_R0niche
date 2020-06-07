
#### Aim of prog: plot the mortality in function of dbh, with the range of dbh parametrisation
## Rebuild the parameters:
#		- Read (G)LMM results
#		- Recombine for each species using stringi in params names
#
## Plot
#		- Using tikzDevice for latex
#

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool functions
## Combine names to make the interaction name automatically
combineNames = function(str1, str2)
{
	n1 = length(str1)
	n2 = length(str2)
	comb = character(length = n1*n2)
	for (i in 1:n1)
		for (j in 1:n2)
			comb[(i - 1)*n2 + j] = paste0(str1[i], ":", str2[j])
	return (comb)
}

## Extract parameters from GLMM results of mortality
extractParams = function(fixef, ppa_status)
{
	climVars = c("min_temperature_of_coldest_month", "I(min_temperature_of_coldest_month^2)",
		"precipitation_of_driest_quarter", "I(precipitation_of_driest_quarter^2)")
	sizeVars = c("dbh", "I(dbh^2)")
	regSlopesName = c("(Intercept)",
		sizeVars,
		climVars)

	# Ref values
	regSlopes = fixef[regSlopesName]

	# Canopy status interactions
	clim_ppa_interact = combineNames("canopy_statusTRUE", climVars)

	# Canopy status modifier
	if (ppa_status)
	{
		regSlopes["(Intercept)"] = regSlopes["(Intercept)"] + fixef["canopy_statusTRUE"]
		regSlopes[climVars] = regSlopes[climVars] + fixef[clim_ppa_interact]
	}
	return (regSlopes)
}

## Mortality function /!\ The dbh is scaled within the function /!\
# Regardless of the estimation algorithm, point estimates are medians computed from simulations.
mortality_fct = function(x, temp, precip, params, scalingMortality, deltaYear = 1)
{
	# Intercept
	beta_0 = params["(Intercept)"] + log(deltaYear)

	# Phi
	beta_1 = params["dbh"]
	beta_2 = params["I(dbh^2)"]

	# Temp
	beta_3 = params["min_temperature_of_coldest_month"]
	beta_4 = params["I(min_temperature_of_coldest_month^2)"]

	# Precip
	beta_5 = params["precipitation_of_driest_quarter"]
	beta_6 = params["I(precipitation_of_driest_quarter^2)"]

	# Z-transform
	x = (x - scalingMortality[var == "dbh", mu])/scalingMortality[var == "dbh", sd]

	# complementary logarithm (p)
	cloglog = unname(beta_0 + beta_1*x + beta_2*x^2 + beta_3*temp + beta_4*temp^2 + beta_5*precip + beta_6*precip^2)

	# 1 - exp(-exp(x)) is the reciprocal function of cloglog
	return ( 1 - exp(-exp(cloglog)) )
}

#### Plot
## Comon variables
nbSpecies = 14

for (i in 1:nbSpecies)
{
	# Get species' names
	loadPath = paste0("array_", i, "/")
	species = list.files(path = loadPath, pattern = "^[0-9]{4,}.*.txt$")
	species = stri_sub(species, to = stri_length(species) - 4)

	# Get parameters (no competition, i.e., canopy status = true)
	fixef = readRDS(paste0(loadPath, "fixef.rds"))
	params_above = extractParams(fixef, TRUE)
	params_below = extractParams(fixef, FALSE)
	scalingMortality = readRDS(paste0(loadPath, "normalisation_mortality_data.rds"))

	# Get dbh limits and dbh min/max from sp-specific dataset
	dbh_lim_m = 0
	dbh_lim_M = 1500

	mortality_dt = readRDS("../createData/mortality_dt.rds")[(species_id == species) & (4 < deltaYear) & (deltaYear < 12), .(dbh)]
	dbh_m = mortality_dt[, min(dbh)]
	dbh_M = mortality_dt[, max(dbh)]

	tikz(paste0(loadPath, species, "_M_dbh.tex"), width = 3.1, height = 3.1) #, standAlone = TRUE)
	op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
	curve(mortality_fct(x, 0, 0, params_below, scalingMortality), from = dbh_lim_m, to = dbh_lim_M,
		xlab = "$ dbh $, in mm", ylab = "$ \\mu(dbh) $, per year", lwd = 2, col = "#3333ff")
	curve(mortality_fct(x, 0, 0, params_above, scalingMortality), from = dbh_lim_m, to = dbh_lim_M,
		lwd = 2, col = "#ff9933", add = TRUE)
	DescTools::Shade(mortality_fct(x, 0, 0, params_below, scalingMortality),
		breaks = c(dbh_m, dbh_M), col = "#A1A1A166", density = NA)
	dev.off()
}
