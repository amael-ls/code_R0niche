
#### Aim of prog: plot the growth in function of dbh, with the range of dbh parametrisation
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

## Extract parameters from GLMM results of growth
extractParams = function(fixef, species, ppa_status)
{
	climVars = c("annual_mean_temperature", "I(annual_mean_temperature^2)",
		"annual_precipitation", "I(annual_precipitation^2)")
	sizeVars = c("dbh", "I(dbh^2)")
	regSlopesName = c("(Intercept)",
		sizeVars,
		climVars,
		combineNames(sizeVars, climVars))

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

## Growth function /!\ The dbh is scaled within the function /!\
growth_fct = function(x, temp, precip, params, species, scalingGrowth)
{
	scaling_G_sd = scalingGrowth$sd[1]
	scaling_G_mu = scalingGrowth$mu[1]

	# Intercept, note that ppa is already included, according to extractParams fct
	beta_0 = params["(Intercept)"]

	# Phi
	beta_1 = params["dbh"]
	beta_2 = params["I(dbh^2)"]

	# Temp
	beta_3 = params["annual_mean_temperature"]
	beta_4 = params["I(annual_mean_temperature^2)"]

	# Precip
	beta_5 = params["annual_precipitation"]
	beta_6 = params["I(annual_precipitation^2)"]

	# Interactions clim:dbh
	beta_7 = params["dbh:annual_mean_temperature"]
	beta_8 = params["dbh:I(annual_mean_temperature^2)"]
	beta_9 = params["dbh:annual_precipitation"]
	beta_10 = params["dbh:I(annual_precipitation^2)"]
	beta_11 = params["I(dbh^2):annual_mean_temperature"]
	beta_12 = params["I(dbh^2):I(annual_mean_temperature^2)"]
	beta_13 = params["I(dbh^2):annual_precipitation"]
	beta_14 = params["I(dbh^2):I(annual_precipitation^2)"]

	# Z-transform
	x = (x - scalingGrowth[var == "dbh", mu])/scalingGrowth[var == "dbh", sd]

	return (exp(scaling_G_sd[1] * (beta_0 +
		(beta_1 + beta_7*temp + beta_8*temp^2 + beta_9*precip + beta_10*precip^2)*x +
		(beta_2 + beta_11*temp + beta_12*temp^2 + beta_13*precip + beta_14*precip^2)*x^2 +
		beta_3*temp + beta_4*temp^2 +
		beta_5*precip + beta_6*precip^2) + scaling_G_mu))
}

#### Plot
## Comon variables
nbSpecies = 14

for (i in 1:nbSpecies)
{
	# Get species' names
	loadPath = paste0("array_", i, "/")
	species = list.files(path = loadPath, pattern = ".txt$")
	species = stri_sub(species, to = stri_length(species) - 4)

	# Get parameters (no competition, i.e., canopy status = true)
	fixef = readRDS(paste0(loadPath, "fixef_growth.rds"))
	params_above = extractParams(fixef, species, TRUE)
	params_below = extractParams(fixef, species, FALSE)
	scalingGrowth = readRDS(paste0(loadPath, "normalisation_growth_data_log_final.rds"))

	# Get dbh limits and dbh min/max from sp-specific dataset
	dbh_lim_m = 0
	dbh_lim_M = 1500

	growth_dt = readRDS("../createData/growth_dt.rds")[species_id == species, .(dbh)]
	dbh_m = growth_dt[, min(dbh)]
	dbh_M = growth_dt[, max(dbh)]

	tikz(paste0(loadPath, species, "_G_dbh.tex"), width = 3.1, height = 3.1) #, standAlone = TRUE)
	op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
	curve(growth_fct(x, 0, 0, params_above, species, scalingGrowth), from = dbh_lim_m, to = dbh_lim_M,
		xlab = "$ dbh $, in mm", ylab = "$ G(dbh) $, in mm/yr", lwd = 2, col = "#ff9933")
	curve(growth_fct(x, 0, 0, params_below, species, scalingGrowth), from = dbh_lim_m, to = dbh_lim_M,
		lwd = 2, col = "#3333ff", add = TRUE)
	DescTools::Shade(growth_fct(x, 0, 0, params_above, species, scalingGrowth),
		breaks = c(dbh_m, dbh_M), col = "#A1A1A166", density = NA)
	dev.off()
}
