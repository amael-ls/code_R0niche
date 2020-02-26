
#### Aim of prog: evaluate uncertainty on gowth:
## Explanations:
# For each species, we sample the slopes estimated in the final glmm from 'posterior distributions',
# which are defined here by the confidence intervals. Given we are no in a Bayesian framework, we
# used quoting marks for posterior distributions. Then, we calculate for each set of parameters the
# average of growth. We keep the sets that generate in average the fastest and slowest growth.
# We set the climate to each species-specific average (i.e., 0).
# We also coded the uncertainty for the centroid of the distribution defined by Little 1971, but because
# the centroid might be far from the parametrisation spatial extent, we do not use it.
#
## Load data
#	- Climate of 2010 (average over 5 years: 2006-2010)
#
## Tool functions
#	- To extract parameters from GLMMs
#
## Growth function
#	- Requires climate to be scaled, but not the dbh
#	- Requires the scaling parameters of dbh and growth
#
## Set of parameters
#	- Load confidence intervals
#	- Sample n times for each parameter
#
## Compute average
#	- For each set, integrate in understorey and overstorey conditions
#	- Keep the min and max
#
## Plot
#	- tikzDevice for latex

#### Load package and clear memory
library(data.table)
library(tikzDevice)
library(merTools)
library(stringi)
library(raster)
library(sf)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data and common variables
## Climate year 2010 (average over 5 years, 2006-2010)
# climate2010_stack = stack("../clim60sec/clim_2010.grd", bands = c(1, 12))
# names(climate2010_stack) = c("annual_mean_temperature", "annual_precipitation")

## Common variables
n = 4000
ls_species = readRDS("../createData/ls_species.rds")
nbSpecies = length(ls_species)
nbParametersPerSpecies = 20
species_info = readRDS("../createData/speciesSpecificInfo.rds")
confInt_output = data.table(species = character(nbSpecies * nbParametersPerSpecies),
	term = character(nbSpecies * nbParametersPerSpecies),
	mean = numeric(nbSpecies * nbParametersPerSpecies),
	median = numeric(nbSpecies * nbParametersPerSpecies),
	sd = numeric(nbSpecies * nbParametersPerSpecies))

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

genParams = function(params, confidence_interval, n) # Params might be useless! Check if mean confInt = params!!!
{
	params_dt = data.table()
	varNames = names(params)
	setDT(confidence_interval)
	for (i in 1:length(varNames))
		params_dt[, c(varNames[i]) := rnorm(n, params[varNames[i]], confidence_interval[term == varNames[i], sd])]

	return (params_dt)
}

#### Growth function /!\ The dbh is scaled within the function /!\
growth_fct = function(x, temp, precip, params, species, scalingGrowth)
{
	scaling_G_sd = scalingGrowth[var == "growth", sd]
	scaling_G_mu = scalingGrowth[var == "growth", mu]

	if (!is.data.table(params))
		print("Error from growth: params type transmitted is not a data table")

	# Intercept, note that ppa is already included, according to extractParams fct
	beta_0 = params[["(Intercept)"]]

	# Phi
	beta_1 = params[["dbh"]]
	beta_2 = params[["I(dbh^2)"]]

	# Temp
	beta_3 = params[["annual_mean_temperature"]]
	beta_4 = params[["I(annual_mean_temperature^2)"]]

	# Precip
	beta_5 = params[["annual_precipitation"]]
	beta_6 = params[["I(annual_precipitation^2)"]]

	# Interactions clim:dbh
	beta_7 = params[["dbh:annual_mean_temperature"]]
	beta_8 = params[["dbh:I(annual_mean_temperature^2)"]]
	beta_9 = params[["dbh:annual_precipitation"]]
	beta_10 = params[["dbh:I(annual_precipitation^2)"]]
	beta_11 = params[["I(dbh^2):annual_mean_temperature"]]
	beta_12 = params[["I(dbh^2):I(annual_mean_temperature^2)"]]
	beta_13 = params[["I(dbh^2):annual_precipitation"]]
	beta_14 = params[["I(dbh^2):I(annual_precipitation^2)"]]

	# Z-transform
	x = (x - scalingGrowth[var == "dbh", mu])/scalingGrowth[var == "dbh", sd]

	return (exp(scaling_G_sd[1] * (beta_0 +
		(beta_1 + beta_7*temp + beta_8*temp^2 + beta_9*precip + beta_10*precip^2)*x +
		(beta_2 + beta_11*temp + beta_12*temp^2 + beta_13*precip + beta_14*precip^2)*x^2 +
		beta_3*temp + beta_4*temp^2 +
		beta_5*precip + beta_6*precip^2) + scaling_G_mu))
}

## Read the climatic scalings
# Temperature
scaling_T = fread("../createMatlabData/growthTempScaling.csv")

# Precipitations
scaling_P = fread("../createMatlabData/growthPrecipScaling.csv")

#### Plot
for (i in 1:nbSpecies)
{
	set.seed(1969-08-18) # Woodstock seed

	# Get species' names
	loadPath = paste0("array_", i, "/")
	sp = list.files(path = loadPath, pattern = "^[0-9]{4}.*.txt$")
	sp = stri_sub(sp, to = stri_length(sp) - 4)

	# Get parameters
	fixef = readRDS(paste0(loadPath, "fixef_growth.rds"))
	params_above = extractParams(fixef, sp, TRUE)
	params_below = extractParams(fixef, sp, FALSE)
	scalingGrowth = readRDS(paste0(loadPath, "normalisation_growth_data_log_final.rds"))

	# Get dbh limits and dbh min/max from sp-specific dataset
	dbh_lim_m = 0
	dbh_lim_M = 1500

	growth_dt = readRDS("../createData/growth_dt.rds")[species_id == sp, .(dbh)]
	dbh_m = growth_dt[, min(dbh)]
	dbh_M = growth_dt[, max(dbh)]

	# Get confidence intervals
	conf_int = readRDS(paste0(loadPath, "fe_final_model.rds"))

	# # Get climate centroid (or closest point to centroid), uncomment L 192-208 if centroid
	# pathCentroid = paste0("../little1971/", sp, "/centroid")
	# centroid = st_read(dsn = pathCentroid)
	#
	# centroid = st_transform(centroid,
	# 	crs = st_crs(climate2010_stack))
	#
	# dist = raster::distanceFromPoints(climate2010_stack, centroid)
	# index = raster::which.min(dist)
	#
	# col_id = ifelse(index %% ncol(climate2010_stack) != 0, index %% ncol(climate2010_stack), ncol(climate2010_stack))
	# row_id = (index - col_id)/ncol(climate2010_stack) + 1
	#
	# climate_closestPoint = climate2010_stack[index]
	# temperature = (climate_closestPoint[, "annual_mean_temperature"] - scalingGrowth[var == "annual_mean_temperature", mu])/scalingGrowth[var == "annual_mean_temperature", sd]
	# precipitation = (climate_closestPoint[, "annual_precipitation"] - scalingGrowth[var == "annual_precipitation", mu])/scalingGrowth[var == "annual_precipitation", sd]
	# coords_closestPoint = coordinates(climate2010_stack)[row_id,]

	temperature = 0
	precipitation = 0

	# Generate sets
	parameters_above = genParams(params_above, conf_int, n)
	parameters_below = genParams(params_below, conf_int, n)

	# Keep conf_int for the .csv appendix
	confInt_output[((i - 1)*nbParametersPerSpecies + 1):(i*nbParametersPerSpecies),
		c("term", "mean", "median", "sd"):= conf_int]
	confInt_output[((i - 1)*nbParametersPerSpecies + 1):(i*nbParametersPerSpecies),
		species := sp]

	# Compute average
	result_integrals = data.table(above = numeric(n), below = numeric(n))
	for (j in 1:n)
	{
		if (j %% 100 == 0)
			print(j)
		result_integrals[j, above := integrate(growth_fct, lower = dbh_lim_m, upper = dbh_lim_M,
			temp = temperature, precip = precipitation,
			params = parameters_above[j, ], species = sp, scalingGrowth = scalingGrowth)$value]
		result_integrals[j, below := integrate(growth_fct, lower = dbh_lim_m, upper = dbh_lim_M,
			temp = temperature, precip = precipitation,
			params = parameters_below[j, ], species = sp, scalingGrowth = scalingGrowth)$value]
	}

	# Select the set that generates the highest and lowest for above and below
	index_max_above = which.max(result_integrals[, above])
	index_min_above = which.min(result_integrals[, above])
	index_max_below = which.max(result_integrals[, below])
	index_min_below = which.min(result_integrals[, below])

	# Set limit max y-axis of plot
	# --- Above max
	lim_M_a = optimize(growth_fct, interval = c(dbh_lim_m, dbh_lim_M), maximum = TRUE,
		temp = temperature, precip = precipitation,
		params = parameters_above[index_max_above, ], species = sp, scalingGrowth = scalingGrowth)$objective
	# --- Above min
	lim_m_a = optimize(growth_fct, interval = c(dbh_lim_m, dbh_lim_M), maximum = TRUE,
		temp = temperature, precip = precipitation,
		params = parameters_above[index_min_above, ], species = sp, scalingGrowth = scalingGrowth)$objective
	# --- Below max
	lim_M_b = optimize(growth_fct, interval = c(dbh_lim_m, dbh_lim_M), maximum = TRUE,
		temp = temperature, precip = precipitation,
		params = parameters_above[index_max_below, ], species = sp, scalingGrowth = scalingGrowth)$objective
	# --- Below min
	lim_m_a = optimize(growth_fct, interval = c(dbh_lim_m, dbh_lim_M), maximum = TRUE,
		temp = temperature, precip = precipitation,
		params = parameters_above[index_min_below, ], species = sp, scalingGrowth = scalingGrowth)$objective

	y_max = max(lim_M_a, lim_m_a, lim_M_b, lim_m_a)
	# Plot
	tikz(paste0(loadPath, sp, "_G_dbh.tex"), width = 3.1, height = 3.1,
	packages = c(getOption("tikzLatexPackages"), "\\usetikzlibrary{arrows}"))
	op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
	plot(1, type = "n", xlab = "$ dbh $, in mm", ylab = "$ G(dbh) $, radial growth in mm/yr",
		xlim = c(dbh_lim_m, dbh_lim_M), ylim = c(0, y_max))
	# Above
	curve(growth_fct(x, temperature, precipitation, parameters_above[index_max_above], sp, scalingGrowth),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#ff9933", add = TRUE, lty = "dashed")
	curve(growth_fct(x, temperature, precipitation, parameters_above[index_min_above], sp, scalingGrowth),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#ff9933", add = TRUE, lty = "dashed")

	# Functions we used
	curve(growth_fct(x, temperature, precipitation, as.data.table(as.list(params_above)), sp, scalingGrowth),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#ff9933", add = TRUE)
	curve(growth_fct(x, temperature, precipitation, as.data.table(as.list(params_below)), sp, scalingGrowth),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#3333ff", add = TRUE)

	# Below
	curve(growth_fct(x, temperature, precipitation, parameters_below[index_max_below], sp, scalingGrowth),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#3333ff", add = TRUE, lty = "dotted")
	curve(growth_fct(x, temperature, precipitation, parameters_below[index_min_below], sp, scalingGrowth),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#3333ff", add = TRUE, lty = "dotted")
	# Polygon of parametrisation
	tikzCoord(dbh_m, 0, 'arrow_start')
	tikzCoord(dbh_M, 0, 'arrow_end')
	tikzAnnotate("\\draw[<->,>=stealth,thick,rounded corners=4pt,line width=1pt] (arrow_start) -- (arrow_end);")
	# Polygon growth overstorey
	x_points = seq(dbh_lim_m, dbh_lim_M, length.out = 3000)
	polygon(c(x_points, rev(x_points)),
		c(growth_fct(x_points, temperature, precipitation, parameters_above[index_max_above], sp, scalingGrowth),
		rev(growth_fct(x_points, temperature, precipitation, parameters_above[index_min_above], sp, scalingGrowth))),
		col = "#ff993355", border = NA)
	# Polygone growth understorey
	polygon(c(x_points, rev(x_points)),
		c(growth_fct(x_points, temperature, precipitation, parameters_below[index_max_below], sp, scalingGrowth),
		rev(growth_fct(x_points, temperature, precipitation, parameters_below[index_min_below], sp, scalingGrowth))),
		col = "#3333ff55", border = NA)
	dev.off()

	print(paste0("--- Species: ", sp, " done"))
}

fwrite(confInt_output, "./confInt_allSpecies_growth.csv")
