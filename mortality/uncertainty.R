
#### Aim of prog: evaluate uncertainty on mortality:
## Explanations:
# For each species, we sample the posterior distributions to generate different sets of parameters.
# Then, we calculate for each set of parameters the average of mortality. We keep the sets that generate
# in average the highest and lowest mortality. We set the climate to each species-specific average (i.e., 0).
# We also coded the uncertainty for the centroid of the distribution defined by Little 1971, but because
# the centroid might be far from the parametrisation spatial extent, we do not use it.
#
## Load data
#	- Climate of 2010 (average over 5 years: 2006-2010)
#
## Mortality function
#	- Requires climate to be scaled, but not the dbh
#	- Requires the scaling parameters of dbh
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

#### Load packages
library(data.table)
library(tikzDevice)
library(rstanarm)
library(stringi)
library(raster)
library(sf)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data and common variables (Uncomment L 43-44 if centroid)
## Climate year 2010 (average over 5 years, 2006-2010)
# climate2010_stack = stack("../clim60sec/clim_2010.grd", bands = c(1, 12))
# names(climate2010_stack) = c("min_temperature_of_coldest_month", "precipitation_of_driest_quarter")

## Common variables
ls_folders = list.files(path = "./", pattern = "array_[0-9]{1,}")
nbFolders = length(ls_folders)
chosenModel = "model_7.rds"
nbParametersPerSpecies = 12
nbParametersMortality_fct = 7
species_info = readRDS("../createData/speciesSpecificInfo.rds")
confInt_output = data.table(species = character(nbFolders * nbParametersPerSpecies),
	term = character(nbFolders * nbParametersPerSpecies),
	min_5p = numeric(nbFolders * nbParametersPerSpecies),
	median = numeric(nbFolders * nbParametersPerSpecies),
	max_95p = numeric(nbFolders * nbParametersPerSpecies))

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
genParams = function(fixef, ppa_status)
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
mortality_fct = function(x, temp, precip, params, scalingMortality)
{
	# Intercept
	beta_0 = params[["(Intercept)"]]

	# Phi
	beta_1 = params[["dbh"]]
	beta_2 = params[["I(dbh^2)"]]

	# Temp
	beta_3 = params[["min_temperature_of_coldest_month"]]
	beta_4 = params[["I(min_temperature_of_coldest_month^2)"]]

	# Precip
	beta_5 = params[["precipitation_of_driest_quarter"]]
	beta_6 = params[["I(precipitation_of_driest_quarter^2)"]]

	# Z-transform
	x = (x - scalingMortality[var == "dbh", mu])/scalingMortality[var == "dbh", sd]

	# logit (p)
	logit_p = unname(beta_0 + beta_1*x + beta_2*x^2 + beta_3*temp + beta_4*temp^2 + beta_5*precip + beta_6*precip^2)

	# 1/(1 + exp(-logit_p)) is the sigmoid, the reciprocal function of logit
	return ( 1/(1 + exp(-logit_p)) )
}

#### Plot
for (i in 1:nbFolders)
{
	# Get species' names
	loadPath = paste0("array_", i, "/")
	sp = list.files(path = loadPath, pattern = "^[0-9]{4,}.*.txt$")
	sp = stri_sub(sp, to = stri_length(sp) - 4)

	# Get parameters (no competition, i.e., canopy status = true)
	fixef = readRDS(paste0(loadPath, "fixef.rds"))
	params_above = genParams(fixef, TRUE)
	params_below = genParams(fixef, FALSE)
	mortalityParamsNames = names(params_below)
	scalingMortality = readRDS(paste0(loadPath, "normalisation_mortality_data.rds"))

	# Get dbh limits and dbh min/max from sp-specific dataset
	dbh_lim_m = 0
	dbh_lim_M = 1500

	mortality_dt = readRDS("../createData/mortality_dt.rds")[species_id == sp, .(dbh)]
	dbh_m = mortality_dt[, min(dbh)]
	dbh_M = mortality_dt[, max(dbh)]

	# Confidence interval
	model = readRDS(paste0(loadPath, chosenModel))
	conf_int = posterior_interval(model, prob = 0.90) # http://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.html

	conf_int = conf_int[stri_detect(str = rownames(conf_int), regex = "^(?!b\\[\\(Intercept\\))"), ]
	conf_int = conf_int[stri_detect(str = rownames(conf_int), regex = "^(?!Sigma)"), ]

	conf_int = as.data.frame(conf_int)
	setDT(conf_int, keep.rownames = TRUE)
	setnames(conf_int, old = "rn", new = "term")

	# Keep conf_int for the .csv appendix
	confInt_output[((i - 1)*nbParametersPerSpecies + 1):(i*nbParametersPerSpecies),
		c("term", "min_5p", "max_95p") := conf_int]
	confInt_output[((i - 1)*nbParametersPerSpecies + 1):(i*nbParametersPerSpecies),
		median := model$coefficients[conf_int[["term"]]]]
	confInt_output[((i - 1)*nbParametersPerSpecies + 1):(i*nbParametersPerSpecies),
		species := sp]

	# # Get climate centroid (or closest point to centroid), uncomment L 170-185 if centroid
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
	# temperature = (climate_closestPoint[, "min_temperature_of_coldest_month"] - scalingMortality[var == "min_temperature_of_coldest_month", mu])/scalingMortality[var == "min_temperature_of_coldest_month", sd]
	# precipitation = (climate_closestPoint[, "precipitation_of_driest_quarter"] - scalingMortality[var == "precipitation_of_driest_quarter", mu])/scalingMortality[var == "precipitation_of_driest_quarter", sd]
	# coords_closestPoint = coordinates(climate2010_stack)[row_id,]

	temperature = 0
	precipitation = 0

	# Generate sets
	parameters = as.matrix(model, pars = conf_int[["term"]])
	parameters_above = setnames(data.table(matrix(data = 0, nrow = nrow(parameters), ncol = nbParametersMortality_fct)),
		mortalityParamsNames)
	parameters_below = setnames(data.table(matrix(data = 0, nrow = nrow(parameters), ncol = nbParametersMortality_fct)),
		mortalityParamsNames)
	for (j in 1:nrow(parameters))
	{
		parameters_above[j, (mortalityParamsNames) := as.list(genParams(parameters[j,], TRUE))]
		parameters_below[j, (mortalityParamsNames) := as.list(genParams(parameters[j,], FALSE))]
	}

	n = nrow(parameters)

	# Compute average
	result_integrals = data.table(above = numeric(n), below = numeric(n))
	for (j in 1:n)
	{
		if (j %% 100 == 0)
			print(j)
		result_integrals[j, above := integrate(mortality_fct, lower = dbh_lim_m, upper = dbh_lim_M,
			temp = temperature, precip = precipitation,
			params = parameters_above[j, ], scalingMortality = scalingMortality)$value]
		result_integrals[j, below := integrate(mortality_fct, lower = dbh_lim_m, upper = dbh_lim_M,
			temp = temperature, precip = precipitation,
			params = parameters_below[j, ], scalingMortality = scalingMortality)$value]
	}

	# Select the set that generates the highest and lowest for above and below
	index_max_above = which.max(result_integrals[, above])
	index_min_above = which.min(result_integrals[, above])
	index_max_below = which.max(result_integrals[, below])
	index_min_below = which.min(result_integrals[, below])

	# Plot
	tikz(paste0(loadPath, sp, "_M_dbh.tex"), width = 3.1, height = 3.1,
	packages = c(getOption("tikzLatexPackages"), "\\usetikzlibrary{arrows}"))
	op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
	plot(1, type = "n", xlab = "$ dbh $, in mm", ylab = "$ \\mu(dbh) $, probability of dying per year",
		xlim = c(dbh_lim_m, dbh_lim_M), ylim = c(0, 1))
	# Above
	curve(mortality_fct(x, temperature, precipitation, parameters_above[index_max_above], scalingMortality),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#ff9933", add = TRUE, lty = "dashed")
	curve(mortality_fct(x, temperature, precipitation, parameters_above[index_min_above], scalingMortality),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#ff9933", add = TRUE, lty = "dashed")

	# Functions we used
	curve(mortality_fct(x, temperature, precipitation, as.data.table(as.list(params_above)), scalingMortality),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#ff9933", add = TRUE)
	curve(mortality_fct(x, temperature, precipitation, as.data.table(as.list(params_below)), scalingMortality),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#3333ff", add = TRUE)

	# Below
	curve(mortality_fct(x, temperature, precipitation, parameters_below[index_max_below], scalingMortality),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#3333ff", add = TRUE, lty = "dotted")
	curve(mortality_fct(x, temperature, precipitation, parameters_below[index_min_below], scalingMortality),
		from = dbh_lim_m, to = dbh_lim_M, lwd = 2, col = "#3333ff", add = TRUE, lty = "dotted")
	# Polygon of parametrisation
	tikzCoord(dbh_m, 0, 'arrow_start')
	tikzCoord(dbh_M, 0, 'arrow_end')
	tikzAnnotate("\\draw[<->,>=stealth,thick,rounded corners=4pt,line width=1pt] (arrow_start) -- (arrow_end);")
	# Polygon mortality overstorey
	x_points = seq(dbh_lim_m, dbh_lim_M, length.out = 3000)
	polygon(c(x_points, rev(x_points)),
		c(mortality_fct(x_points, temperature, precipitation, parameters_above[index_max_above], scalingMortality),
		rev(mortality_fct(x_points, temperature, precipitation, parameters_above[index_min_above], scalingMortality))),
		col = "#ff993355", border = NA)
	# Polygone mortality understorey
	polygon(c(x_points, rev(x_points)),
		c(mortality_fct(x_points, temperature, precipitation, parameters_below[index_max_below], scalingMortality),
		rev(mortality_fct(x_points, temperature, precipitation, parameters_below[index_min_below], scalingMortality))),
		col = "#3333ff55", border = NA)
	dev.off()

	print(paste0("--- Species: ", sp, " done"))
}

fwrite(confInt_output, "./confInt_allSpecies_mortality.csv")
