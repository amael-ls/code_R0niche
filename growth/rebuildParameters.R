
#### Aims of prog:
## Rebuild parameters from LMM results
#		- Read (G)LMM results
#		- Recombine for each species using stringi in params names
#		- Verification using predict and analytic function
#
## Write latex table in a txt file
#		- Writing function
#		- Copy-paste the result in latex appendix app_glmm

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)
library(xtable)
library(lme4)

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

## Extract canopy status only
extract_cs = function(array_id)
{
	(loadPath = paste0("./array_", array_id, "/"))
	model = readRDS(paste0(loadPath, "final_model.rds"))
	fixef_growth = readRDS(paste0(loadPath, "fixef_growth.rds"))

	return (unname(fixef_growth["canopy_statusTRUE"]))
}

## Growth function
growth_fct = function(x, temp, precip, params, species, ppa_status, scalingGrowth)
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

	return (exp(scaling_G_sd[1] * (beta_0 +
		(beta_1 + beta_7*temp + beta_8*temp^2 + beta_9*precip + beta_10*precip^2)*x +
		(beta_2 + beta_11*temp + beta_12*temp^2 + beta_13*precip + beta_14*precip^2)*x^2 +
		beta_3*temp + beta_4*temp^2 +
		beta_5*precip + beta_6*precip^2) + scaling_G_mu))
}

## Growth over- and under- storey, for average clim, and +/- x sd
G_above_below = function(array_id, xsd = 1)
{
	(loadPath = paste0("./array_", array_id, "/"))

	## Load model
	model = readRDS(paste0(loadPath, "final_model.rds"))
	fixef_growth = readRDS(paste0(loadPath, "fixef_growth.rds"))
	scaling = readRDS(paste0(loadPath, "normalisation_growth_data_log.rds"))
	species = stri_sub(list.files(path = loadPath, pattern = ".txt"),
		to = stri_locate_last(list.files(path = loadPath, pattern = ".txt"), regex = ".txt")[1] - 1)

	params_above = extractParams(fixef_growth, species, TRUE)
	params_below = extractParams(fixef_growth, species, FALSE)

	## Climate conditions
	# Average dbh, precipitation and temperature
	dbh0 = 0
	temp0 = 0
	precip0 = 0

	# Average dbh, but precipitation and temperature at average + xsd standard dev (h = high)
	hprecip = xsd
	htemp = xsd

	# Average dbh, but precipitation and temperature at average - xsd standard dev (l = low)
	lprecip = -xsd
	ltemp = -xsd

	## Reference growth (everything at average)
	G_ref_above = growth_fct(dbh0, temp0, precip0, params_above, species, TRUE, scaling)
	G_ref_below = growth_fct(dbh0, temp0, precip0, params_below, species, FALSE, scaling)

	## Precipitation and temperature at +/- 1sd
	# High-high
	G_hh_above = growth_fct(dbh0, htemp, hprecip, params_above, species, TRUE, scaling)
	G_hh_below = growth_fct(dbh0, htemp, hprecip, params_below, species, FALSE, scaling)

	# High-low
	G_hl_above = growth_fct(dbh0, ltemp, hprecip, params_above, species, TRUE, scaling)
	G_hl_below = growth_fct(dbh0, ltemp, hprecip, params_below, species, FALSE, scaling)

	# Low-high
	G_lh_above = growth_fct(dbh0, htemp, lprecip, params_above, species, TRUE, scaling)
	G_lh_below = growth_fct(dbh0, htemp, lprecip, params_below, species, FALSE, scaling)

	# Low-low
	G_ll_above = growth_fct(dbh0, ltemp, lprecip, params_above, species, TRUE, scaling)
	G_ll_below = growth_fct(dbh0, ltemp, lprecip, params_below, species, FALSE, scaling)

	return (list(species,
		G_ref_above, G_ref_below,
		G_hh_above, G_hh_below,
		G_hl_above, G_hl_below,
		G_lh_above, G_lh_below,
		G_ll_above, G_ll_below))
}

## Sorting species (to write their parameters with latex)
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

#### Relative gradient of change
## Get parameters
growth_dt = data.table(array_id = 1:14)
colsName = c("species_id", "ref_a", "ref_b", "hh_a", "hh_b", "hl_a", "hl_b", "lh_a", "lh_b", "ll_a", "ll_b")
growthsName = colsName[colsName != "species_id"]
growth_dt[, c(colsName) := G_above_below(array_id, 1.96), by = array_id]

## At equal temperature (either low or high), but different precipitation
# precip gradient = from low to high, temp = high, i.e., lh -> hh
growth_dt[, lh_hh_x := hh_b - lh_b]
growth_dt[, lh_hh_y := hh_a - lh_a]

# precip gradient = from low to high, temp = low, i.e., ll -> hl
growth_dt[, ll_hl_x := hl_b - ll_b]
growth_dt[, ll_hl_y := hl_a - ll_a]

## At equal precipitation (either low or high), but different temperature
# temp gradient = from low to high, precip = high, i.e., hl -> hh
growth_dt[, hl_hh_x := hh_b - hl_b]
growth_dt[, hl_hh_y := hh_a - hl_a]

# temp gradient = from low to high, precip = low, i.e., ll -> lh
growth_dt[, ll_lh_x := lh_b - ll_b]
growth_dt[, ll_lh_y := lh_a - ll_a]

## Averaged response
# The x- and y- directions are for below and above canopy respectively
growth_dt[, lapply(.SD, function(x) {return (mean(x))}),
	.SDcols = c("lh_hh_x", "lh_hh_y", "ll_hl_x", "ll_hl_y", "hl_hh_x", "hl_hh_y", "ll_lh_x", "ll_lh_y")]

#### Graphics growth overstorey vs understorey
## Common variables
myOrange = rgb(1, 153/255, 0)
myOrange30 = rgb(1, 153/255, 0, 0.3)

myBlue = rgb(51/255, 51/255, 204/255)
myBlue30 = rgb(51/255, 51/255, 204/255, 0.3)

## High/low precip, high temperature
# Blue = High precip, orange = low
axe_min = min(growth_dt[, lapply(.SD, function(x) {return (min(x))}),
	.SDcols = c("hh_b", "hh_a", "lh_b", "lh_a")]) - 0.1
axe_max = max(growth_dt[, lapply(.SD, function(x) {return (max(x))}),
	.SDcols = c("hh_b", "hh_a", "lh_b", "lh_a")]) + 0.1

pdf("./G_over-under-storey_hh-lh.pdf", height = 8, width = 8)
plot(growth_dt[, hh_b], growth_dt[, hh_a], pch = 19, col = myBlue,
	cex = 1.5, cex.axis = 1.5, cex.lab = 1.5,
	xlab = "Understorey growth (mm/yr)", ylab = "overstorey growth (mm/yr)",
	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

points(growth_dt[, lh_b], growth_dt[, lh_a],
	pch = 19, cex = 1.5, col = myOrange)

segments(x0 = growth_dt[, lh_b], y0 = growth_dt[, lh_a],
	x1 = growth_dt[, hh_b], y1 = growth_dt[, hh_a])

title("Warm")

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## High/low precip, low temperature
# Blue = High precip, orange = low
axe_min = min(growth_dt[, lapply(.SD, function(x) {return (min(x))}),
	.SDcols = c("hl_b", "hl_a", "ll_b", "ll_a")]) - 0.1
axe_max = max(growth_dt[, lapply(.SD, function(x) {return (max(x))}),
	.SDcols = c("hl_b", "hl_a", "ll_b", "ll_a")]) + 0.1

pdf("./G_over-under-storey_hl-ll.pdf", height = 8, width = 8)
plot(growth_dt[, ll_b], growth_dt[, ll_a], pch = 19, col = myOrange,
	cex = 1.5, cex.axis = 1.5, cex.lab = 1.5,
	xlab = "Understorey growth (mm/yr)", ylab = "overstorey growth (mm/yr)",
	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

points(growth_dt[, hl_b], growth_dt[, hl_a],
	pch = 19, cex = 1.5, col = myBlue)

segments(x0 = growth_dt[, ll_b], y0 = growth_dt[, ll_a],
	x1 = growth_dt[, hl_b], y1 = growth_dt[, hl_a])

title("Cold")

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## Overcrowded, but all the combination of climate
# Square = low precip, triangle = high precip, blue = cold, orange = warm
axe_min = min(growth_dt[, lapply(.SD, function(x) {return (min(x))}), .SDcols = growthsName]) - 0.1
axe_max = max(growth_dt[, lapply(.SD, function(x) {return (max(x))}), .SDcols = growthsName]) + 0.1

pdf("./G_over-under-storey_allComb.pdf", height = 8, width = 8)
plot(growth_dt[, ref_b], growth_dt[, ref_a], pch = 19, cex = 1.5, cex.axis = 1.5,
	xlab = "Understorey growth (mm/yr)", ylab = "overstorey growth (mm/yr)",
	cex.lab = 1.5, xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

points(growth_dt[, hh_b], growth_dt[, hh_a],
	pch = 17, col = myOrange)
points(growth_dt[, hl_b], growth_dt[, hl_a],
	pch = 17, col = myBlue)
points(growth_dt[, lh_b], growth_dt[, lh_a],
	pch = 15, col = myOrange)
points(growth_dt[, ll_b], growth_dt[, ll_a],
	pch = 15, col = myBlue)

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## Averaged conditions
axe_min = min(growth_dt[, lapply(.SD, function(x) {return (min(x))}),
	.SDcols = c("ref_a", "ref_b")]) - 0.1
axe_max = max(growth_dt[, lapply(.SD, function(x) {return (max(x))}),
	.SDcols = c("ref_a", "ref_b")]) + 0.1

pdf("./G_over-under-storey_averaged.pdf", height = 8, width = 8)
plot(growth_dt[, ref_b], growth_dt[, ref_a], pch = 19, cex = 1.5, cex.axis = 1.5,
	xlab = "Understorey growth (mm/yr)", ylab = "overstorey growth (mm/yr)",
	cex.lab = 1.5, xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## Tikz version of average conditions
# For article
tikz('./G_over-under-storey_averaged.tex', width = 3.1, height = 3.1) #, standAlone = TRUE)
op <- par(mar = c(2.5, 2.5, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(growth_dt[, ref_b], growth_dt[, ref_a], pch = 19, cex = 1.5,
	xlab = "Understorey growth (mm/yr)", ylab = "Overstorey growth (mm/yr)",
	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

# For beamer
# tikz('./G_over-under-storey_averaged.tex', width = 3.1, height = 3.1, standAlone = TRUE)
# op <- par(mar = c(3, 3, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
# plot(growth_dt[, ref_b], growth_dt[, ref_a], pch = 19, cex = 1.5,
# 	cex.lab = 1.5, cex.axis = 1.5,
# 	xlab = "Understorey growth (mm/yr)", ylab = "Overstorey growth (mm/yr)",
# 	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))
#
# # Add identity line
# abline(a = 0, b = 1, lwd = 2)
# dev.off()

#### GLM shade tolerance vs predicted beta canopy status
## Data table tsn
tsn_dt = sortingSpecies(growth_dt[, species_id])
tsn_dt[, c("latin", "vernacular") := .(
	c("Abies balsamea",
		"Acer rubrum", "Acer saccharum",
		"Betula alleghaniensis", "Betula papyrifera",
		"Fagus grandifolia",
		"Picea glauca", "picea mariana", "picea rubens",
		"Pinus banksiana", "Pinus strobus",
		"Populus tremuloides",
		"Thuja occidentalis",
		"Tsuga canadensis"),
	c("Balsam fir",
		"Red maple", "Sugar maple",
		"Yellow birch", "White birch",
		"American beech",
		"White spruce", "Black spruce", "Red spruce",
		"Jack pine", "Eastern white pine",
		"Quaking aspen",
		"Northern white cedar",
		"Eastern hemlock"))]

## Info from Silvics of North America: Hardwoods and Conifers, Burns & Honkala 1990
# Shade tolerance = TRUE, shade intolerant = FALSE. For tol, it can be High, Medium or Low
# Note that the shade tolerance I took is for adult trees, I have very little saplings
tsn_dt[, c("shadeTol", "pageInfo", "tolLevel") := .(
	c(TRUE,
		TRUE, TRUE,
		TRUE, FALSE,
		TRUE,
		TRUE, TRUE, TRUE,
		FALSE,
		TRUE,
		FALSE,
		TRUE, TRUE),
	c("37",
		"171-172", "204",
		"302", "349",
		"661",
		"415", "454", "504",
		"569",
		"987-988",
		"1095",
		"1199", "1247"),
	c("H",
		"M", "H",
		"M", "L",
		"H",
		"M", "M", "H",
		"L",
		"M",
		"L",
		"H", "H"))]

## Collect the slopes coefficients of canopy status (cs)
tsn_dt[, species_id := paste(tsn, species, sep = "-")]
tsn_dt = tsn_dt[growth_dt[, .(array_id, species_id)], on = "species_id"]
setorderv(tsn_dt, "species")
tsn_dt[, cs_coeff := extract_cs(array_id), by = array_id]

## Data plots
# Convert tolerance level tolLevel to a factor
tsn_dt[, tolLevel := factor(tolLevel, levels = c("L", "M", "H"))]

# Plot for the three groups
pdf("./groups.pdf", height = 8, width = 8)
plot(tsn_dt$tolLevel, tsn_dt$cs_coeff,
	xlab = "Shade tolerance level", ylab = "Response of G to light",
	cex.lab = 1.5, cex.axis = 1.5)
dev.off()

## Tikz version
# For article
tikz('groups.tex', width = 3.1, height = 3.1)
op <- par(mar = c(2.5, 2.5, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(tsn_dt$tolLevel, tsn_dt$cs_coeff,
	xlab = "Shade tolerance level", ylab = "Response of $ G $ to light")
dev.off()

# For beamer
# tikz('./groups.tex', width = 3.1, height = 3.1, standAlone = TRUE)
# op <- par(mar = c(3, 3, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
# plot(tsn_dt$tolLevel, tsn_dt$cs_coeff, cex = 1.5,
# 	cex.lab = 1.5, cex.axis = 1.5,
# 	xlab = "Shade tolerance level", ylab = "Response of $ G $ to light")
# dev.off()


## Output tsn_dt
saveRDS(tsn_dt, "./tsn.rds")

tsn_xt = xtable(setcolorder(tsn_dt[, -c("species_id", "shadeTol", "tsn", "pageInfo", "array_id")],
	c("species", "latin", "vernacular", "cs_coeff", "tolLevel")))

print(tsn_xt, file = "./tsn.tex", include.rownames = FALSE, booktabs = TRUE)

#### Extract the species-specific parameters and save them in a data table (to create matlab data)
nbSpecies = 14
parameters_above = data.table(array_id = 1:nbSpecies, species_id = character(nbSpecies))
parameters_below = data.table(array_id = 1:nbSpecies, species_id = character(nbSpecies))
i = 1

for (i in 1:nbSpecies)
{
	(loadPath = paste0("./array_", i, "/"))

	fixef_growth = readRDS(paste0(loadPath, "fixef_growth.rds"))
	species = stri_sub(list.files(path = loadPath, pattern = ".txt"),
		to = stri_locate_last(list.files(path = loadPath, pattern = ".txt"), regex = ".txt")[1] - 1)

	params_above = extractParams(fixef_growth, species, TRUE)
	params_below = extractParams(fixef_growth, species, FALSE)

	parameters_above[i, names(params_above) := as.list(params_above)]
	parameters_above[i, species_id := species]
	parameters_below[i, names(params_below) := as.list(params_below)]
	parameters_below[i, species_id := species]
}

## Change the name for intercept, because of incompatibilit with data table
setnames(parameters_above, "(Intercept)", "Intercept")
setnames(parameters_below, "(Intercept)", "Intercept")

setnames(parameters_above, "I(dbh^2)", "dbh2")
setnames(parameters_below, "I(dbh^2)", "dbh2")

## Save the results
saveRDS(parameters_above, "../createMatlabData/parameters_above_growth.rds")
saveRDS(parameters_below, "../createMatlabData/parameters_below_growth.rds")

# #################
# #### Verification
# array_id = 3
# (loadPath = paste0("./array_", array_id, "/"))
#
# #### Load LMM result
# model = readRDS(paste0(loadPath, "final_model.rds"))
# fixef_growth = readRDS(paste0(loadPath, "fixef_growth.rds"))
# species = stri_sub(list.files(path = loadPath, pattern = ".txt"),
# 	to = stri_locate_last(list.files(path = loadPath, pattern = ".txt"), regex = ".txt")[1] - 1)
#
# print(paste0("species: ", species))
#
# ## Test
# # Using predict lme4
# scaling = readRDS(paste0(loadPath, "normalisation_growth_data_log.rds"))
# newData_above = data.table(species_id = species, dbh = 0.5, annual_mean_temperature = 0.9,
# 	annual_precipitation = -0.2, canopy_status = TRUE)
# newData_below = data.table(species_id = species, dbh = 0.5, annual_mean_temperature = 0.9,
# 	annual_precipitation = -0.2, canopy_status = FALSE)
#
# exp(scaling[var == "growth", sd]*lme4:::predict.merMod(model, newData_above, re.form = NA)
# 	+ scaling[var == "growth", mu])
# exp(scaling[var == "growth", sd]*lme4:::predict.merMod(model, newData_below, re.form = NA)
# 	+ scaling[var == "growth", mu])
#
# # Using growth function
# params_above = extractParams(fixef_growth, species, TRUE)
# params_below = extractParams(fixef_growth, species, FALSE)
# growth_fct(0.5, 0.9, -0.2, params_above, species, TRUE, scaling)
# growth_fct(0.5, 0.9, -0.2, params_below, species, FALSE, scaling)
#
# # To be sure, calculate difference (for whatever reason, it is a vector named "1")
# exp(scaling[var == "growth", sd]*lme4:::predict.merMod(model, newData_above, re.form = NA) +
# 	scaling[var == "growth", mu]) -
# 	growth_fct(0.5, 0.9, -0.2, params_above, species, TRUE, scaling)
#
# exp(scaling[var == "growth", sd]*lme4:::predict.merMod(model, newData_below, re.form = NA) +
# 	scaling[var == "growth", mu]) -
# 	growth_fct(0.5, 0.9, -0.2, params_below, species, FALSE, scaling)
