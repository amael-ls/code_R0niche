
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
library(rstanarm)
library(stringi)
library(xtable)

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

## Extract canopy status only
extract_cs = function(array_id)
{
	fixef = readRDS(paste0("array_", array_id, "/fixef.rds"))
	return (unname(fixef["canopy_statusTRUE"]))
}

## Regardless of the estimation algorithm, point estimates are medians computed from simulations.
## Mortality function, note that ppa is already included, according to extractParams fct
mortality_fct = function(x, temp, precip, params, deltaYear = 1)
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

	# complementary logarithm (p)
	cloglog = unname(beta_0 + beta_1*x + beta_2*x^2 + beta_3*temp + beta_4*temp^2 + beta_5*precip + beta_6*precip^2)

	# 1 - exp(-exp(x)) is the reciprocal function of cloglog
	return ( 1 - exp(-exp(cloglog)) )
}

## Growth over- and under- storey, for average clim, and +/- x sd
M_above_below = function(array_id, xsd = 1)
{
	loadPath = paste0("./array_", array_id, "/")

	## Load fixef
	fixef_mortality = readRDS(paste0(loadPath, "fixef.rds"))
	species = stri_sub(list.files(path = loadPath, pattern = ".txt"),
		to = stri_locate_last(list.files(path = loadPath, pattern = ".txt"), regex = ".txt")[1] - 1)

	species = species[stri_detect(species, regex = "^[0-9]{4,}")]

	params_above = extractParams(fixef_mortality, TRUE)
	params_below = extractParams(fixef_mortality, FALSE)

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

	## Reference mortality (everything at average)
	M_ref_above = mortality_fct(dbh0, temp0, precip0, params_above)
	M_ref_below = mortality_fct(dbh0, temp0, precip0, params_below)

	## Precipitation and temperature at +/- 1sd
	# High-high
	M_hh_above = mortality_fct(dbh0, htemp, hprecip, params_above)
	M_hh_below = mortality_fct(dbh0, htemp, hprecip, params_below)

	# High-low
	M_hl_above = mortality_fct(dbh0, ltemp, hprecip, params_above)
	M_hl_below = mortality_fct(dbh0, ltemp, hprecip, params_below)

	# Low-high
	M_lh_above = mortality_fct(dbh0, htemp, lprecip, params_above)
	M_lh_below = mortality_fct(dbh0, htemp, lprecip, params_below)

	# Low-low
	M_ll_above = mortality_fct(dbh0, ltemp, lprecip, params_above)
	M_ll_below = mortality_fct(dbh0, ltemp, lprecip, params_below)

	return (list(species,
		M_ref_above, M_ref_below,
		M_hh_above, M_hh_below,
		M_hl_above, M_hl_below,
		M_lh_above, M_lh_below,
		M_ll_above, M_ll_below))
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
mortality_dt = data.table(array_id = 1:14)
colsName = c("species_id", "ref_a", "ref_b", "hh_a", "hh_b", "hl_a", "hl_b", "lh_a", "lh_b", "ll_a", "ll_b")
mortalityName = colsName[colsName != "species_id"]
mortality_dt[, c(colsName) := M_above_below(array_id, 1.96), by = array_id]

## At equal temperature (either low or high), but different precipitation
# precip gradient = from low to high, temp = high, i.e., lh -> hh
mortality_dt[, lh_hh_x := hh_b - lh_b]
mortality_dt[, lh_hh_y := hh_a - lh_a]

# precip gradient = from low to high, temp = low, i.e., ll -> hl
mortality_dt[, ll_hl_x := hl_b - ll_b]
mortality_dt[, ll_hl_y := hl_a - ll_a]

## At equal precipitation (either low or high), but different temperature
# temp gradient = from low to high, precip = high, i.e., hl -> hh
mortality_dt[, hl_hh_x := hh_b - hl_b]
mortality_dt[, hl_hh_y := hh_a - hl_a]

# temp gradient = from low to high, precip = low, i.e., ll -> lh
mortality_dt[, ll_lh_x := lh_b - ll_b]
mortality_dt[, ll_lh_y := lh_a - ll_a]

## Averaged response
# The x- and y- directions are for below and above canopy respectively
mortality_dt[, lapply(.SD, function(x) {return (mean(x))}),
	.SDcols = c("lh_hh_x", "lh_hh_y", "ll_hl_x", "ll_hl_y", "hl_hh_x", "hl_hh_y", "ll_lh_x", "ll_lh_y")]

#### Graphics mortality overstorey vs understorey
## Common variables
myOrange = rgb(1, 153/255, 0)
myOrange30 = rgb(1, 153/255, 0, 0.3)

myBlue = rgb(51/255, 51/255, 204/255)
myBlue30 = rgb(51/255, 51/255, 204/255, 0.3)

## High/low precip, high temperature
# Blue = High precip, orange = low
axe_min = min(mortality_dt[, lapply(.SD, function(x) {return (min(x))}),
	.SDcols = c("hh_b", "hh_a", "lh_b", "lh_a")]) - 0.01
axe_max = max(mortality_dt[, lapply(.SD, function(x) {return (max(x))}),
	.SDcols = c("hh_b", "hh_a", "lh_b", "lh_a")]) + 0.01

pdf("./M_over-under-storey_hh-lh.pdf", height = 8, width = 8)
plot(mortality_dt[, hh_b], mortality_dt[, hh_a], pch = 19, col = myBlue,
	cex = 1.5, cex.axis = 1.5, cex.lab = 1.5,
	xlab = "Understorey mortality", ylab = "overstorey mortality",
	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

points(mortality_dt[, lh_b], mortality_dt[, lh_a],
	pch = 19, cex = 1.5, col = myOrange)

segments(x0 = mortality_dt[, lh_b], y0 = mortality_dt[, lh_a],
	x1 = mortality_dt[, hh_b], y1 = mortality_dt[, hh_a])

title("Warm")

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## High/low precip, low temperature
# Blue = High precip, orange = low
axe_min = min(mortality_dt[, lapply(.SD, function(x) {return (min(x))}),
	.SDcols = c("hl_b", "hl_a", "ll_b", "ll_a")]) - 0.01
axe_max = max(mortality_dt[, lapply(.SD, function(x) {return (max(x))}),
	.SDcols = c("hl_b", "hl_a", "ll_b", "ll_a")]) + 0.01

pdf("./M_over-under-storey_hl-ll.pdf", height = 8, width = 8)
plot(mortality_dt[, ll_b], mortality_dt[, ll_a], pch = 19, col = myOrange,
	cex = 1.5, cex.axis = 1.5, cex.lab = 1.5,
	xlab = "Understorey mortality", ylab = "overstorey mortality",
	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

points(mortality_dt[, hl_b], mortality_dt[, hl_a],
	pch = 19, cex = 1.5, col = myBlue)

segments(x0 = mortality_dt[, ll_b], y0 = mortality_dt[, ll_a],
	x1 = mortality_dt[, hl_b], y1 = mortality_dt[, hl_a])

title("Cold")

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## Overcrowded, but all the combination of climate
# Square = low precip, triangle = high precip, blue = cold, orange = warm
axe_min = min(mortality_dt[, lapply(.SD, function(x) {return (min(x))}), .SDcols = mortalityName]) - 0.01
axe_max = max(mortality_dt[, lapply(.SD, function(x) {return (max(x))}), .SDcols = mortalityName]) + 0.01

pdf("./M_over-under-storey_allComb.pdf", height = 8, width = 8)
plot(mortality_dt[, ref_b], mortality_dt[, ref_a], pch = 19, cex = 1.5, cex.axis = 1.5,
	xlab = "Understorey mortality", ylab = "overstorey mortality",
	cex.lab = 1.5, xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

points(mortality_dt[, hh_b], mortality_dt[, hh_a],
	pch = 17, col = myOrange)
points(mortality_dt[, hl_b], mortality_dt[, hl_a],
	pch = 17, col = myBlue)
points(mortality_dt[, lh_b], mortality_dt[, lh_a],
	pch = 15, col = myOrange)
points(mortality_dt[, ll_b], mortality_dt[, ll_a],
	pch = 15, col = myBlue)

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## Averaged conditions
axe_min = ifelse(min(mortality_dt[, lapply(.SD, function(x) {return (min(x))}), .SDcols = c("ref_a", "ref_b")]) - 0.005 < 0,
	0, min(mortality_dt[, lapply(.SD, function(x) {return (min(x))}), .SDcols = c("ref_a", "ref_b")]) - 0.005)
axe_max = max(mortality_dt[, lapply(.SD, function(x) {return (max(x))}),
	.SDcols = c("ref_a", "ref_b")]) + 0.005

pdf("./M_over-under-storey_averaged.pdf", height = 8, width = 8)
plot(mortality_dt[, ref_b], mortality_dt[, ref_a], pch = 19, cex = 1.5, cex.axis = 1.5,
	xlab = "Understorey mortality", ylab = "overstorey mortality",
	cex.lab = 1.5, xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

## Tikz version
# For article
tikz('./M_over-under-storey_averaged.tex', width = 3.1, height = 3.1)
op <- par(mar = c(2.5, 2.5, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(mortality_dt[, ref_b], mortality_dt[, ref_a], pch = 19,
	xlab = "Understorey mortality $ \\mu $ (yr\\textsuperscript{-1})", ylab = "Overstorey mortality $ \\mu $ (yr\\textsuperscript{-1})",
	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))

# Add identity line
abline(a = 0, b = 1, lwd = 2)
dev.off()

# For beamer
# tikz('./M_over-under-storey_averaged.tex', width = 3.1, height = 3.1, standAlone = TRUE)
# op <- par(mar = c(3, 3, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
# plot(mortality_dt[, ref_b], mortality_dt[, ref_a], pch = 19,
# 	cex.lab = 1.5, cex.axis = 1.5,
# 	xlab = "Understorey mortality", ylab = "Overstorey mortality",
# 	xlim = c(axe_min, axe_max), ylim = c(axe_min, axe_max))
#
# # Add identity line
# abline(a = 0, b = 1, lwd = 2)
# dev.off()

#### GLM shade tolerance vs predicted beta canopy status
## Data table tsn
tsn_dt = sortingSpecies(mortality_dt[, species_id])
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
tsn_dt = tsn_dt[mortality_dt[, .(array_id, species_id)], on = "species_id"]
setorderv(tsn_dt, "species")
tsn_dt[, cs_coeff := extract_cs(array_id), by = array_id]

## Data plots
# Convert tolerance level tolLevel to a factor
tsn_dt[, tolLevel := factor(tolLevel, levels = c("L", "M", "H"))]

# Plot for the three groups
pdf("./groups.pdf", height = 8, width = 8)
plot(tsn_dt$tolLevel, tsn_dt$cs_coeff,
	xlab = "Shade tolerance level", ylab = "Response of mu to light",
	cex.lab = 1.5, cex.axis = 1.5)
dev.off()

# Tikz version
# For article
tikz('groups_mortality.tex', width = 3.1, height = 3.1)
op <- par(mar = c(2.5, 2.5, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(tsn_dt$tolLevel, tsn_dt$cs_coeff,
	xlab = "Shade tolerance level", ylab = "Response of mortality ($ \\mu $) to light")
dev.off()

# For beamer
# tikz('groups_mortality.tex', width = 3.1, height = 3.1, standAlone = TRUE)
# op <- par(mar = c(3, 3, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
# plot(tsn_dt$tolLevel, tsn_dt$cs_coeff,
# 	cex.lab = 1.5, cex.axis = 1.5,
# 	xlab = "Shade tolerance level", ylab = "Response of $ \\mu $ to light")
# dev.off()

## Output tsn_dt
tsn_xt = xtable(setcolorder(tsn_dt[, -c("species_id", "shadeTol", "tsn", "pageInfo", "array_id")],
	c("species", "latin", "vernacular", "cs_coeff", "tolLevel")))

print(tsn_xt, file = "./tsn.tex", include.rownames = FALSE, booktabs = TRUE)

#### Extract the species-specific parameters and save them in a data table (to create matlab data)
nbSpecies = 14
parameters_above = data.table(array_id = 1:nbSpecies, species_id = character(nbSpecies))
parameters_below = data.table(array_id = 1:nbSpecies, species_id = character(nbSpecies))

for (i in 1:nbSpecies)
{
	(loadPath = paste0("./array_", i, "/"))

	fixef_mortality = readRDS(paste0(loadPath, "fixef.rds"))
	species = stri_sub(list.files(path = loadPath, pattern = ".txt"),
		to = stri_locate_last(list.files(path = loadPath, pattern = ".txt"), regex = ".txt")[1] - 1)

	species = species[stri_detect(species, regex = "^[0-9]{4,}")]

	params_above = extractParams(fixef_mortality, TRUE)
	params_below = extractParams(fixef_mortality, FALSE)

	parameters_above[i, names(params_above) := as.list(params_above)]
	parameters_above[i, species_id := species]
	parameters_below[i, names(params_below) := as.list(params_below)]
	parameters_below[i, species_id := species]
}

## Change colnames because of incompatibility with data table
setnames(parameters_above, "(Intercept)", "Intercept")
setnames(parameters_below, "(Intercept)", "Intercept")

setnames(parameters_above, "I(dbh^2)", "dbh2")
setnames(parameters_below, "I(dbh^2)", "dbh2")

## Save the results
saveRDS(parameters_above, "../createMatlabData/parameters_above_mortality.rds")
saveRDS(parameters_below, "../createMatlabData/parameters_below_mortality.rds")

# #################
# #### Verification
# array_id = 3
# (loadPath = paste0("./array_", array_id, "/"))
#
# #### Load LMM result
# model = readRDS(paste0(loadPath, "/model_7.rds"))
# fixef = readRDS(paste0(loadPath, "./fixef.rds"))
# species = stri_sub(list.files(path = loadPath, pattern = ".txt"),
# 	to = stri_locate_last(list.files(path = loadPath, pattern = ".txt"), regex = ".txt")[1] - 1)
#
# print(paste0("species: ", species))
#
# ## Test
# # Using posterior_predict rstanarm 200 times and check if mean converge to what I want
# newData_above = data.table(species_id = species, dbh = 0.5, min_temperature_of_coldest_month = 0.9,
# 	precipitation_of_driest_quarter = -0.2, deltaYear = 6, canopy_status = as.factor(TRUE))
# newData_below = data.table(species_id = species, dbh = 0.5, min_temperature_of_coldest_month = 0.9,
# 	precipitation_of_driest_quarter = -0.2, deltaYear = 6, canopy_status = as.factor(FALSE))
#
# aa = numeric(2e2)
# bb = aa
# for (i in 1:2e2)
# {
# 	pred_above = posterior_predict(model, newData_above, offset = log(newData_above[, deltaYear]), re.form = NA)
# 	pred_below = posterior_predict(model, newData_below, offset = log(newData_below[, deltaYear]), re.form = NA)
# 	aa[i] = sum(pred_above)/length(pred_above)
# 	bb[i] = sum(pred_below)/length(pred_below)
# 	if (i %% 10 == 0)
# 		print(i)
# }
#
# # Using mortality function
# params_above = extractParams(fixef, TRUE)
# params_below = extractParams(fixef, FALSE)
# mortality_fct(0.5, 0.9, -0.2, params_above, deltaYear = 6)
# mortality_fct(0.5, 0.9, -0.2, params_below, deltaYear = 6)
#
# mean(aa)
# mean(bb)
#
# (mortality_fct(0.5, 0.9, -0.2, params_above, deltaYear = 6) - mean(aa))/mean(aa)*100
# (mortality_fct(0.5, 0.9, -0.2, params_below, deltaYear = 6) - mean(bb))/mean(bb)*100
