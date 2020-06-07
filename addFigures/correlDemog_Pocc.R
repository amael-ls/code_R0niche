
#### Aim of prog: Compute correlation between vital rates and Proba occurrence
## Growth -- Proba occurrence
#	- Load growth data
#	- Load random forest
#	- Compute species-specific proba occurence
#
## Mortality -- Proba occurrence
#	- Load mortality data
#	- Compute species-specific proba occurence
#
## Remark
# Given the high variability of growth and mortality all along the distribution,
# I do not expect any trend here.

#### Load package and clear memory
library(randomForest)
library(data.table)
library(tikzDevice)
library(stringi)
library(xtable)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool functions
## To get the name of the Species Distribution model
getSDMname = function(str)
{
	# Replace underscores by scores
	species = stri_replace_all(str = str, replacement = "_", regex = "-")

	# Switch the Taxonomic Serial Number (tsn) with species code
	tsn_pos = stri_locate_first(str = species, regex = "_")[, 1]
	tsn = stri_sub(str = species, to = tsn_pos - 1)
	sp_code = stri_sub(str = species, from = tsn_pos + 1)
	species = paste0(sp_code, "_", tsn)

	results = data.table(sp_code = stri_replace_all(str = sp_code, replacement = "-", regex = "_"),
	sdm = species, tsn = tsn)

	return (results)
}

## Function to plot (quick and dirty, unsing global access variable)
correlPlot_fct = function(nameFig, rangeFig = 1:3, stand = FALSE)
{
	spPerFig = length(rangeFig)
	tikz(nameFig, width = 6, height = 5, standAlone = stand)
	op = par(mar = c(0, 2.5, 2, 0), mgp = c(1.5, 0.75, 0),
		oma = c(0,2,0,0), tck = -0.015, las = 1)

	tikzAnnotate(paste0("\\definecolor{shade_under}{RGB}{", paste(shade_under, collapse = ","), "}"))
	tikzAnnotate(paste0("\\definecolor{shade_over}{RGB}{", paste(shade_over, collapse = ","), "}"))

	plot(x = NULL, y = NULL, xlim = c(0, spPerFig*speciesSpace + (spPerFig - 1)*interSpecies + 0.1),
		ylim = c(-1, 1 + 0.2), axes = FALSE, xlab = "",
		ylab = "")

	count = 0
	for (i in rangeFig)
	{
		## Species line and name
		sp = correlation[i, sp_code]
		species_coords = halfSpeciesSpace + count*(speciesSpace + interSpecies)

		# Horizontal line species
		tikzCoord(species_coords - halfSpeciesSpace, 1 + 0.2, paste0("pt_start_", i))
		tikzCoord(species_coords + halfSpeciesSpace, 1 + 0.2, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(species_coords, 1 + 0.2, 'species_coords')

		# Species name
		tikzAnnotate(paste0("\\node[above] (spPos) at (species_coords) {", sp, "};"))

		# Vertical line
		tikzCoord(species_coords, -1, paste0("vert_start_", i))
		tikzCoord(species_coords, 1 + 0.1, paste0("vert_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (vert_start_", i, ") -- (vert_end_", i, ");"))

		## Categories 15 cm (left) and 60 cm (right)
		tikzCoord(species_coords - interDot - interDot/2, 1 + 0.1, '15cm_coords')
		tikzAnnotate(paste0("\\node (15cmPos) at (15cm_coords) {15 cm};"))

		tikzCoord(species_coords + interDot + interDot/2, 1 + 0.1, '60cm_coords')
		tikzAnnotate(paste0("\\node (60cmPos) at (60cm_coords) {60 cm};"))

		## Add points correlation demography <--> Pocc
		# Mortality 15
		points(x = rep(species_coords - 2*interDot, 2),
			y = correlation[i, .(correl_15cm_mortality_above, correl_15cm_mortality_below)],
			col = c("#0E55FF", "#010166"), pch = 4)

		# Growth 15
		points(x = rep(species_coords - interDot, 2),
			y = correlation[i, .(correl_15cm_growth_above, correl_15cm_growth_below)],
			col = c("#0E55FF", "#010166"), pch = 16)

		# Mortality 60
		points(x = rep(species_coords + interDot, 2),
			y = correlation[i, .(correl_60cm_mortality_above, correl_60cm_mortality_below)],
			col = c("#0E55FF", "#010166"), pch = 4)

		# Growth 60
		points(x = rep(species_coords + 2*interDot, 2),
			y = correlation[i, .(correl_60cm_growth_above, correl_60cm_growth_below)],
			col = c("#0E55FF", "#010166"), pch = 16)

		## Add points correlation R0 <--> Pocc
		points(x = rep(species_coords, 2),
			y = correlation[i, .(correl_0m, correl_10m)],
			col = c("#0E55FF", "#010166"), pch = 17)

		count = count + 1
	}

	axis(side = 2, at = seq(-1, 1, by = 0.25),
		labels = c("-1", "", "-0.5", "", "0", "", "0.5", "", "1"))

	# segments(x0 = 0, y0 = 0, x1 = spPerFig*speciesSpace + (spPerFig - 1)*interSpecies + 0.1, y1 = 0)
	abline(h = 0)

	mtext(text = "Correlations", side = 2, outer = TRUE, las = 0)

	# Legend
	tikzAnnotate("
		\\matrix [below right] at (current bounding box.north west) {
			\\node [shape = rectangle, fill = shade_under, label = right:With competition] {}; &
			\\node [shape = rectangle, fill = shade_over, label = right:Without competition] {}; \\\\
		};
	")
	dev.off()
}

#### Demographic functions
## Growth
growth_fct = function(dbh, c0, c1, c2, scmu_g, scsd_g, scmu_dbh, scsd_dbh)
	return (exp(scsd_g * (c0 + c1*(dbh - scmu_dbh)/scsd_dbh + c2*((dbh - scmu_dbh)/scsd_dbh)^2) + scmu_g));

## Mortality
mortality_fct = function(dbh, c0, c1, c2, scmu_dbh, scsd_dbh)
	return (1 - exp(-exp(c0 + c1*(dbh - scmu_dbh)/scsd_dbh + c2*((dbh - scmu_dbh)/scsd_dbh)^2)));

#### Load data
## Coordinates of the plot
response = readRDS("../createData/treeData_presence_absence.rds")
coords = response[, .(longitude, latitude)]

## Climate database average of data over 5 years, from 2006 to 2010
database = readRDS("../createData/clim_2010.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")
setnames(database, old = names(database[, !c("longitude", "latitude")]), new = climaticVariables)
database = database[coords, on = c("latitude", "longitude"), nomatch = 0]

print(paste0("Random forest trained with ", database[, .N], " data"))

# Update the coords, loss of fiew data
coords = database[, .(longitude, latitude)]
response = response[coords, on = c("latitude", "longitude"), nomatch = 0]

## Species list
ls_species = readRDS("../createData/ls_species.rds")
n = length(ls_species)

## Load scalings
scalingGrowth = fread("../createMatlabData/growthScaling.csv");

dbh_scalingGrowth = fread("../createMatlabData/growthDbhScaling.csv");
dbh_scalingMortality = fread("../createMatlabData/mortalityDbhScaling.csv");

## data table
correlation = data.table(species = ls_species,
	correl_15cm_growth_above = numeric(length = n), correl_15cm_growth_below = numeric(length = n),
	correl_15cm_mortality_above = numeric(length = n), correl_15cm_mortality_below = numeric(length = n),
	correl_60cm_growth_above = numeric(length = n), correl_60cm_growth_below = numeric(length = n),
	correl_60cm_mortality_above = numeric(length = n), correl_60cm_mortality_below = numeric(length = n))

#### Calculate correlation
for (i in 1:n)
{
	## Species-specific variables
	species = ls_species[i]
	sdm_name = getSDMname(species)[, sdm]
	print(paste0("species: ", species))

	mu_g = scalingGrowth[species_id == species, mu];

	sd_g = scalingGrowth[species_id == species, sd];
	mu_dbh_g = dbh_scalingGrowth[species_id == species, mu];
	sd_dbh_g = dbh_scalingGrowth[species_id == species, sd];
	mu_dbh_m = dbh_scalingMortality[species_id == species, mu];
	sd_dbh_m = dbh_scalingMortality[species_id == species, sd];

	## Load the rf + vital rates models
	# SDM"s results
	sdm = readRDS(paste0("../randomForest/results/calibration/", sdm_name, "_cal.rds"))

	# Growth data
	growth_above_dt = fread(paste0("../validation/Matlab_data/", species, "/matlabGrowth_above.csv"))
	growth_below_dt = fread(paste0("../validation/Matlab_data/", species, "/matlabGrowth_below.csv"))

	# Mortality data
	mortality_above_dt = fread(paste0("../validation/Matlab_data/", species, "/matlabMortality_above.csv"))
	mortality_below_dt = fread(paste0("../validation/Matlab_data/", species, "/matlabMortality_below.csv"))

	## Prediction sdm on the plots
	# Response (boolean)
	set.seed(1969-08-18) # Woodstock seed
	pred_resp = as.logical(predict(sdm, new = database, type = "response", OOB = TRUE))

	# Probability of presence
	set.seed(1969-08-18) # Woodstock seed
	pred_prob = predict(sdm, new = database, type = "prob", OOB = TRUE)

	## Prediction growth and mortality on the plots
	# Growth
	growth_above_dt[, pred_g15 := growth_fct(150, beta0, beta1, beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g)]
	growth_below_dt[, pred_g15 := growth_fct(150, beta0, beta1, beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g)]

	growth_above_dt[, pred_g60 := growth_fct(600, beta0, beta1, beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g)]
	growth_below_dt[, pred_g60 := growth_fct(600, beta0, beta1, beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g)]

	# Mortality
	mortality_above_dt[, pred_m15 := mortality_fct(150, beta0, beta1, beta2, mu_dbh_m, sd_dbh_m)]
	mortality_below_dt[, pred_m15 := mortality_fct(150, beta0, beta1, beta2, mu_dbh_m, sd_dbh_m)]

	mortality_above_dt[, pred_m60 := mortality_fct(600, beta0, beta1, beta2, mu_dbh_m, sd_dbh_m)]
	mortality_below_dt[, pred_m60 := mortality_fct(600, beta0, beta1, beta2, mu_dbh_m, sd_dbh_m)]

	# Correlation and performance (matrix confusion)
	correlation[i, correl_15cm_growth_above := cor(pred_prob[, "TRUE"], growth_above_dt[, pred_g15])]
	correlation[i, correl_15cm_growth_below := cor(pred_prob[, "TRUE"], growth_below_dt[, pred_g15])]
	correlation[i, correl_60cm_growth_above := cor(pred_prob[, "TRUE"], growth_above_dt[, pred_g60])]
	correlation[i, correl_60cm_growth_below := cor(pred_prob[, "TRUE"], growth_below_dt[, pred_g60])]

	correlation[i, correl_15cm_mortality_above := cor(pred_prob[, "TRUE"], mortality_above_dt[, pred_m15])]
	correlation[i, correl_15cm_mortality_below := cor(pred_prob[, "TRUE"], mortality_below_dt[, pred_m15])]
	correlation[i, correl_60cm_mortality_above := cor(pred_prob[, "TRUE"], mortality_above_dt[, pred_m60])]
	correlation[i, correl_60cm_mortality_below := cor(pred_prob[, "TRUE"], mortality_below_dt[, pred_m60])]
}

#### Save the results in rds and tex format
## Add species code for tex format, sort by alphabetical order
correlation[, sp_code := getSDMname(species)[, "sp_code"]]
setorderv(correlation, "sp_code", +1)

saveRDS(correlation, "./correlation.rds")

correlation_xt = xtable(correlation[, !c("species")])

print(correlation_xt, file = "./corTab_G_mu_Pocc.tex", include.rownames = FALSE, booktabs = TRUE)

#### Plot parameters
## Colours for Northern, Middle and Southern regions
shade_under = col2rgb("#010166")[,1]
shade_over = col2rgb("#0E55FF")[,1]

## Space and width
interDot = 0.35
interSpecies = 1

speciesSpace = 4*interDot
halfSpeciesSpace = speciesSpace/2

## Load correlation R0 and SDM
correl_R0_SDM = readRDS("../validation/correlation.rds")[, .(sp_code, correl_0m, correl_10m)]

correlation = correlation[correl_R0_SDM, on = "sp_code"]
rm(correl_R0_SDM)

#### Plot
correlPlot_fct("demog_Pocc1-4.tex", rangeFig = 1:4)
correlPlot_fct("demog_Pocc5-8.tex", rangeFig = 5:8)
correlPlot_fct("demog_Pocc9-12.tex", rangeFig = 9:12)
correlPlot_fct("demog_Pocc13-14.tex", rangeFig = 13:14)
