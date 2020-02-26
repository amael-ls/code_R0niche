
#### Aims of prog: Validate the PDEs model on the plots
# Correlation between R0 and the random forest prediction (on the plots)
#

#### Load packages
library(randomForest)
library(data.table)
library(stringi)
library(xtable)
# library(ROCR)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool function
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

#### Create directory
if (!dir.exists("./graphs"))
	dir.create("./graphs")

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

# Update the coords, lose of fiew data
coords = database[, .(longitude, latitude)]
response = response[coords, on = c("latitude", "longitude"), nomatch = 0]

## Species list
ls_species = readRDS("../createData/ls_species.rds")
n = length(ls_species)

## data table
correlation = data.table(species = ls_species,
	correl_0m = numeric(length = n), correl_10m = numeric(length = n),
	false_error = numeric(length = n),
	true_error = numeric(length = n),
	R2_tjur = numeric(length = n))

#### Calculate correlation
for (i in 1:n)
{
	## Species-specific variables
	species = ls_species[i]
	sdm_name = getSDMname(species)[, sdm]
	print(paste0("species: ", species))

	## Load the two models
	# SDM's results
	sdm = readRDS(paste0("../randomForest/results/calibration/", sdm_name, "_cal.rds"))

	# Matlab's results R0
	R0 = unlist(fread(paste0("results/", species, "/R0_10m.csv")))
	R0_0m = unlist(fread(paste0("results/", species, "/R0_0m.csv")))

	## Prediction sdm on the plots
	# Response (boolean)
	set.seed(1969-08-18) # Woodstock seed
	pred_resp = as.logical(predict(sdm, new = database, type = "response", OOB = TRUE))

	# Probability of presence
	set.seed(1969-08-18) # Woodstock seed
	pred_prob = predict(sdm, new = database, type = "prob", OOB = TRUE)

	# Correlation and performance (matrix confusion)
	correlation[i, correl_0m := cor(pred_prob[, "TRUE"], R0_0m)]
	correlation[i, correl_10m := cor(pred_prob[, "TRUE"], R0)]
	correlation[i, false_error := sdm["confusion"][[1]]["FALSE", "class.error"]]
	correlation[i, true_error := sdm["confusion"][[1]]["TRUE", "class.error"]]

	## Plots
	# Response (boolean)
	jpeg(paste0("graphs/", species, "_resp.jpg"))
	par(mar = c(1.5, 1.5, 0.5, 0.5), mgp = c(1.5, 0.3, 0), tck = -0.015)
	plot(R0, pred_resp)
	dev.off()

	# Probability of presence
	jpeg(paste0("graphs/", species, "_prob.jpg"))
	par(mar = c(1.5, 1.5, 0.5, 0.5), mgp = c(1.5, 0.3, 0), tck = -0.015)
	plot(R0, pred_prob[, "TRUE"])
	dev.off()

	## R² Tjur 2009
	# Create data
	R2_data = response[, ..sdm_name]
	R2_data[, prediction := pred_prob[, "TRUE"]]

	# Change the name of the species to a generic name (quick an dirty cheating solution)
	setnames(R2_data, old = sdm_name, new = "focus_sp")

	# Calculate R² Tjur
	R2 = R2_data[(focus_sp), mean(prediction)] - R2_data[!(focus_sp), mean(prediction)]

	correlation[i, R2_tjur := R2]
}

#### Save the results in rds and tex format
## Add species code for tex format, sort by alphabetical order
correlation[, sp_code := getSDMname(species)[, "sp_code"]]
setorderv(correlation, "sp_code", +1)

saveRDS(correlation, "./correlation.rds")

correlation_xt = xtable(correlation[, .(sp_code, correl_0m, correl_10m, R2_tjur, true_error, false_error)])

print(correlation_xt, file = "./validation.tex", include.rownames = FALSE, booktabs = TRUE)
