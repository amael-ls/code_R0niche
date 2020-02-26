
#### Aims of prog: Validate the PDEs model on the plots
# Correlation between R0 and the data from permanent plots
#

#### Load packages
library(data.table)
library(stringi)
library(xtable)

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

#### Load data and common variables
## Coordinates of the plot
presAbs_data = readRDS("../createData/treeData_presence_absence.rds")
coords = presAbs_data[, .(longitude, latitude)]

## Climate database used to compute R0 (useful to subset presAbs_data on the same coordinates than R0)
climData = readRDS("../createData/clim_2010.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")
setnames(climData, old = names(climData[, !c("longitude", "latitude")]), new = climaticVariables)
climData = climData[coords, on = c("latitude", "longitude"), nomatch = 0]

print(paste0("R0 computed on ", climData[, .N], " data"))

# Update the coords, lose of fiew data (but match R0 data)
coords = climData[, .(longitude, latitude)]
presAbs_data = presAbs_data[coords, on = c("latitude", "longitude"), nomatch = 0]

## Species list
ls_species = readRDS("../createData/ls_species.rds")
n = length(ls_species)

## data table results
correlation = data.table(species = ls_species,
	correlR0_presAbs_0m = numeric(length = n),
	correlR0_presAbs_10m = numeric(length = n))

#### Calculate correlation
for (i in 1:n)
{
	## Species-specific variables
	species = ls_species[i]
	sdmName = getSDMname(species)[, sdm]
	print(paste0("species: ", species))

	if (species != correlation[i, species])
		print(paste0("*** Error: Species mismatch for ", species))

	## Load Matlab's results R0
	R0_0m = unlist(fread(paste0("results/", species, "/R0_0m.csv")))
	R0_10m = unlist(fread(paste0("results/", species, "/R0_10m.csv")))

	correlation[i, correlR0_presAbs_0m := cor(presAbs_data[, ..sdmName], R0_0m)[1]]
	correlation[i, correlR0_presAbs_10m := cor(presAbs_data[, ..sdmName], R0_10m)[1]]
}

#### Save the results in rds and tex format
## Add species code for tex format, sort by alphabetical order
correlation[, sp_code := getSDMname(species)[, "sp_code"]]
setorderv(correlation, "sp_code", +1)

saveRDS(correlation, "./correlation_R0_presAbsData.rds")

correlation_xt = xtable(correlation[, .(sp_code, correlR0_presAbs_0m, correlR0_presAbs_10m)])

print(correlation_xt, file = "./validation_R0_presAbsData.tex", include.rownames = FALSE, booktabs = TRUE)
