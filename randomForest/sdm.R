
#### Aims of prog: Species Distribution Model (SDM) using a random forest
## SDM:
#	- Prepare the data (created by sdmData.R):
#		* Load climate 2010 (averaged over 5 years!)
#		* Load the presence/absence data from forest inventories (1999-2010)
#		* Merge the two dataset
#	- SDM without space:
#		* Variables: all the climate variables
#		* Calibration, predictions and graphs are saved
#

#### Load packages
library(randomForest)
library(RColorBrewer)
library(data.table)
library(doParallel)
library(stringi)
library(sf)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Create the cluster
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

#### Tool function
getSpeciesLittle = function(str)
{
	# Replace underscores by scores
	species = stri_replace_all(str = str, replacement = "-", regex = "_")

	# Switch the Taxonomic Serial Number (tsn) with species code
	tsn_pos = stri_locate_last(str = species, regex = "-")[1]
	tsn = stri_sub(str = species, from = tsn_pos + 1)
	sp_code = stri_sub(str = species, to = tsn_pos - 1)
	species = paste0(tsn, "-", sp_code)
	return (species)
}

#### Load data
dataCalibration = readRDS("../createData/treeData_presence_absence.rds")
clim_2010 = readRDS("../createData/clim_2010.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")
setnames(clim_2010, old = names(clim_2010[, !c("longitude", "latitude")]), new = climaticVariables)

## Subset for the target species (using col names)
(varToCalibrate = names(dataCalibration)[stri_detect(names(dataCalibration), regex = "_")][array_id])
dataCalibration = dataCalibration[, c("longitude", "latitude", varToCalibrate),
	with = FALSE]

## Transform the TRUE/FALSE to factor
dataCalibration[, c(varToCalibrate) := lapply(.SD, function(x) (as.factor(x))), .SDcols = varToCalibrate]

## Add climate (of 2010) to dataCalibration
dataCalibration = dataCalibration[clim_2010, on = c("latitude", "longitude"), nomatch = 0]

## Load climate on Little's 1971 map; for crs, check createMatlabData.R
# crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
(little_id = getSpeciesLittle(varToCalibrate))
dataEvaluation = readRDS(paste0("../R0/Matlab_data/", little_id, "/", little_id, "_clim2010.rds"))

#### Run calibration, prediction, and maps
## Calibration
set.seed(1969-08-18) # Woodstock seed
sdm_formula = paste0(varToCalibrate, " ~ ", paste0(climaticVariables, collapse = " + "))
SDM_cal = randomForest(as.formula(sdm_formula), mtry = 12,
	data = dataCalibration, ntree = 2000)

## Preciction on Little's 1971 map
set.seed(1969-08-18) # Woodstock seed
pred_resp = predict(SDM_cal, new = dataEvaluation, type = "response", OOB = TRUE)

set.seed(1969-08-18) # Woodstock seed
pred_prob = predict(SDM_cal, new = dataEvaluation, type = "prob", OOB = TRUE)

dataEvaluation[, c("pres", "prob_pres", "prob_abs") := .(pred_resp, pred_prob[, "TRUE"], pred_prob[, "FALSE"])]
dataEvaluation[, pres := as.logical(pres)]

#### Save results and plot predicted map distribution
## Save
saveRDS(SDM_cal, paste0("results/calibration/",varToCalibrate, "_cal.rds"))
saveRDS(dataEvaluation, paste0("results/prediction/", varToCalibrate, "_pred.rds"))

## Map probability of presence (continuous variable, between 0 and 1)
# Load shapefiles of Northern America
canada = readRDS("~/database/shapefiles/North_America_Shapefile/canada_union.rds")
usa = readRDS("~/database/shapefiles/North_America_Shapefile/usa_union.rds")

erie = st_read("~/database/shapefiles/North_America_Shapefile/lakes/erie")
huron = st_read("~/database/shapefiles/North_America_Shapefile/lakes/huron")
michigan = st_read("~/database/shapefiles/North_America_Shapefile/lakes/michigan")
ontario = st_read("~/database/shapefiles/North_America_Shapefile/lakes/ontario")
stClair = st_read("~/database/shapefiles/North_America_Shapefile/lakes/stClair")
superior = st_read("~/database/shapefiles/North_America_Shapefile/lakes/superior")

canada = st_transform(canada,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
usa = st_transform(usa,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
erie = st_transform(erie,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
huron = st_transform(huron,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
michigan = st_transform(michigan,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
ontario = st_transform(ontario,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
stClair = st_transform(stClair,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
superior = st_transform(superior,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

canada = st_geometry(canada)
usa = st_geometry(usa)
erie = st_geometry(erie)
huron = st_geometry(huron)
michigan = st_geometry(michigan)
ontario = st_geometry(ontario)
stClair = st_geometry(stClair)
superior = st_geometry(superior)

# Load the Little's 1971 shapefile
ls_little = list.files(path = paste0("../little1971/", little_id), pattern = ".shp$")
little1971 = st_read(paste0("../little1971/", little_id, "/", ls_little))
little1971 = st_transform(little1971,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
bounding_little1971 = st_bbox(little1971)

little1971 = st_geometry(little1971)

## Bounding box of Little's 1971 envelop
lonMin = unname(bounding_little1971["xmin"]) - 100
lonMax = unname(bounding_little1971["xmax"]) + 100
latMin = unname(bounding_little1971["ymin"]) - 100
latMax = unname(bounding_little1971["ymax"]) + 100

# rbPal <- colfunc<-colorRampPalette(c("royalblue", "springgreen", "yellow", "red"))
# dataEvaluation[, colourProb_pres := rbPal(30)[as.numeric(cut(dataEvaluation$prob_pres, breaks = 30))]]
dataEvaluation = st_as_sf(x = dataEvaluation, coords = c("longitude", "latitude"),
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

jpeg(paste0("results/graphs/", varToCalibrate, "_prob.jpg"), width = 1000, height = 1000, quality = 100)
plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
	xlab = "", ylab = "")
plot(dataEvaluation["prob_pres"], cex = 0.25, pch = 20, add = TRUE,
	ylab = "", pal = grey.colors(13, start = 0.9, end = 0))
plot(little1971, col = NA, border = "#DD925C", add = TRUE)
plot(canada, col = NA, add = TRUE)
plot(usa, col = NA, add = TRUE)
plot(erie, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(huron, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(michigan, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(ontario, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(stClair, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(superior, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
dev.off()

## map presence/absence (boolean)
jpeg(paste0("results/graphs/", varToCalibrate, ".jpg"), width = 1000, height = 1000, quality = 100)
plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
	xlab = "", ylab = "")
plot(dataEvaluation["pres"], cex = 0.25, pch = 20, add = TRUE,
	ylab = "", pal = grey.colors(2, start = 0.9, end = 0))
plot(little1971, col = NA, border = "#DD925C", add = TRUE)
plot(canada, col = NA, add = TRUE)
plot(usa, col = NA, add = TRUE)
plot(erie, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(huron, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(michigan, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(ontario, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(stClair, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
plot(superior, col = rgb(92/255, 167/255, 221/255, 0.4), add = TRUE)
dev.off()
