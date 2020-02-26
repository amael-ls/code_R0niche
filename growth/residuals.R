
#### Aim of prog: GLMMs assume normality of the residuals among others, we check them
# Check https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#calculating-scaled-residuals
# to analyse the graphical results

#### Load package and clear memory
library(data.table)
library(doParallel)
library(merTools)
library(DHARMa)
library(lme4)

rm(list = ls())
graphics.off()
options(max.print = 500)

# source("../myPlotFEsim.R")
source("../toolFunctions.R")

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

#### Saving folder and load model
(savePath = paste0("./array_", array_id, "/"))
if (!dir.exists(savePath))
	print(paste0("*** ERROR: folder ", savePath, " does not exist"))

model = readRDS(paste0(savePath, "final_model.rds"))

#### Load data and transform them
growth_data = readRDS("../createData/growth_dt.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")

## Subset for the species of interest /!\ Same order as growth_selection /!\
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

(species = ls_14species[array_id])
growth_data = growth_data[species_id == species]

## Check species name in the folder (cf growth_selection.R) to check it corresponds to the good species
txt_file = list.files(path = savePath, pattern = paste0(species, ".txt"))
if (length(txt_file) != 1)
	print(paste0("*** ERROR: check species name ", species, " versus folder id ", array_id))

## Log the growth, and normalise log(growth) and explanatory variables
growth_data[, growth := log(growth)]
growth_data_norm = normalisation(df = growth_data,
	colnames = c("growth", "dbh", climaticVariables),
	filename = "./rubbishFromResiduals.rds")

range_df(growth_data_norm)

## Convert plot_id, species_id and canopy_status as a factor
growth_data_norm[, plot_id := as.factor(plot_id)]
growth_data_norm[, species_id := as.factor(species_id)]
growth_data_norm[, canopy_status := as.factor(canopy_status)]
growth_data_norm[, year_measured := as.factor(year_measured)]

#### Check the residuals of the selected model
simu = simulateResiduals(model, n = 2500, refit = FALSE, integerResponse = FALSE, plot = FALSE)
saveRDS(simu, paste0(savePath, "resTot_final_model.rds"))

testU = testUniformity(simu, plot = FALSE) # alternative = c("two.sided", "less", "greater")
saveRDS(testU, paste0(savePath, "testUniformity_final_model.rds"))

resPerGroup_plot = recalculateResiduals(simu, group = growth_data_norm$plot_id)
saveRDS(resPerGroup_plot, paste0(savePath, "resPerPlot_final_model.rds"))

resPerGroup_year = recalculateResiduals(simu, group = growth_data_norm$year_measured)
saveRDS(resPerGroup_year, paste0(savePath, "resPerYear_final_model.rds"))

jpeg(paste0(savePath, "res_finalModel.jpg"), width = 1000, height = 1000, quality = 100)
plot(simu)
dev.off()

jpeg(paste0(savePath, "testUnif_finalModel.jpg"), width = 1000, height = 1000, quality = 100)
plotQQunif(simu)
dev.off()

jpeg(paste0(savePath, "resRecalcPlot_finalModel.jpg"), width = 1000, height = 1000, quality = 100)
plot(resPerGroup_plot)
dev.off()

jpeg(paste0(savePath, "resRecalcYear_finalModel.jpg"), width = 1000, height = 1000, quality = 100)
plot(resPerGroup_year)
dev.off()

#### Confidence intervals
fe = FEsim(model, n.sims = 500, oddsRatio = FALSE, seed = NULL)
saveRDS(fe, paste0(savePath, "fe_final_model.rds"))

myPlotFEsim(data = fe, level = 0.95, stat = "median", sd = TRUE,
	filenameGraph = paste0(savePath, "./FEsim_final_model.pdf"), nbParamsPerSheet = 20)

#### Summary
sumUp = merTools::fastdisp(model)
print(sumUp)
saveRDS(sumUp, paste0(savePath, "sumUp.rds"))
