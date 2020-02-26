
#### Aims of prog: To select the best model for growth (last step)
## To do it, we follow Zuur et. al (2009, section 5.7 p. 121-122)
# 	- Create the beyond optimal model (OM), ReML = TRUE, to find best random structure
#	- Then search for best fixed structure, ReML = FALSE
# 	- Rerun the selected model with ReML = TRUE

#### Load packages
library(data.table)
library(doParallel)
library(stringi)
library(MuMIn)
library(lme4)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

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

#### Create saving folder
(savePath = paste0("./array_", array_id, "/"))
if (!dir.exists(savePath))
	print(paste0("*** ERROR: folder ", savePath, " does not exist"))

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
	filename = paste0(savePath, "normalisation_growth_data_log_final.rds"))

range_df(growth_data_norm)

## Convert plot_id, species_id and canopy_status as a factor
growth_data_norm[, plot_id := as.factor(plot_id)]
growth_data_norm[, species_id := as.factor(species_id)]
growth_data_norm[, canopy_status := as.factor(canopy_status)]
growth_data_norm[, year_measured := as.factor(year_measured)]

#### Run final model, ReML = TRUE
## Formula = model 2, cf glmm_ReML=FALSE.R
final_formula = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(annual_mean_temperature + I(annual_mean_temperature^2) +
		annual_precipitation + I(annual_precipitation^2))"

## Run this time with ReML = TRUE
model = lmer(formula = as.formula(final_formula), data = growth_data_norm,
	REML = TRUE)

## Save
saveRDS(model, paste0(savePath, "final_model.rds"))

#### R² marginal and conditional
r_squared = r.squaredGLMM(model)
saveRDS(r_squared, paste0(savePath, "rsq_final.rds"))

#### Get fixed effect (for matlab)
fixef_growth = getME(model, "fixef")
saveRDS(fixef_growth, paste0(savePath, "fixef_growth.rds"))
