
#### Aims of prog: To select the best model for growth
## To do it, we follow Zuur et. al (2009, section 5.7 p. 121-122)
# 	- Create the beyond optimal model (OM), ReML = TRUE, to find best random structure
#	- Then search for best fixed structure, ReML = FALSE
# 	- Rerun the selected model with ReML = TRUE

#### Load packages
library(data.table)
library(AICcmodavg)
library(doParallel)
library(stringi)
library(MuMIn)
library(lme4)
library(car)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

source("../toolFunctions.R")

#### Common variables
## Number of formulae
nbFormulae = 29

## Indices of switch for formula_id
limBeyond_OM = 4
limReML = 20

#### Create the cluster
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

(formula_id = ifelse(array_id %% nbFormulae == 0, nbFormulae, array_id %% nbFormulae))
(species_id = array_id %/% nbFormulae + ifelse(array_id %% nbFormulae == 0, 0, 1))

nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

#### Create saving folder
(savePath = paste0("array_", species_id, "/"))
if (!dir.exists(savePath))
	dir.create(savePath)

#### Load data and transform them
growth_data = readRDS("../createData/growth_dt.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")

## Subset for the species of interest
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

(species = ls_14species[species_id])
growth_data = growth_data[species_id == species]

## Write species name in the folder in case of
if (!file.exists(paste0(savePath, species, ".txt")))
	write(x = species, file = paste0(savePath, species, ".txt"))

## Log the growth, and normalise log(growth) and explanatory variables
growth_data[, growth := log(growth)]
growth_data_norm = normalisation(df = growth_data,
	colnames = c("growth", "dbh", climaticVariables),
	filename = paste0(savePath, "normalisation_growth_data_log.rds"))

range_df(growth_data_norm)

## Convert plot_id, species_id and canopy_status as a factor
growth_data_norm[, plot_id := as.factor(plot_id)]
growth_data_norm[, species_id := as.factor(species_id)]
growth_data_norm[, canopy_status := as.factor(canopy_status)]
growth_data_norm[, year_measured := as.factor(year_measured)]

#### Look for best random structure of the beyond OM
## Run the three models (same fixed structure, different random, ReML = TRUE)
if (formula_id <= limBeyond_OM)
{
	source("./glmm_beyondOM.R")

	# Selection of the random structure by ΔAICc
	ls_files_OM = list.files(path = savePath, pattern = "OM.rds")
	ls_files_OM = Filter(function(x) grepl("model", x), ls_files_OM)
	n = length(ls_files_OM)
	if (n == limBeyond_OM)
	{
		model_ls = vector("list", length = n)

		OM_names = stri_sub(ls_files_OM, to = stri_locate_first(ls_files_OM, regex = "_")[,1] - 1)

		for (i in 1:n)
			model_ls[[i]] = readRDS(paste0(savePath, ls_files_OM[i]))

		aic = AICcmodavg::aictab(model_ls, modnames = OM_names)
		r_squared = lapply(model_ls, MuMIn::r.squaredGLMM)
		r_squaredLR = lapply(model_ls, MuMIn::r.squaredLR)

		saveRDS(aic, paste0(savePath, "aic_OM.rds"))
		saveRDS(r_squared, paste0(savePath, "rsq_OM.rds"))
		saveRDS(r_squaredLR, paste0(savePath, "rsq_OM_LR.rds"))

		print(aic)
	}
}

## Get AICc for each species, check if it is always the same best model
# Everything is done in the program AICc_analyses.R. Model3 is the best model for
# all the species, Model2 is the second best model.

#### Selection of the fixed structure, ReML = FALSE
## Search for best climate variable
if (formula_id > limBeyond_OM & formula_id <= limReML)
{
	formula_id = formula_id - limBeyond_OM
	source("./glmm_ReML=FALSE.R")
}
# Selected model: 1 (with K and PP), see remark at the end of the program

# ## Alternative models, nested fixed struct, climate var (K, PP) & random struct settled
if (formula_id > limReML)
{
	formula_id = formula_id - limReML
	source("./glmm_fixef.R")
}

# Model chosen: formula2 from glmm_ReML=FALSE.R, no satisfying alternative

#### AICc in the article
## To create the tex files, I used the program final_rsq_AICc.R
