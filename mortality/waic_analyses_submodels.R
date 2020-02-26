
#### Aim of prog: Summarise waic
## Average position:
#		- Read the waic files for each species and model
#		- Calculate their average position among species according to waic
#		- Data table for .tex file
#
## Acer saccharum (submodels):
#		- Read the waic files for each submodel + selected model (previous step)
#		- Data table for .tex file
#

#### Load packages and clear memory
library(data.table)
library(rstanarm)
library(stringi)
library(xtable)

rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool functions
## Gives the position of each model
indices = function(ref, str) # Gives NA if the model is not listed, i.e., did not satisfied some conditions
{
	n = length(ref)
	ind = integer(length = n)
	for (i in 1:n)
		ind[i] = ifelse(ref[i] %in% str, which(str == ref[i]), NA)
	dt = data.table(Modnames = ref, pos = ind)
	return (dt);
}

position = function(dt, constraint = NULL)
{
	if (!is.null(constraint))
		dt = dt[eval(parse(text = constraint))]

	dt[, mean_pos := mean(pos), by = Modnames]
	dt[, min_pos := min(pos), by = Modnames]
	dt[, max_pos := max(pos), by = Modnames]

	return(dt)
}

#### Read the data
## Common variables
(ls_folders = list.files(path = ".", pattern = "array_"))
n = length(ls_folders)
nbFormulae = 7

ls_waic = vector(mode = "list", length = n)

## Gather all the species and models' statistics together
for (i in 1:n)
{
	loadPath = paste0("./", ls_folders[i], "/")
	waic_files = list.files(path = loadPath, pattern = "waic[0-9]{1,}_submodels.rds")
	nbModels = length(waic_files) + 1 # + 1 for the selected 'mother' model
	species_txt = list.files(path = loadPath, pattern = "^[0-9]{4,}")
	species_txt = species_txt[stri_detect(str = species_txt, regex = ".txt")]
	(species = stri_sub(str = species_txt,
		to = stri_locate_last(str = species_txt, regex = ".txt")[1] - 1))
	sp_specific_data = data.table(species_id = rep(species, nbModels),
		waic = numeric(nbModels), se_waic = numeric(nbModels),
		p_waic = numeric(nbModels), se_p_waic = numeric(nbModels),
		model_id = character(nbModels))

	# Load the submodels
	for (j in 1:(nbModels - 1))
	{
		waic_vals = readRDS(paste0(loadPath, waic_files[j]))$estimates
		sp_specific_data[j, c("waic", "se_waic") := .(waic_vals["waic", "Estimate"], waic_vals["waic", "SE"])]
		sp_specific_data[j, c("p_waic", "se_p_waic") := .(waic_vals["p_waic", "Estimate"], waic_vals["p_waic", "SE"])]
	}

	# Load the 'mother' model, from which all the submodels are derived
	waic_vals = readRDS(paste0(loadPath, "waic7.rds"))$estimates
	sp_specific_data[nbModels, c("waic", "se_waic") := .(waic_vals["waic", "Estimate"], waic_vals["waic", "SE"])]
	sp_specific_data[nbModels, c("p_waic", "se_p_waic") := .(waic_vals["p_waic", "Estimate"], waic_vals["p_waic", "SE"])]

	# Create species- and model- specific names
	sp_specific_data[, model_id := paste0(i, "_model", 1:nbModels)]
	sp_specific_data[, Modnames := paste0("model", 1:nbModels)]

	# Perform the delta WAIC and sort by increasing waic (i.e., from best to worst) within species
	sp_specific_data[, delta_waic := waic - min(waic)]
	setorderv(sp_specific_data, "delta_waic", +1)

	# Store data in list
	ls_waic[[i]] = sp_specific_data
}

ls_waic = rbindlist(ls_waic)

ls_waic_pos = vector(mode = "list", length = n)

#### Calculate average position
## Position of each model within species
for (sp in 1:n)
{
	current_spModnames = paste0(sp, "_", paste0("model", 1:nbModels))
	dt = ls_waic[model_id %in% current_spModnames]
	dt_waic = setorderv(dt, "waic", +1)
	ls_waic_pos[[sp]] = indices(paste0("model", 1:nbModels), dt_waic[, Modnames])
	ls_waic_pos[[sp]] = ls_waic_pos[[sp]][dt_waic, on = "Modnames"]
}

waic_pos = rbindlist(ls_waic_pos)

## Average position (data table modified by reference!)
waic_pos = position(waic_pos) # , constraint = "se_waic < 800"
waic_pos[, modelCounter := .N, by = Modnames]

average_waic_pos = setorderv(unique(waic_pos[, c("Modnames", "mean_pos", "modelCounter")]), "mean_pos", +1)
average_waic_pos[, rel_pos := mean_pos/modelCounter]

#### Performance, 'mother' model chosen, which corresponds to model 8 in this data.table
ls_waic[Modnames == "model8", !c("model_id", "Modnames")]

#### Save results
## Format: .rds
saveRDS(average_waic_pos, "./average_waic_pos-submodels.rds")

## Format: .tex (table used for article)
waic_average_pos_xt = xtable(average_waic_pos[, .(Modnames, mean_pos, rel_pos)])
ls_waic = ls_waic[species_id == "28731-ACE-SAC", c("Modnames", "delta_waic")]
waic_acsa_xt = xtable(ls_waic)

print(waic_average_pos_xt, file = "./average_waic_pos-submodels.tex", include.rownames = FALSE, booktabs = TRUE)
print(waic_acsa_xt, file = "./waic_acsa-submodels.tex", include.rownames = FALSE, booktabs = TRUE)
