
#### Aim of prog: treating AICc and R² for every species
## Check the first three best random structure
#		- Load AICc OM
#		- Create a data table species x best 3 models
#
## Check the first three best fixed structure (random structure already selected)
#		- Load AICc OM
#		- Create a data table species x best 3 models
#

#### Load package and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(AICcmodavg)
library(stringi)
library(xtable)
library(MuMIn)

#### Tool functions
extract_3Models = function(array_id, which_aic)
{
	loadPath = paste0("./array_", array_id, "/")
	aic = readRDS(paste0(loadPath, "aic_", which_aic, ".rds"))
	setDT(aic)
	setorderv(x = aic, cols = c("Delta_AICc"))
	aic_names = as.list(as.character(aic[1:3, Modnames]))
	aic_values = as.list(aic[1:3, Delta_AICc])
	return (list(names = aic_names, deltaAICc = aic_values))
}

extract_3Models_v2 = function(array_id, qty = 3, which_criteria = c("R2_m", "Delta_AICc"),
	condition = NULL, subsetModels = NULL)
{
	loadPath = paste0("./array_", array_id, "/")
	dt = readRDS(paste0(loadPath, "rsq_vif_aic.rds"))

	if (!is.null(subsetModels))
		dt = dt[Modnames %in% subsetModels]

	if (!is.null(condition))
		dt = dt[eval(parse(text = condition))]

	setorderv(x = dt, cols = c(which_criteria), ifelse(which_criteria == "R2_m", -1, 1))

	dt[, sp_Modnames := paste0(array_id, "_", Modnames)]
	return (dt[1:qty])
}

indices = function(ref, str) # Gives NA if the model is not listed, i.e., did not satisfied some conditions
{
	n = length(ref)
	ind = integer(length = n)
	for (i in 1:n)
		ind[i] = ifelse(ref[i] %in% str, which(str == ref[i]), NA)
	dt = data.table(Modnames = ref, pos = ind)
	return (dt);
}

rel_mean = function(pos, freq_rsq, na.rm = TRUE)
{
	freq_rsq = unname(freq_rsq)
	relMean = mean(pos, na.rm = ifelse(na.rm, TRUE, FALSE))/freq_rsq
	return (relMean);
}

#### Create species x best 3 models for random structure
## List of the species /!\ Must be in the same order than growth_selection.R /!\
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

n = length(ls_14species)
randomStruct = data.table(array_id = 1:n, species_id = ls_14species,
	first = character(n), second = character(n), third = character(n))

## Random structure
randomStruct[, c("first", "second", "third") := extract_3Models(array_id, "OM")[["names"]], by = array_id]
randomStruct[, table(first)]
randomStruct[, table(second)]
randomStruct[, table(third)]

## Create the .tex file to select best random structure
speciesTable = readRDS("../createData/speciesTable.rds")
randomStruct = randomStruct[speciesTable[, .(species, species_id)], on = "species_id"]

randomStruct_xt = xtable(randomStruct[, .(species, first, second, third)])

print(randomStruct_xt, file = "./randomStruct.tex", include.rownames = FALSE, booktabs = TRUE)

# I chose the best model, which is the same for all the species as mentionned in the GLMM suppl. inf.

#### Selection of the climatic variables, taking VIF into account
## Merge rsq, vif and aic from: ReML=FALSE (clim variables) and fixef (submodels)
for (sp in 1:n)
{
	rsq_i = readRDS(paste0("array_", sp, "/rsq_mega.rds"))
	vif_i = readRDS(paste0("array_", sp, "/vif_mega.rds"))
	aic_i = setDT(readRDS(paste0("array_", sp, "/aic_mega.rds")))

	vif_i = vif_i[rsq_i, on = "Modnames"]
	vif_i = vif_i[aic_i[, .(Modnames, Delta_AICc)], on = "Modnames", nomatch = 0]

	setorderv(vif_i, "Delta_AICc", +1)
	saveRDS(vif_i, paste0("./array_", sp, "/rsq_vif_aic.rds"))
}

#### Among ReML models (to select best climate predictors) based on R² and then AIC
## Criteria = R2_m
# Create list to store the results for all the species
liszt_rsq = vector(mode = "list", length = n)
liszt_rsq20 = vector(mode = "list", length = n)

## Get models' name from ReML (sp = n, not a problem, equal nb of models among species)
loadPath = paste0("./array_", sp, "/")
ls_files_ReML = list.files(path = loadPath, pattern = "ReML=FALSE_vif.rds")
ReML_names = stri_sub(ls_files_ReML,
	to = stri_locate_first(ls_files_ReML, regex = "_")[,1] - 1)
n_ReML = length(ls_files_ReML)

## All the models based on marginal R², max vif < 20
for (i in 1:n)
	liszt_rsq20[[i]] = extract_3Models_v2(i, n_ReML, "R2_m", "max < 20", ReML_names)

# Test taking all the models
for (i in 1:n)
	liszt_rsq[[i]] = extract_3Models_v2(i, n_ReML, "R2_m", "max < Inf", ReML_names)

liszt_rsq = rbindlist(liszt_rsq)
liszt_rsq20 = rbindlist(liszt_rsq20)

## All the models based on Delta_AICc, max vif < 20
liszt_aic = vector(mode = "list", length = n)
liszt_aic20 = vector(mode = "list", length = n)

for (i in 1:n)
	liszt_aic20[[i]] = extract_3Models_v2(i, n_ReML, "Delta_AICc", "max < 20", ReML_names)

# Test taking all the models
for (i in 1:n)
	liszt_aic[[i]] = extract_3Models_v2(i, n_ReML, "Delta_AICc", "max < Inf", ReML_names)

liszt_aic = rbindlist(liszt_aic)
liszt_aic20 = rbindlist(liszt_aic20)

## Frequences
freq_rsq = table(liszt_rsq[, Modnames]) # Should be the number of species
freq_aic = table(liszt_aic[, Modnames]) # Should be the number of species

freq_rsq20 = table(liszt_rsq20[, Modnames])
freq_aic20 = table(liszt_aic20[, Modnames])

#### Calculate the average position for each model
## List all the models
rsq_models = unique(liszt_rsq[!is.na(Modnames), Modnames])
aic_models = unique(liszt_aic[!is.na(Modnames), Modnames])

ls_rsq_pos = vector(mode = "list", length = n)
ls_aic_pos = vector(mode = "list", length = n)

for (sp in 1:n)
{
	current_spModnames = paste0(sp, "_", rsq_models)
	dt = liszt_rsq[sp_Modnames %in% current_spModnames]
	dt_rsq = setorderv(dt[, .(Modnames, R2_m)], "R2_m", -1)
	ls_rsq_pos[[sp]] = indices(rsq_models, dt_rsq[, Modnames])

	dt_aic = setorderv(dt[Modnames %in% aic_models, .(Modnames, Delta_AICc)], "Delta_AICc", +1)
	ls_aic_pos[[sp]] = indices(aic_models, dt_aic[, Modnames])
}

ls_rsq_pos = rbindlist(ls_rsq_pos)
ls_aic_pos = rbindlist(ls_aic_pos)

## List all the models for which at least there is one vif < threshold (= 20 in my case)
rsq_models20 = unique(liszt_rsq20[!is.na(Modnames), Modnames])
aic_models20 = unique(liszt_aic20[!is.na(Modnames), Modnames])

ls_rsq_pos20 = vector(mode = "list", length = n)
ls_aic_pos20 = vector(mode = "list", length = n)

for (sp in 1:n)
{
	current_spModnames = paste0(sp, "_", rsq_models)
	dt = liszt_rsq20[sp_Modnames %in% current_spModnames]
	dt_rsq = setorderv(dt[, .(Modnames, R2_m)], "R2_m", -1)
	ls_rsq_pos20[[sp]] = indices(rsq_models20, dt_rsq[, Modnames])

	dt_aic = setorderv(dt[Modnames %in% aic_models, .(Modnames, Delta_AICc)], "Delta_AICc", +1)
	ls_aic_pos20[[sp]] = indices(aic_models20, dt_aic[, Modnames])
}

ls_rsq_pos20 = rbindlist(ls_rsq_pos20)
ls_aic_pos20 = rbindlist(ls_aic_pos20)

## Mean position R² based, na.rm = FALSE. Set freq = 1 to have non relative positions
ls_rsq_pos[, meanPos := rel_mean(pos, 1, na.rm = FALSE), by = Modnames]
rsq_average_pos = setorderv(unique(ls_rsq_pos[, .(Modnames, meanPos)]), "meanPos", +1)
print(rsq_average_pos)

## Mean position Δ AICc based
ls_aic_pos[, meanPos := rel_mean(pos, 1, na.rm = FALSE), by = Modnames]
aic_average_pos = setorderv(unique(ls_aic_pos[, .(Modnames, meanPos)]), "meanPos", +1)
print(aic_average_pos)

## Mean position R² based, na.rm = FALSE, vif < 20
ls_rsq_pos20[, meanPos := rel_mean(pos, freq_rsq20[Modnames], na.rm = FALSE), by = Modnames]
rsq_average_pos20 = setorderv(unique(ls_rsq_pos20[, .(Modnames, meanPos)]), "meanPos", +1)
print(rsq_average_pos20)

## Mean position Δ AICc based
ls_aic_pos20[, meanPos := rel_mean(pos, freq_aic20[Modnames], na.rm = FALSE), by = Modnames]
aic_average_pos20 = setorderv(unique(ls_aic_pos20[, .(Modnames, meanPos)]), "meanPos", +1)
print(aic_average_pos20)

#### Get the VIF for all the species/models and average them for each model
ls_vif = vector(mode = "list", length = n)
delMod = paste0("model", 17:25)
for (sp in 1:n)
{
	ls_vif[[sp]] = readRDS(paste0("array_", sp, "/vif_mega.rds"))
	ls_vif[[sp]] = ls_vif[[sp]][!(Modnames %in% delMod)]
}

ls_vif = rbindlist(ls_vif)
ls_vif[, meanMaxVIF := mean(max), by = Modnames]
ls_vif[, maxMaxVIF := max(max), by = Modnames]
vif_average = unique(ls_vif[, .(Modnames, meanMaxVIF, maxMaxVIF)])

#### Save .tex file climate structure
rsq_average_pos = na.omit(rsq_average_pos) # Should not have any NA
aic_average_pos = na.omit(aic_average_pos) # Should not have any NA

rsq_average_pos = rsq_average_pos[aic_average_pos, on = "Modnames"]
rsq_average_pos = rsq_average_pos[vif_average, on = "Modnames"]

setnames(rsq_average_pos, old = colnames(rsq_average_pos),
	new = c("Model", "meanPosR2", "meanPosAICc", "meanMaxVIF", "maxMaxVIF"))

setorderv(rsq_average_pos, "meanPosR2")
rsq_average_xt = xtable(rsq_average_pos)
print(rsq_average_xt, file = "./averagePosClim.tex", include.rownames = FALSE, booktabs = TRUE)

# Model 12 is the only one that respect for all the species max(vif) < 6.
# It is, in average, the 5th best model in AIC and R² (mean = 4.64 and 5.21 respectively)
# However, model 1 and 2 are much better, although their vif are less good. The model 2 also has all the max vif below 20 (6 in average)
# while model 1 has one vif at 24 (11 in average). Model 2 is also much simpler.
# Hence I CHOSE model 2 for its simplicity, interpretability and performance.

#### And then between model 2 (ReML) and its submodels (fixef), to show the impact of climate
## Criteria = R2_m
liszt_rsq = vector(mode = "list", length = n)

loadPath = paste0("./array_", sp, "/")
ls_files_fixef = list.files(path = loadPath, pattern = "fixef_vif.rds")
fixef_names = stri_sub(ls_files_fixef,
	to = stri_locate_first(ls_files_fixef, regex = "_")[,1] - 1)
rplc = as.integer(stri_sub(fixef_names,
	from = stri_locate_first(fixef_names, regex = "\\d")[,1])) + n_ReML

fixef_names = stri_replace_all(str = fixef_names, replacement = rplc, regex = "\\d")

addModels = c(fixef_names, "model2", "model19", "model21") # no vif for 19 and 21, only one variable!
for (i in 1:n)
	liszt_rsq[[i]] = extract_3Models_v2(i, 10, "R2_m", subsetModels = addModels)

liszt_rsq = rbindlist(liszt_rsq)

## Criteria = Delta_AICc
liszt_aic = vector(mode = "list", length = n)
for (i in 1:n)
	liszt_aic[[i]] = extract_3Models_v2(i, 10, "Delta_AICc", subsetModels = addModels)

liszt_aic = rbindlist(liszt_aic)

## Frequence
freq_rsq = table(liszt_rsq[, Modnames])
freq_aic = table(liszt_aic[, Modnames])

#### Calculate the average position for each model
rsq_models = unique(liszt_rsq[!is.na(Modnames), Modnames])
aic_models = unique(liszt_aic[!is.na(Modnames), Modnames])

ls_rsq_pos = vector(mode = "list", length = n)
ls_aic_pos = vector(mode = "list", length = n)

for (sp in 1:n)
{
	current_spModnames = paste0(sp, "_", rsq_models)
	dt = liszt_rsq[sp_Modnames %in% current_spModnames]
	dt_rsq = setorderv(dt[, .(Modnames, R2_m)], "R2_m", -1)
	ls_rsq_pos[[sp]] = indices(rsq_models, dt_rsq[, Modnames])

	dt_aic = setorderv(dt[Modnames %in% aic_models, .(Modnames, Delta_AICc)], "Delta_AICc", +1)
	ls_aic_pos[[sp]] = indices(aic_models, dt_aic[, Modnames])
}

ls_rsq_pos = rbindlist(ls_rsq_pos)
ls_aic_pos = rbindlist(ls_aic_pos)

## Mean position R² based
ls_rsq_pos[, meanPos := rel_mean(pos, freq_rsq[Modnames], na.rm = TRUE), by = Modnames]
rsq_average_pos = setorderv(unique(ls_rsq_pos[, .(Modnames, meanPos)]), "meanPos", +1)
print(rsq_average_pos)

## Mean position Δ AICc based
ls_aic_pos[, meanPos := rel_mean(pos, freq_aic[Modnames], na.rm = TRUE), by = Modnames]
aic_average_pos = setorderv(unique(ls_aic_pos[, .(Modnames, meanPos)]), "meanPos", +1)
print(aic_average_pos)

## Mean position R² based
ls_rsq_pos[, meanPos := rel_mean(pos, 1, na.rm = FALSE), by = Modnames]
rsq_average_pos = setorderv(unique(ls_rsq_pos[, .(Modnames, meanPos)]), "meanPos", +1)
print(rsq_average_pos)

## Mean position Δ AICc based
ls_aic_pos[, meanPos := rel_mean(pos, 1, na.rm = FALSE), by = Modnames]
aic_average_pos = setorderv(unique(ls_aic_pos[, .(Modnames, meanPos)]), "meanPos", +1)
print(aic_average_pos)

# Here in both cases, model 2 is the best

#### Sum-up perfomance of model 12
liszt_rsq[Modnames == "model2"]
liszt_aic[Modnames == "model2"]

ls_rsq_pos[Modnames == "model2", .(Modnames, pos)]
ls_aic_pos[Modnames == "model2", .(Modnames, pos)]

#### Table for ACSA
## Get the models
acsa_dt = liszt_rsq[sp_Modnames %in% paste0("3_model", c(2, 17:25))]

## Rescale the Delta_AICc on 2nd model (rather than the 7th that have insane VIF)
acsa_dt[, Delta_AICc := Delta_AICc - min(Delta_AICc)]

## Keep only columns of interest
acsa_dt[, c("max", "mean", "sp_Modnames") := NULL]

## Sort by Δ AICc
setorderv(acsa_dt, "Delta_AICc", +1)

## Save the .tex file
acsa_xt = xtable(acsa_dt)

print(acsa_xt, file = "./acsa.tex", include.rownames = FALSE, booktabs = TRUE)
