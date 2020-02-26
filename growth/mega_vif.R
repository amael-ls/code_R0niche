
####
library(data.table)
library(stringi)

#### Tool function
vif_combine = function(array_id)
{
	loadPath = paste0("./array_", array_id, "/")
	count = 1

	# List files
	ls_files_ReML = list.files(path = loadPath, pattern = "ReML=FALSE_vif.rds")
	ls_files_ReML = Filter(function(x) grepl("model", x), ls_files_ReML)
	n_ReML = length(ls_files_ReML)

	ls_files_fixef = list.files(path = loadPath, pattern = "fixef_vif.rds")
	ls_files_fixef = Filter(function(x) grepl("model", x), ls_files_fixef)
	n_fixef = length(ls_files_fixef)

	n = n_ReML + n_fixef
	ls_vif = vector(mode = "list", length = n)

	# Read files
	for (file in ls_files_ReML)
	{
		current_vif = readRDS(paste0(loadPath, file))
		model = stri_sub(str = file, to = stri_locate_first(str = file, regex = "_")[1] - 1)
		ls_vif[[count]] = data.table(Modnames = model,
			max = max(current_vif), mean = mean(current_vif))
		count = count + 1
	}

	for (file in ls_files_fixef)
	{
		current_vif = readRDS(paste0(loadPath, file))
		model = stri_sub(str = file, to = stri_locate_first(str = file, regex = "_")[1] - 1)

		rplc = as.integer(stri_sub(model,
			from = stri_locate_first(model, regex = "\\d")[,1])) + n_ReML

		model = stri_replace_all(str = model, replacement = rplc, regex = "\\d")

		ls_vif[[count]] = data.table(Modnames = model,
			max = max(current_vif), mean = mean(current_vif))
		count = count + 1
	}
	ls_vif = rbindlist(ls_vif)
	return(ls_vif)
}

for (array_id in 1:14)
{
	mega_vif = vif_combine(array_id)
	saveRDS(mega_vif, paste0("./array_", array_id, "/vif_mega.rds"))
}


#### vif giga!
ls_vif = vector(mode = "list", length = 14)

for (array_id in 1:14)
{
	ls_vif[[array_id]] = readRDS(paste0("./array_", array_id, "/vif_mega.rds"))
	ls_vif[[array_id]][, Modnames := paste0(Modnames, "_", array_id)]
}

vif_giga = rbindlist(ls_vif)
vif_giga[stri_detect(str = Modnames, regex = "model2_")]
