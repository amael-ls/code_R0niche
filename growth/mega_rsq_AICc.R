
#### Aim of prog: Calculate R² and AICc of ReML = FALSE and fixef models
## Comments:
# There are 5 models fitted in glmm_ReML=FALSE.R and 9 simpler alternatives fitted
# in glmm_fixef.R, also using Maximum Likelihood method (=> 14 comparable models).
#
## Calculate AICc & R²:
#		- List and load all the concerned models
#		- Calculate AICc using AICcmodavg package
#		- Calculate R² using MuMIn
#		- Save the results in rds (and tex format if uncommented)

#### Load package and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(AICcmodavg)
library(doParallel)
library(stringi)
library(xtable)
library(MuMIn)

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

(loadPath = paste0("./array_", array_id, "/"))

#### Load data
## List files
ls_files_ReML = list.files(path = loadPath, pattern = "ReML=FALSE.rds")
ls_files_ReML = Filter(function(x) grepl("model", x), ls_files_ReML)
(n_ReML = length(ls_files_ReML))

ls_files_fixef = list.files(path = loadPath, pattern = "fixef.rds")
ls_files_fixef = Filter(function(x) grepl("model", x), ls_files_fixef)
(n_fixef = length(ls_files_fixef))

n = n_ReML + n_fixef

model_ls = vector("list", length = n)

## Create models' name
ReML_names = stri_sub(ls_files_ReML,
	to = stri_locate_first(ls_files_ReML, regex = "_")[,1] - 1)

fixef_names = stri_sub(ls_files_fixef,
	to = stri_locate_first(ls_files_fixef, regex = "_")[,1] - 1)

# The digits in fixef_names must be replaced to distinguish them from ReML names
rplc = as.integer(stri_sub(fixef_names,
	from = stri_locate_first(fixef_names, regex = "\\d")[,1])) + n_ReML

fixef_names = stri_replace_all(str = fixef_names, replacement = rplc, regex = "\\d")

## Read files
for (i in 1:n_ReML)
{
	print(i)
	model_ls[[i]] = readRDS(paste0(loadPath, ls_files_ReML[i]))
}

if (n_fixef > 0)
{
	for (i in (n_ReML + 1):n)
	{
		print(i)
		model_ls[[i]] = readRDS(paste0(loadPath, ls_files_fixef[i - n_ReML]))
	}
}

#### AICc and R²
aic = AICcmodavg::aictab(model_ls, modnames = c(ReML_names, fixef_names))
r_squared = lapply(model_ls, MuMIn::r.squaredGLMM)

## Reshape r_squared, to make it easier to read
rsq = as.data.table(matrix(data = unlist(r_squared), nrow = n, ncol = 2, byrow = TRUE))
setnames(rsq, old = c("V1", "V2"), new = c("R2_m", "R2_c"))
modnames_ReML = stri_sub(str = ls_files_ReML,
	to = stri_locate_first(str = ls_files_ReML, regex = "_")[, 1] - 1)
rsq[1:n_ReML, Modnames := modnames_ReML]

if (n_fixef > 0)
{
	modnames_fixef = stri_sub(str = ls_files_fixef,
		to = stri_locate_first(str = ls_files_fixef, regex = "_")[, 1] - 1)
	rsq[(n_ReML + 1):n, Modnames := paste0("model", (n_ReML + 1):n)]
}

#### Save results
## rds format
saveRDS(aic, paste0(loadPath, "aic_mega.rds"))
saveRDS(rsq, paste0(loadPath, "rsq_mega.rds"))

# ## latex format
# # Transform to data table
# setDT(aic)
# setDT(r_squared)
# r_squared = transpose(r_squared)
# r_squared[, Modnames := c(ReML_names, fixef_names)]
#
# # Merge and sort by Delta_AICc for ease
# aic = aic[r_squared, on = "Modnames"]
# aic = aic[order(Delta_AICc)]
# aic[, percentageError := Delta_AICc/min(AICc)*100]
# print(paste0("Min AICc: ", aic[, min(AICc)]))
#
# aic_xt = xtable(aic[, .(Modnames, Delta_AICc)])
# r_squared_xt = xtable(setcolorder(aic[, .(Modnames, V1, V2)], c("Modnames", "V1", "V2")))
# aic_rsq_xt = xtable(setcolorder(aic[, .(Modnames, Delta_AICc, V1, V2, percentageError)],
# 	c("Modnames", "Delta_AICc", "V1", "V2", "percentageError")))
#
# print(aic_xt, file = paste0(loadPath, "aic_mega.tex"), include.rownames = FALSE, booktabs = TRUE)
# print(r_squared_xt, file = paste0(loadPath, "rsq_mega.tex"), include.rownames = FALSE, booktabs = TRUE)
# print(aic_rsq_xt, file = paste0(loadPath, "aic_rsq_mega.tex"), include.rownames = FALSE, booktabs = TRUE)
