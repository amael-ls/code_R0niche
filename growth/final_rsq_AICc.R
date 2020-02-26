
#### Aim of prog: Calculate R² and AICc chosen model (the 2) and fixef models
## Comments:
# The climatic selection was done ReML=FALSE
#
## Calculate AICc & R²:
#		- List and load all the concerned models
#		- Calculate AICc using AICcmodavg package
#		- Calculate R² using MuMIn
#		- Save the results in rds and tex format

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
## List files fixef
ls_files_fixef = list.files(path = loadPath, pattern = "fixef.rds")
ls_files_fixef = Filter(function(x) grepl("model", x), ls_files_fixef)
(n_fixef = length(ls_files_fixef))

n = n_fixef + 1 # (+1 for model 2 from ReML=FALSE.R prog)

model_ls = vector("list", length = n)

## Create models' name (correspond to glmm_fixef.R names)
fixef_names = stri_sub(ls_files_fixef,
	to = stri_locate_first(ls_files_fixef, regex = "_")[,1] - 1)

## Read files
# Fixef
for (i in 1:n_fixef)
	model_ls[[i]] = readRDS(paste0(loadPath, ls_files_fixef[i]))

# Model 2 from ReML=FALSE.R, /!\ cannot load final model as REML = TRUE in the estimation /!|
model_ls[[n]] = readRDS(paste0(loadPath, "model2_ReML=FALSE.rds"))

#### AICc and R²
aic = AICcmodavg::aictab(model_ls, modnames = c(fixef_names, "chosenModel"))
r_squared = lapply(model_ls, MuMIn::r.squaredGLMM)

## Reshape r_squared, to make it easier to read
rsq = as.data.table(matrix(data = unlist(r_squared), nrow = n, ncol = 2, byrow = TRUE))
setnames(rsq, old = c("V1", "V2"), new = c("R2_m", "R2_c"))
modnames_fixef = stri_sub(str = ls_files_fixef,
	to = stri_locate_first(str = ls_files_fixef, regex = "_")[, 1] - 1)
rsq[, Modnames := c(modnames_fixef, "chosenModel")]

#### Save results
## rds format
saveRDS(aic, paste0(loadPath, "aic_final.rds"))
saveRDS(rsq, paste0(loadPath, "rsq_final.rds"))

## latex format
# Transform to data table
setDT(aic)

# Merge and sort by Delta_AICc for ease
aic = aic[rsq, on = "Modnames"]
aic = aic[order(Delta_AICc)]
aic[, percentageError := Delta_AICc/min(AICc)*100]
print(paste0("Min AICc: ", aic[, min(AICc)]))

aic_xt = xtable(aic[, .(Modnames, Delta_AICc)])
rsq_xt = xtable(setcolorder(aic[, .(Modnames, R2_m, R2_c)], c("Modnames", "R2_m", "R2_c")))
aic_rsq_xt = xtable(setcolorder(aic[, .(Modnames, Delta_AICc, R2_m, R2_c, percentageError)],
	c("Modnames", "Delta_AICc", "R2_m", "R2_c", "percentageError")))

print(aic_xt, file = paste0(loadPath, "aic_final.tex"), include.rownames = FALSE, booktabs = TRUE)
print(rsq_xt, file = paste0(loadPath, "rsq_final.tex"), include.rownames = FALSE, booktabs = TRUE)
print(aic_rsq_xt, file = paste0(loadPath, "aic_rsq_final.tex"), include.rownames = FALSE, booktabs = TRUE)
