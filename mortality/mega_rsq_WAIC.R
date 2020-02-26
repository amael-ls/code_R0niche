
#### Aim of prog: Merge R² and WAIC calculated
## Comments:
# There are 7 models fitted in rstanarm.R and 7 simpler alternatives fitted
# in rstanarm_submodels.R (=> 14 comparable models).
#
## Merge WAIC & R²:
#		- Load all the results from waic_mortality(-submodels).R
#		- Bind and save

#### Load package and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(rstanarm)
library(stringi)
library(xtable)

#### Common variables
nbSpecies = 14
nbModels = 7 # Got lucky, I have the same number of models for both rstanarm.R and rstanarm_submodels.R

#### Bind (sub)models for each sp
for (i in 1:nbSpecies)
{
	(loadPath = paste0("./array_", i, "/"))
	ls_waic = vector("list", length = nbModels)
	ls_waic_submodels = vector("list", length = nbModels)
	ls_rsq = vector("list", length = nbModels)
	for (j in 1:nbModels)
	{
		waic = readRDS(paste0(loadPath, "waic", j, ".rds"))
		waic_submodels = readRDS(paste0(loadPath, "waic", j, "_submodels.rds"))
		rsq = readRDS(paste0(loadPath, "rsq_", j, ".rds"))

		ls_waic[[j]] = data.table(model = paste0("model", j),
			waic = waic$estimates["waic", "Estimate"], se = waic$estimates["waic", "SE"])
		ls_waic_submodels[[j]] = data.table(model = paste0("submodel", j),
			waic = waic_submodels$estimates["waic", "Estimate"], se = waic_submodels$estimates["waic", "SE"])
		ls_rsq[[j]] = c(mean(rsq), sd(rsq))
	}
	ls_waic = rbindlist(ls_waic)
	ls_waic_submodels = rbindlist(ls_waic_submodels)
	ls_merge = rbindlist(list(ls_waic, ls_waic_submodels))
	saveRDS(ls_merge, paste0(loadPath, "waic_mega.rds"))
}
