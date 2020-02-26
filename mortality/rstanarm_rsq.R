
#### Aims of prog: Calculate the R^2 of the 7 bayesian regression models
## 7 models are compared:
#

#### Load packages
library(doParallel)
library(rstanarm)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Common variables
nbFormulae = 7

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

#### Load model
(loadPath = paste0("./array_", species_id, "/"))

model = readRDS(paste0(loadPath, "model_", formula_id, ".rds"))
rsq = bayes_R2(model)

saveRDS(rsq, paste0(loadPath, "rsq_", formula_id, ".rds"))
