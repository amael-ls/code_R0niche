
#### Aims of prog: Calculate the waic (Gelman2013, p. 181-183)
## 7 models are compared

#### Load packages
library(doParallel)
library(rstanarm)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Common variables
## Number of formula (which is lukily the same for models and submodels, otherwise it should have been an argument to provide)
nbFormulae = 7

## Get the extra argument to distinguish model/submodel
bash_args = commandArgs(trailingOnly = TRUE)

if (length(bash_args) == 0)
{
	print(paste0("No argument provided. By default the program compute waic of models"))
	bash_args = ""
}

if (bash_args == "submodel")
	print(paste0("The argument is ", bash_args, ". The program will run on the submodels"))

if ((bash_args != "") & (bash_args != "submodel"))
  stop("Error, only no argument or `submodel' argument allowed", call.=FALSE)

print(paste0("bash_args = <", bash_args, ">"))

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

#### Load data
## Model
(loadPath = paste0("./array_", species_id, "/"))
(ls_files = list.files(path = loadPath, pattern = ifelse(bash_args == "submodel", "^submodel_", "^model_")))

model = readRDS(paste0(loadPath, ls_files[formula_id]))

#### Run waic
start_time = Sys.time()
waic_model = waic(x = model)
waic_name = paste0(loadPath, "waic", formula_id, ifelse(bash_args == "", "", "_"), bash_args, ".rds")
saveRDS(waic_model, waic_name)
end_time = Sys.time()

print(end_time - start_time)
