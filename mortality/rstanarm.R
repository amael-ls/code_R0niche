
#### Aims of prog: Parametrise mortality using rstanarm (bayesian method)
## 7 models are compared:
#	- random structure: plot_id
#		* k, pp
#		* k, pp (quadratic)
#		* cs
#		* dbh
#		* dbh (quadratic)
#		* k, pp, dbh, cs (quadratic)
#		* cs:(k, pp), dbh (quadratic)
#		* cs:(k, pp), (k, pp):dbh (quadratic)
#

#### Load packages
library(data.table)
library(doParallel)
library(rstanarm)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

source("../toolFunctions.R")

#### Common variables
## Number of formulae
nbFormulae = 7

## Get the extra argument to distinguish model/submodel
bash_args = commandArgs(trailingOnly = TRUE)

if (length(bash_args) == 0)
{
	print(paste0("No argument provided. By default the program compute the models to select climatic variables"))
	bash_args = ""
}

if (bash_args == "submodel")
	print(paste0("The argument is ", bash_args, ". The program will compute the submodels"))

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

#### Make directory and saving path
(savePath = paste0("./array_", species_id, "/"))

if (bash_args == "")
{
	if (!dir.exists(savePath))
		dir.create(savePath)
} else {
	if (!dir.exists(savePath))
		stop(paste0("*** Error: folder ", savePath, " does not exist"), call.=FALSE)
}
#### Load data and transform them
# mortality_data = readRDS("../createData/mortality_dt.rds")
mortality_data = readRDS("../createData/mortality_dt.rds")
climaticVariables = readRDS("../createData/climaticVariables.rds")

## Subset for the species of interest and keep deltaYear \in [5,11]
ls_species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

(species = ls_species[species_id])

mortality_data = mortality_data[(4 < deltaYear) & (deltaYear < 12)]

## Write species name in the folder in case of
if (!file.exists(paste0(savePath, species, ".txt")))
	write(x = species, file = paste0(savePath, species, ".txt"))

mortality_data = mortality_data[species_id == species]
print(paste0("number of levels plot_id: ", length(unique(mortality_data[, plot_id]))))
print(paste0("number of data: ", mortality_data[, .N]))

print("Distribution of the deltaYear:")
print(mortality_data[, table(deltaYear)])

## Normalise explanatory variables
mortality_data_norm = normalisation(df = mortality_data,
	colnames = c("dbh", climaticVariables),
	filename = ifelse(bash_args == "", paste0(savePath, "normalisation_mortality_data.rds"), "./rubbish_normalisation.rds"))

range_df(mortality_data_norm)

## Convert plot_id, species_id and canopy_status as a factor
mortality_data_norm[, plot_id := as.factor(plot_id)]
mortality_data_norm[, species_id := as.factor(species_id)]
mortality_data_norm[, canopy_status := as.factor(canopy_status)]
mortality_data_norm[, year_measured := as.factor(year_measured)]

## Get some information
print("Number of tree measurement alive and dead in the under/over-storey:")
mortality_data[, .N, by = c("canopy_status", "deltaState")]

setDF(mortality_data_norm)

#### Formulae, I assume the random structure is the same as growth
## Select the best climatic variables
if (bash_args == "")
	source("./selectClimaticVariables_inter.R")

if (bash_args == "submodel")
	source("./submodels.R")

print(paste0("glmm rstan ", array_id, " done"))
stopCluster(cl)
