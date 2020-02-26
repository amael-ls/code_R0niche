
#### Aim of prog: create table for article
#	- Load correlation R0 -- distance
#	- Save table in a texfile

#### Load package and clear memory
library(data.table)
library(stringi)
library(xtable)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool function
sortingSpecies = function(ls_species)
{
	tsn_dt = data.table(species = character(length(ls_species)),
		tsn = character(length(ls_species))) # tsn = Taxonomic Serial Number

	tsn_dt[, species := stri_sub(str = ls_species,
		from = stri_locate_first(str = ls_species, regex = "-")[,1] + 1)]

	tsn_dt[, tsn := stri_sub(str = ls_species, from = 1,
		to = stri_locate_first(str = ls_species, regex = "-")[,1] - 1)]
	tsn_dt = tsn_dt[order(species),]
	return (tsn_dt)
}

#### Common variables
## List species
ls_species = dir(path = "./results", pattern = "[0-9]{4,}")
n = length(ls_species)

## Data table results
results = sortingSpecies(ls_species) # data.table(species = ls_species, correl = numeric(n))
results[, correl := numeric(n)]

#### Data table resutls
## Load data
for (sp in ls_species)
	results[stri_detect(str = sp, regex = tsn),
		correl := readRDS(paste0("./results/", sp, "/correl_R0_distEdge.rds"))]

saveRDS(results, "./correlation.rds")

## Save tex file
tsn_xt = xtable(results[, !"tsn"])
print(tsn_xt, file = "./correlation.tex", include.rownames = FALSE, booktabs = TRUE)
