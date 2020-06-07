
#### Aim of prog: create table for article
#	- Load correlation R0 -- distance
#	- Save table in a texfile
#
## Remark:
# correl_*_proj (where * is either north or south) are correlations computed
# over the points closer to a border belonging to the northern region, or
# closer to a border belonging to the southern region

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
results[, c("correl_tot_10", "correl_north_10", "correl_south_10", "correl_north_proj_10", "correl_south_proj_10",
	"correl_tot_0", "correl_north_0", "correl_south_0", "correl_north_proj_0", "correl_south_proj_0") := numeric(n)]

#### Data table resutls
## Load data
for (sp in ls_species)
{
	results[stri_detect(str = sp, regex = tsn),
		correl_tot_10 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_north_10 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_north.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_north_proj_10 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_north_proj.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_south_10 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_south.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_south_proj_10 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_south_proj.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_tot_0 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_0m.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_north_0 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_north_0m.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_north_proj_0 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_north_proj_0m.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_south_0 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_south_0m.rds"))]

	results[stri_detect(str = sp, regex = tsn),
		correl_south_proj_0 := readRDS(paste0("./results/", sp, "/correl_R0_distEdge_south_proj_0m.rds"))]
}

saveRDS(results, "./correlation.rds")

## Save tex file, only *proj correls (cf Remark) and tot correls
# Delete columns
results[, c("tsn", "correl_north_10", "correl_south_10", "correl_north_0", "correl_south_0") := NULL]

# Set col order
setcolorder(x = results, neworder = c("species", "correl_tot_0", "correl_tot_10",
	"correl_north_proj_0", "correl_north_proj_10",
	"correl_south_proj_0", "correl_south_proj_10"))

# Print xtable .tex file
tsn_xt = xtable(results)
print(tsn_xt, file = "./correlation.tex", include.rownames = FALSE, booktabs = TRUE)
