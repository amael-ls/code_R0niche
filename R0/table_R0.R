
#### Aim of prog: create a table of R0 in .tex format

#### Load libraries
library(data.table)
library(xtable)

options(max.print = 500)
rm(list = ls())

#### Common variables
## List folders
ls_species = list.files(path = "./results/", pattern = "^[0-9]{4,}")

## Data table R0 min max
n = length(ls_species)
R0_min_max = data.table(species_id = ls_species,
	maxR0_10 = numeric(n), minR0_10 = numeric(n), prop_10 = numeric(n),
	maxR0_0 = numeric(n), minR0_0 = numeric(n), prop_0 = numeric(n))

## Load tsn
tsn = readRDS("../growth/tsn.rds")[, .(species_id, species)]

## Merge tsn and R0_min_max
R0_min_max = R0_min_max[tsn, on = "species_id"]
setorderv(R0_min_max, "species", +1)

#### Run
## Compute data
for (sp in ls_species)
{
	# Load R0 with competition (s* = 10 m) and without
	R0_10 = fread(paste0("./results/", sp, "/R0_10m.csv"))
	R0_0 = fread(paste0("./results/", sp, "/R0_0m.csv"))

	R0_min_max[species_id == sp, maxR0_10 := max(R0_10)]
	R0_min_max[species_id == sp, minR0_10 := min(R0_10)]
	R0_min_max[species_id == sp, prop_10 := R0_10[V1 < 1, .N]*100/R0_10[, .N]]

	R0_min_max[species_id == sp, maxR0_0 := max(R0_0)]
	R0_min_max[species_id == sp, minR0_0 := min(R0_0)]
	R0_min_max[species_id == sp, prop_0 := R0_0[V1 < 1, .N]*100/R0_0[, .N]]

	print(paste0("species: ", sp, " done"))
}

## Save results
saveRDS(R0_min_max, "./R0_min_max.rds")

R0_min_max_xt = xtable(R0_min_max[, c(8, 2:7)], digits = 2)
print(R0_min_max_xt, file = "./R0_min_max.tex", include.rownames = FALSE, booktabs = TRUE)
