
#### Aim of prog: Create the species name table
#

#### Load packages and clear memory
library(data.table)
library(stringi)
library(xtable)

rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool function
## Sorting species (to write their parameters with latex)
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

#### Create the species table
## Common variables
ls_species = readRDS("ls_species.rds")

## Read species code
species_dt = sortingSpecies(ls_species)
species_dt[, c("latin", "vernacular") := .(
	c("Abies balsamea",
		"Acer rubrum", "Acer saccharum",
		"Betula alleghaniensis", "Betula papyrifera",
		"Fagus grandifolia",
		"Picea glauca", "picea mariana", "picea rubens",
		"Pinus banksiana", "Pinus strobus",
		"Populus tremuloides",
		"Thuja occidentalis",
		"Tsuga canadensis"),
	c("Balsam fir",
		"Red maple", "Sugar maple",
		"Yellow birch", "White birch",
		"American beech",
		"White spruce", "Black spruce", "Red spruce",
		"Jack pine", "Eastern white pine",
		"Quaking aspen",
		"Northern white cedar",
		"Eastern hemlock"))]

species_dt[, species_id := paste(tsn, species, sep = "-")]

#### Save the results
## File .rds, for other programs
saveRDS(species_dt, "speciesTables.rds")

## File .tex for article
species_dt_xt = xtable(species_dt[, .(species, latin, vernacular)])

print(species_dt_xt, file = "./speciesTable.tex", include.rownames = FALSE, booktabs = TRUE)
