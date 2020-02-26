
#### Aims of prog: Clean the tree data that will be later used to parametrise demography
## 	Clean data
#		- Keep species that have dbh--crown allometries (109 - 3 = 106 species)
#		- Insane dbh
#		- Unknown leaving states
#		- No climate data available (they stop in 2010)
#		- Height above 85 m of living trees
#		- Same individual has different species id
#		- Resurrected trees (the last stage is considered the true one)
#		- Insane diametral growth (considered insane after 25mm/yr for our database) of living trees
#
# /!\ It is important to have the data ordered by decreasing years for resurrected checks /!\
#
# Removing unknown crown--dbh allometries:
# 183336-PIN-EDU 194812-JUN-ASH 194859-JUN-OST
#		   31759          10935          46065

#### Load packages and tool functions; clear memory and graphs
library(data.table)
library(dplyr)

rm(list = ls())
graphics.off()
options(max.print = 500)

source("../toolFunctions.R")
source("./parametersAllometries.R")

#### Load data
treeData = readRDS("~/database/trees/north-america_QUICCFOR/treeData.rds")
setnames(x = treeData, old = "id_spe", new = "species_id")

parametrisedSpecies = readRDS("./ls_speciesParametrised.rds")$species_id # i.e., have allometries
parametrisedSpecies = parametrisedSpecies[-which(parametrisedSpecies %in%
	c("183336-PIN-EDU", "194812-JUN-ASH", "194859-JUN-OST"))]

class(treeData)
names(treeData)
range_df(treeData)

#### Discard insane dbh, unknown leaving state, post 2010 mesures, and no dbh--crown allometries
treeData = treeData[dbh < 9999, ]
treeData = treeData[is_dead != "",]
treeData = treeData[year_measured < 2011, ]
treeData = treeData[species_id %in% parametrisedSpecies, ]

#### Remove height > 85
## Calculate height
treeData[, height := dbhToHeight(dbh, purves2007_allometries[species == species_id, a],
		purves2007_allometries[species == species_id, b], mm = TRUE), by = species_id]

## Remove h > 85
# Why 85? Most of the trees above this height are Pseudotsuga menziesii, according
# to the silvics (Russel M. Burns, 1990), it is one of the biggest tree of my database.
# Here is a description:
# Trees 150 to 180 cm (60 to 72 in) in diameter and 76 m (250 ft) in height are
# common in old-growth forests (22). The tallest tree on record, found near Little Rock, WA,
# was 100.5 m (330 ft) tall and had a diameter of 182 cm (71.6 in).
#
# treeData[height > 85 & is_dead == "f", table(species_id)]
#
#  18032-ABI-BAL  18034-PIC-RUB  18044-THU-PLI 181824-ABI-AMA 183327-PIN-CON
#              1              1              1              2              4
# 183400-TSU-HET 183417-LAR-OCC 183424-PSE-MEN
#              1              3             20

treeData = treeData[(is_dead == "t") | ((is_dead == "f") & (height <= 85))]

#### Check individual that have multi species id
checkMultiSp = function(vec_speciesId)
{
	l = length(unique(vec_speciesId))
	if (l > 1)
		return(TRUE) # i.e., there is a species conflict

	return(FALSE) # i.e., there is no species conflict
}

treeData[, toNotKeep := checkMultiSp(species_id), by = tree_id]
treeData = treeData[!(toNotKeep), ]
treeData = treeData[, toNotKeep := NULL]

#### Check resurrected, /!\ treeData must be ordered by decreasing years /!\
## Function to modify life state
modifLifeState = function(vec_lifeState) # Last state is the truth
{
	ind_alive = which(vec_lifeState == "f")
	ind_dead = which(vec_lifeState == "t")
	if (vec_lifeState[1] == "t")
		vec_lifeState[ind_alive[which(ind_alive < max(ind_dead))]] = "t"

	if (vec_lifeState[1] == "f")
		vec_lifeState[ind_dead[which(min(ind_alive) < ind_dead)]] = "f"

	if (length(vec_lifeState) == 2)
		if (isTRUE(all.equal(vec_lifeState, c("f", "t"))))
			vec_lifeState[2] = "f"

	return(vec_lifeState)
}

## Correcting resurrected trees, the last state is considered as the truth
# /!\ treeData must be ordered by decreasing years /!\
setorderv(treeData, cols = "year_measured", order = -1)
ls_deadIndiv = treeData[is_dead == "t", unique(tree_id)]
treeData[tree_id %in% ls_deadIndiv, is_dead := modifLifeState(is_dead), by = .(tree_id)]

# Set to 0 the height of dead trees (their trunk can be big due to decomposition)
treeData[is_dead == "t", height := 0] # Has to be done after correction of dead status

#### Delete insane growth
## Find individuals that have at least two measures, but do not delete single-measured trees!
treeData[, nbMeasure := .N, by = tree_id]

## Calculate growth, /!\ treeData must be ordered by decreasing years /!\
delta = function(x)
	return (c(x[1:(length(x) - 1)] - x[2:length(x)], NA)) # NA for the last row which will be discarded

setorderv(treeData, cols = "year_measured", order = -1)
treeData[nbMeasure > 1, c("dbhIncr", "deltaYear") := lapply(.SD, delta), by = .(tree_id),
	.SDcols = c("dbh", "year_measured")]

treeData[, growth := dbhIncr/deltaYear]

## Delete growth above the (arbitrary) threshold of 25mm
# obviously, I will remove the dead for the growth data; here, the unrealistic growth are
# deleted because it is a proof of dbh measurement mistake for living trees. Note that I keep negative
# growth of dead trees as they might shrink or expand due to breaking/decaying. It is not a problem since it is
# the dbh at time t - 1 that will be the predictor of mortality, not dbh at time t
# I also keep NA measure for now, since they are either single-measured trees or the first measure

treeData = treeData[is.na(growth) | is_dead == "t" | ((is_dead == "f") & (growth >=0) & (growth <= 25)),]

## Verif that insane growth (i.e., < 0 or > 25 mm/yr) contains only dead trees
range(treeData[growth < 0 | growth > 25, is_dead])

#### Check plots that have multi area within a year
checkMultiArea = function(vec_plotArea)
{
	l = length(unique(vec_plotArea))
	if (l > 1)
		return(TRUE) # i.e., there is a plot area conflict

	return(FALSE) # i.e., there is no plot area conflict
}

treeData[, toNotKeep := checkMultiArea(plot_size), by = .(plot_id, year_measured)]
treeData = treeData[!(toNotKeep), ]
treeData = treeData[, toNotKeep := NULL]

saveRDS(treeData, "treeData_cleaned.rds")
