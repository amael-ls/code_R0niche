
#### Aims of prog
## Merge climate and tree data
# I assume the climate database exists, and that competition is already calculated
#
## Create growth and mortality data set
#	- Growth
#		Keep species that have enough measurement (discard single-measured trees)
#		Delete dead trees (cannot grow)
#		Keep all dbh, although crown/height allometries might not be right below 90
#
#	- mortality
#		Keep species that have enough measurement (discard single-measured trees)
#
## Schubert impromptu Maria Joao-Pirez https://www.youtube.com/watch?v=5yVZu05WZ9o

#### Load packages and tool functions; clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

library(data.table)
library(magrittr)
# library(dplyr)

# source("../toolFunctions.R")
source("./parametersAllometries.R")

#### Load data
## Tree data
treeData = readRDS("./tree_sStar.rds")

## Climate data
clim = readRDS("./averagedClim_5yearsAllVar.rds")

#### Merge data table tree data and 5-year-averaged climate
## Loss of data due to merging (clim had NA values for some years)
treeData = treeData[clim, on = c("latitude", "longitude", "year_measured"), nomatch = 0]
saveRDS(treeData, "treeClim_sStar.rds")

#### Subset valid for both growth and mortality
## Keep only 14 species
ls_14species = c("19481-BET-ALL", "19489-BET-PAP", "28731-ACE-SAC", "28728-ACE-RUB",
	"19462-FAG-GRA", "183295-PIC-GLA", "183397-TSU-CAN", "195773-POP-TRE", "183385-PIN-STR",
	"18032-ABI-BAL", "505490-THU-OCC", "183302-PIC-MAR", "18034-PIC-RUB", "183319-PIN-BAN")

treeData = treeData[species_id %in% ls_14species]

## Keep only non single-measured trees
treeData[, nbMeasure := .N, by = tree_id]
treeData = treeData[nbMeasure > 1]

sort(treeData[, table(species_id)])

#### -------------------- ####
##			GROWTH			##
#### -------------------- ####
#### Create growth
## Keep only living trees
growth_dt = treeData[is_dead == "f"]

## Trees that have been measured twice, alive and then dead cannot be in growth_dt
growth_dt[, nbMeasure := .N, by = tree_id]
growth_dt = growth_dt[nbMeasure > 1]

## Shift the dbhIncr, growth and deltaYear of one rank (to keep event at t0 and not t1)
# /!\ growth_dt must be ordered by increasing years /!\
shiftData = function(vec)
	return (c(vec[2:length(vec)], NA))

growth_dt = setorderv(growth_dt, "year_measured", 1)
growth_dt[, c("dbhIncr", "deltaYear", "growth") := lapply(.SD, shiftData), by = tree_id,
	.SDcols = c("dbhIncr", "deltaYear", "growth")]

## Remove the NA of dbhIncr (i.e., the last measure of the remaining trees)
growth_dt = na.omit(growth_dt)

## Calculate canopy status TRUE/FALSE
growth_dt[, canopy_status := height > s_star]

## Keep only species that have more than 1e4 individuals (should be the case with the 14 sp)
table_species = sort(unique(growth_dt[, .(species_id, tree_id)]) %>%
    .[, table(species_id)])
lim = 1e4

ls_species = sort(names(table_species[table_species > lim]))
saveRDS(ls_species, "ls_species.rds")

growth_dt = growth_dt[species_id %in% ls_species]

#K Keep positive growth (~ 5% of them are zeros)
growth_dt = growth_dt[growth > 0]

## Save growth_dt
saveRDS(growth_dt, "growth_dt.rds")

#### ------------------------ ####
##			MORTALITY			##
#### ------------------------ ####
#### Mortality, the subset for non-single-measured tree is already done
mortality_dt = treeData[species_id %in% ls_14species]

## Rewrite mortality as a transition
# Sort mortality by increasing year
setorderv(mortality_dt, cols = "year_measured", order = +1)

# Change is_dead by TRUE/FALSE rather than characters "t", "f"
mortality_dt[, is_dead := ifelse(is_dead == "t", TRUE, FALSE)]

# Keep only the first death record, i.e. F[...]FT[...]T becomes F[...]FT
#		/!\ database sorted by increasing years /!\
trackFirstDeath = function(vec_lifeState)
{
    ind_alive = which(vec_lifeState == FALSE)
    ind_dead = which(vec_lifeState == TRUE)

    toKeep = rep(TRUE, length(vec_lifeState))

    if (length(ind_dead) > 1)
        toKeep[ind_dead[2:length(ind_dead)]] = FALSE

    return(toKeep)
}

mortality_dt[, keepTree := trackFirstDeath(is_dead), by = tree_id]

mortality_dt = mortality_dt[keepTree == TRUE]
mortality_dt[, keepTree := NULL]

# Remove single-measured trees (should be only trees that have always been recorded dead)
mortality_dt[, nbMeasure := .N, by = tree_id]
mortality_dt = mortality_dt[nbMeasure > 1]
mortality_dt[, nbMeasure := NULL]

## Calculate mortality event
#		/!\ treeData must be ordered by decreasing years /!\
delta = function(x)
	return (c(x[1:(length(x) - 1)] - x[2:length(x)], NA)) # NA for the last row which will be discarded

# Sort by decreasing year order
setorderv(mortality_dt, cols = "year_measured", order = -1)
mortality_dt[, deltaState := delta(is_dead), by = tree_id]

## Shift the dbhIncr, growth, deltaYear and deltaState of one rank (to keep event at t0 and not t1)
# /!\ mortality_dt must be ordered by decreasing years (already the case) /!\
shiftData = function(vec)
	return (c(NA, vec[1:(length(vec) - 1)]))

mortality_dt[, c("dbhIncr", "deltaYear", "growth", "deltaState") := lapply(.SD, shiftData),
	by = tree_id, .SDcols = c("dbhIncr", "deltaYear", "growth", "deltaState")]

## Canopy status
mortality_dt[, canopy_status := height > s_star]
mortality_dt[, .N, by = canopy_status]
mortality_dt[, .N, by = c("canopy_status", "deltaState")]

## Delete measure with a too big/small time interval, as mortality event is dependent of freq. obs.
mortality_dt[, .N, by = c("deltaState", "deltaYear")]
mortality_dt = mortality_dt[is.na(deltaYear) | ((4 < deltaYear) & (deltaYear < 11))]
mortality_dt[, .N, by = c("deltaState", "deltaYear")]

## Shall I delete small trees due to biased sample, threshold = 91 mm?
mortality_dt[, table(dbh)][1:120]
mortality_dt[dbh < 91, table(deltaState)]
# mortality_dt = mortality_dt[dbh > 90] # I decided to keep all trees

# Remove single-measured trees (trees that had two records with one deltaYear > 10)
mortality_dt[, nbMeasure := .N, by = tree_id]
mortality_dt = mortality_dt[nbMeasure > 1]
mortality_dt[, nbMeasure := NULL]

## Keep states at time 0 (i.e., year_measured) and t1, keep time t1
# Obviously, is_dead_t1 is useless as it is a duplicate of deltaState
mortality_dt[, is_dead_yr_measured := is_dead]
mortality_dt[, is_dead_t1 := deltaState]
mortality_dt[, t1 := year_measured + deltaYear]

mortality_dt = mortality_dt[!is.na(deltaState)]
mortality_dt[, deltaState := as.logical(deltaState)]
mortality_dt[, is_dead_t1 := as.logical(is_dead_t1)]

## Delete is_dead column (duplicate of is_dead_yr_measured)
mortality_dt[, is_dead := NULL]

mortality_dt[, .N, by = c("canopy_status", "deltaState")]

saveRDS(mortality_dt, "mortality_dt.rds")
