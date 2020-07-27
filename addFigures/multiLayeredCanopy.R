
#### Aim of prog: Check how many layers are in the canopy
#

#### Load package and clear memory
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

source("../createData/parametersAllometries.R")

#### Tool functions
## Compute s* for a given patch and year
compute_sStar_patch = function(height, supBound, a, b, T_param, plotArea, id, tolSize = 0.01)
{
	supBound = max(supBound)
	infBound = 0
	sumArea = 0
	s_star = 0

	beyondSupBound_ind = supBound < height # Exclude trees taller than supBound (useful for 2nd, 3rd, ... layer computation)

	if (length(plotArea) != 1)
		print(id)

	while (supBound - infBound > tolSize)
	{
		s_star = (supBound + infBound)/2
		understorey_ind = s_star <= height
		indices = understorey_ind & !beyondSupBound_ind
		
		if (!any(indices))
		{
			supBound = s_star
			next;
		}

		distToTop = height[indices] - s_star
		sumArea = sum(heightToCrownArea(height[indices], distToTop, a[indices], b[indices], T_param[indices], C0_C1))
		
		if (sumArea < plotArea)
			supBound = s_star
		if (sumArea >= plotArea)
			infBound = s_star
	}
	return (s_star)
}

## Check that sumArea(s*) equals plotArea (i.e., not too far) if s* != 0
check_sStar = function(height, supBound, s_star, a, b, T_param, plotArea, id)
{
	if (length(plotArea) != 1)
		print(id)

	if (length(s_star) != 1)
		print(id)

	supBound = max(supBound)
	
	# Index of heights to keep
	beyondSupBound_ind = supBound < height # Exclude trees taller than supBound (useful for 2nd, 3rd, ... layer computation)
	overstorey_ind = s_star <= height
	indices = overstorey_ind & !beyondSupBound_ind

	if (!any(indices))
		return (0);

	distToTop = height[indices] - s_star
	sumArea = sum(heightToCrownArea(height[indices], distToTop, a[indices],
		b[indices], T_param[indices], C0_C1))

	return (sumArea)
}

#### Load data and subset them
## Load tree data
treeData = readRDS("../createData/treeData_cleaned.rds")

## Keep only leaving trees and eastern north american trees
treeData = treeData[(is_dead == "f") & (longitude > -123)]

## Keep only variable of interest
treeData = treeData[, .(plot_id, plot_size, year_measured, species_id, dbh, height)]

## Unique Id
treeData[, id := paste(plot_id, year_measured, sep = "-")]

## Sort by height (no need to do it within id)
setorder(x = treeData, -height)

## Set a, b, and T
# a
treeData[, a := purves2007_allometries[species == species_id, a], by = species_id]

# b
treeData[, b := purves2007_allometries[species == species_id, b], by = species_id]

# T
treeData[, T_param := purves2007_allometries[species == species_id, T], by = species_id]

#### Compute s*
treeData[, s_star := compute_sStar_patch(height, max(height), a, b, T_param, unique(plot_size), id), by = id]
treeData[, sumArea1 := check_sStar(height, max(height), unique(s_star), a, b, T_param, unique(plot_size), id), by = id]

treeData[, s_star2 := compute_sStar_patch(height, s_star, a, b, T_param, unique(plot_size), id), by = id]
treeData[, sumArea2 := check_sStar(height, unique(s_star), unique(s_star2), a, b, T_param, unique(plot_size), id), by = id]

treeData[, s_star3 := compute_sStar_patch(height, s_star2, a, b, T_param, unique(plot_size), id), by = id]
treeData[, sumArea3 := check_sStar(height, unique(s_star2), unique(s_star3), a, b, T_param, unique(plot_size), id), by = id]

saveRDS(treeData, "./treeWith_s_star.rds")

#### Crash test zone
# treeData[, diff := plot_size - sumArea]
# aa = treeData[, .(plot_id, plot_size, year_measured, s_star, diff)]
# aa[, relDiff := diff*100/plot_size]
# aa = aa[diff < -5]

# pb = aa[relDiff < -10, unique(plot_id)]

# uu = treeData[plot_id %in% pb]
# dt = setorder(uu[, .N, by = c("plot_id", "year_measured")], plot_id)


# bb = treeData[id == "662742-2008"]

# height = bb$height
# supBound = bb$s_star
# a = bb$a
# b = bb$b
# T_param = bb$T_param
# plotArea = unique(bb$plot_size)
