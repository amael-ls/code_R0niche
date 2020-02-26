
#### Aim of prog: Slice the data into N subset to fasten the competition computation
#

#### Load packages and tool functions; clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

library(data.table)

#### Load data
## Tree data
treeData = readRDS("./treeData_cleaned.rds")

#### Subset the data
## List plot/year combinations
treeData[, plot_year := paste0(plot_id, "_", year_measured)]
plot_year_combination = unique(treeData[, plot_year])

## Number of arrays you are going to use for competition
nb_arrays = 5000

## Calculate rows to deal with
deltaStep = (length(plot_year_combination) - 1)/nb_arrays
rows = ceiling(seq(1, length(plot_year_combination), deltaStep))
rows[length(rows)] = rows[length(rows)] + 1 # (for the -1 below)

for (array_id in 1:nb_arrays)
{
	slice = rows[array_id:(array_id + 1)]
	# Remove the last value that is the starting of the next array_id
	slice[2] = slice[2] - 1 # Works even for the last array_id thanks to my +1 above

	# Subset
	selectedCombinations = plot_year_combination[slice[1]:slice[2]]
	slicedData = treeData[plot_year %in% selectedCombinations]
	saveRDS(slicedData, paste0("./slice/", array_id, ".rds"))
	print(paste0("array ", array_id, " done"))
}
