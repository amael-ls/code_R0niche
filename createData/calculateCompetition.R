
#### Aim of prog: Calculate s* for any plot-year combination
# It is assumed the program slicer.R has alredy been called.

#### Load packages and tool functions; clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

library(data.table)

# source("../toolFunctions.R")
source("./parametersAllometries.R")

#### Access parallel variables
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

# nb_arrays = as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
# print(paste0("Number of arrays = ", nb_arrays))

#### Load data
## Tree data
(file = paste0("./slice/", array_id, ".rds"))
treeData = readRDS(file)

#### Calculate competition
## Function to calculate s*
# dbh and species_id are vectors; plotArea a scalar
# Trees are sorted by decreasing height
canopyHeight = function(height, species_id, plotArea, C0_C1, plot_id, plot_year)
{
	print(unique(plot_year))
	sumArea = 0
	i = 1
	n = length(height)
	plot_size = plotArea[1]

	while ((sumArea < plot_size) & (i <= n))
	{
		distToTop = height[1:i] - height[i]
		paramsAllometries = getAllometries(species_id[1:i])
		sumArea = sum(heightToCrownArea(height[1:i], distToTop, paramsAllometries$a, paramsAllometries$b,
			paramsAllometries$T_param, C0_C1)) # For i = 1, this is always 0, but cannot always start at i = 2
		i = i + 1
	}

	if (sumArea < plot_size) # i.e., gap in the canopy
		return (0)

	return (height[i - 1])
}

# Sort by plot id, year (within plot_id), and then decreasing height (within plot AND years)
setorderv(treeData, cols = c("plot_id", "year_measured", "height"), order = c(1, 1, -1))

treeData[, s_star := canopyHeight(height, species_id, plot_size, C0_C1, plot_id, plot_year),
	by = plot_year]

saveRDS(treeData, paste0("./competition/", array_id, ".rds"))
