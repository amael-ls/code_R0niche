
#### Aim of prog: Plot three species-specific correations related to R0:
## R0 -- distrib, s* = 0m
#	- Load correlation data
#	- Plot with triangles for each species
#
## R0 -- distrib, s* = 10m
#	- Plot with squares for each species
#
## R0 -- distance to closest edge.
#	- Load data from ../R0 folder (assume the programs corrR0_distEdge.R and correlTable.R ran)
#

#### Load package and clear memory
library(data.table)
library(tikzDevice)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data
## Correlation R0 -- competition
corrR0_comp = readRDS("./correlation.rds")[, .(sp_code, correl_0m, correl_10m)]

## Correlation nearest edge
corrR0_distEdge = readRDS("../R0/correlation.rds")
setnames(corrR0_distEdge, old = "species", new = "sp_code")

## Merge
corrR0_comp = corrR0_comp[corrR0_distEdge, on = "sp_code"]

## Common climaticVariables
nbSpecies = corrR0_comp[, .N]

#### Plot
tikz("3correlations.tex", width = 6.5, height = 5)
op = par(mar = c(0, 3, 0, 2), mgp = c(1.5, 0.75, 0),
	oma = c(0,0,0.9,0), tck = -0.015, las = 1)

plot(x = NULL, y = NULL, xlim = c(0.65, nbSpecies + 0.35),
	ylim = c(-1.3, 1.3), axes = FALSE, xlab = "",
	ylab = "Correlation")

abline(h = 0)

for (i in 1:nbSpecies)
{
	# Species line and name
	if (i %% 4 == 0)
	{
		segments(x0 = i, y0 = -1.15, x1 = i, y1 = 1, lwd = 0.25, lty = "dashed")
		tikzCoord(i, -1.25, 'species_coords')
	}

	if (i %% 4 == 1)
	{
		segments(x0 = i, y0 = -1, x1 = i, y1 = 1, lwd = 0.25, lty = "dashed")
		tikzCoord(i, 1.1, 'species_coords')
	}

	if (i %% 4 == 2)
	{
		segments(x0 = i, y0 = -1, x1 = i, y1 = 1, lwd = 0.25, lty = "dashed")
		tikzCoord(i, -1.1, 'species_coords')
	}

	if (i %% 4 == 3)
	{
		segments(x0 = i, y0 = -1, x1 = i, y1 = 1.15, lwd = 0.25, lty = "dashed")
		tikzCoord(i, 1.25, 'species_coords')
	}

	tikzAnnotate(paste0("\\node (spPos) at (species_coords) {", corrR0_comp[i, sp_code], "};"))

	# Add points
	c0 = corrR0_comp[i, correl_0m]
	points(x = i - 0.15, y = c0, pch = 17, # Triangle
		col = ifelse(i %% 2 == 0, "#ff9933", "#3333ff"))

	c_dist = corrR0_comp[i, correl]
	points(x = i, y = c_dist, pch = 15, # Square
		col = ifelse(i %% 2 == 0, "#ff9933", "#3333ff"))

	correl_10 = corrR0_comp[i, correl_10m]
	points(x = i + 0.15, y = correl_10, pch = 16, # Circle
		col = ifelse(i %% 2 == 0, "#ff9933", "#3333ff"))
}

axis(side = 2, at = seq(-1, 1, by = 0.25))

dev.off()
