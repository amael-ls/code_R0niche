
#### Aim of prog: Plot two species-specific correations related to R0:
## R0 -- distrib, s* = 0m
#	- Load correlation data
#	- Plot with triangles for each species
#
## R0 -- distrib, s* = 10m
#	- Plot with squares for each species
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
## Correlation R0 -- competition, random forest
corrR0_comp = readRDS("./correlation.rds") #[, .(sp_code, correl_0m, correl_10m)]

## Common climaticVariables
nbSpecies = corrR0_comp[, .N]

## Load tsn for shade tolerance
tsn = readRDS("../growth/tsn.rds")[, .(species, tolLevel)]
setnames(tsn, old = "species", new = "sp_code")
corrR0_comp = corrR0_comp[tsn, on = "sp_code"]

## Sort by shade tolerance and by species name within shade tolerance, prepare colour
setorderv(x = corrR0_comp, "correl_10m") # cols = c("tolLevel", "sp_code"))
corrR0_comp[tolLevel == "L", colour := "#00F9FF"]
corrR0_comp[tolLevel == "M", colour := "#0E51FF"]
corrR0_comp[tolLevel == "H", colour := "#010120"]

shadeL = col2rgb("#00F9FF")[,1]
shadeM = col2rgb("#0E51FF")[,1]
shadeH = col2rgb("#010120")[,1]

#### Plot
tikz("3correlations.tex", width = 6.5, height = 5)
op = par(mar = c(3, 3, 3, 2), mgp = c(1.5, 0.75, 0),
	oma = c(0,0,2,0), tck = -0.015, las = 1)

plot(x = NULL, y = NULL, xlim = c(0.65, nbSpecies + 0.35),
	ylim = c(-1.3, 1.3), axes = FALSE, xlab = "",
	ylab = "Correlation")

abline(h = 0)

for (i in 1:nbSpecies)
{
	# Species line and name
	if (i %% 4 == 0)
	{
		tikzCoord(i, -1.15, paste0("pt_start_", i))
		tikzCoord(i, 1, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(i, -1.25, 'species_coords')
	}

	if (i %% 4 == 1)
	{
		tikzCoord(i, -1, paste0("pt_start_", i))
		tikzCoord(i, 1, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(i, 1.1, 'species_coords')
	}

	if (i %% 4 == 2)
	{
		tikzCoord(i, -1, paste0("pt_start_", i))
		tikzCoord(i, 1, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(i, -1.1, 'species_coords')
	}

	if (i %% 4 == 3)
	{
		tikzCoord(i, -1, paste0("pt_start_", i))
		tikzCoord(i, 1.15, paste0("pt_end_", i))
		tikzAnnotate(paste0("\\draw[loosely dashed] (pt_start_", i, ") -- (pt_end_", i, ");"))
		tikzCoord(i, 1.25, 'species_coords')
	}

	tikzAnnotate(paste0("\\node (spPos) at (species_coords) {", corrR0_comp[i, sp_code], "};"))

	# Add points
	c0 = corrR0_comp[i, correl_0m]
	points(x = i - 0.15, y = c0, pch = 17, # Triangle
		col = corrR0_comp[i, colour])

	correl_10 = corrR0_comp[i, correl_10m]
	points(x = i + 0.15, y = correl_10, pch = 16, # Circle
		col = corrR0_comp[i, colour])
}

axis(side = 2, at = seq(-1, 1, by = 0.25),
	labels = c("-1", "", "-0.5", "", "0", "", "0.5", "", "1"))

## Legend
tikzAnnotate(paste0("\\definecolor{shadeL}{RGB}{", paste(shadeL, collapse = ","), "}"))
tikzAnnotate(paste0("\\definecolor{shadeM}{RGB}{", paste(shadeM, collapse = ","), "}"))
tikzAnnotate(paste0("\\definecolor{shadeH}{RGB}{", paste(shadeH, collapse = ","), "}"))

tikzAnnotate("
	\\matrix [below right = 1cm and 0cm] at (current bounding box.north west) {
		\\node [shape = rectangle, fill = shadeL, label = right:Low] {}; &
		\\node [shape = rectangle, fill = shadeM, label = right:Medium] {}; &
		\\node [shape = rectangle, fill = shadeH, label = right:High] {}; \\\\
	};
")

tikzAnnotate("\\node [below right = 0.5cm and 0cm] at (current bounding box.north west) {Shade tolerance:};")

dev.off()
