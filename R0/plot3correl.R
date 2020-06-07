
#### Aim of prog: Plot three species-specific correations related to R0:
## R0 -- distance, s* = height_canopy (0 or 10, to set in the program)
#	- Load correlation data
#	- Plot with triangles up/down (north/south), circles (total) for each species

#### Load package and clear memory
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data
## Correlation nearest edge
corrR0_distEdge = readRDS("./correlation.rds")

height_canopy = "10" # "0" or "10"

## Common climaticVariables
nbSpecies = corrR0_distEdge[, .N]

## Load tsn for shade tolerance
tsn = readRDS("../growth/tsn.rds")[, .(species, tolLevel)]
corrR0_distEdge = corrR0_distEdge[tsn, on = "species"]

## Sort by shade tolerance and by species name within shade tolerance, prepare colour
setorderv(x = corrR0_distEdge, cols = "species")
corrR0_distEdge[tolLevel == "L", colour := "#00F9FF"]
corrR0_distEdge[tolLevel == "M", colour := "#0E51FF"]
corrR0_distEdge[tolLevel == "H", colour := "#010120"]

shadeL = col2rgb("#00F9FF")[,1]
shadeM = col2rgb("#0E51FF")[,1]
shadeH = col2rgb("#010120")[,1]

#### Plot
tikz(paste0("3correlations_proj", height_canopy, ".tex"), width = 6.5, height = 5)
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

	tikzAnnotate(paste0("\\node (spPos) at (species_coords) {", corrR0_distEdge[i, species], "};"))

	# Add points
	correl_north_proj = paste0("correl_north_proj_", height_canopy)
	north_proj = corrR0_distEdge[i, ..correl_north_proj]
	points(x = i - 0.15, y = north_proj, pch = 24, lwd = 0, # Triangle up, filled, no border
		bg = corrR0_distEdge[i, colour])

	correl_tot = paste0("correl_tot_", height_canopy)
	corr_tot = corrR0_distEdge[i, ..correl_tot]
	points(x = i, y = corr_tot, pch = 21, lwd = 0, # Circle, filled, no border
		bg = corrR0_distEdge[i, colour])

	correl_south_proj = paste0("correl_south_proj_", height_canopy)
	south_proj = corrR0_distEdge[i, ..correl_south_proj]
	points(x = i + 0.15, y = south_proj, pch = 25, lwd = 0, # Triangle down, filled, no border
		bg = corrR0_distEdge[i, colour])
}

axis(side = 2, at = seq(-1, 1, by = 0.25),
	labels = c("-1.00", "", "-0.5", "", "0", "", "0.5", "", "1.00"))

## Legend
tikzAnnotate(paste0("\\definecolor{shadeL}{RGB}{", paste(shadeL, collapse = ","), "}"))
tikzAnnotate(paste0("\\definecolor{shadeM}{RGB}{", paste(shadeM, collapse = ","), "}"))
tikzAnnotate(paste0("\\definecolor{shadeH}{RGB}{", paste(shadeH, collapse = ","), "}"))

tikzAnnotate("
	\\matrix [below right] at (current bounding box.north west) {
		\\node [shape = rectangle, fill = shadeL, label = right:Low] {}; &
		\\node [shape = rectangle, fill = shadeM, label = right:Medium] {}; &
		\\node [shape = rectangle, fill = shadeH, label = right:High] {}; \\\\
	};
")
dev.off()
