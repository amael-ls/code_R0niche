
#### Aim of prog: plot correlations R0 <--> random forest versus R0 <--> presAbs

#### Load package and clear memory
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data
## Correlation R0 -- competition, randomForest
corrR0_rf = readRDS("./correlation.rds") #[, .(sp_code, correl_0m, correl_10m)]

## Correlation R0 -- competition, presAbs
corrR0_pa = readRDS("correlation_R0_presAbsData.rds")

## tsn
tsn = readRDS("../growth/tsn.rds")[, .(species, tolLevel)]
setnames(tsn, old = "species", new = "sp_code")

corrR0_rf = corrR0_rf[tsn, on = "sp_code"]
corrR0_pa = corrR0_pa[tsn, on = "sp_code"]

#### Plot
## Set colours
corrR0_rf[tolLevel == "L", colour := "#00F9FF"]
corrR0_rf[tolLevel == "M", colour := "#0E51FF"]
corrR0_rf[tolLevel == "H", colour := "#010120"]

shadeL = col2rgb("#00F9FF")[,1]
shadeM = col2rgb("#0E51FF")[,1]
shadeH = col2rgb("#010120")[,1]

## Tikz plot
tikz('./correl_rf_vs_pa.tex', width = 5, height = 5)
op = par(mar = c(3, 3, 2, 2), mgp = c(2, 0.75, 0),
	oma = c(0,0,0.9,0), tck = -0.015, las = 1)

plot(x = NULL, y = NULL, xlim = c(-1, 1),
	ylim = c(-1, 1), axes = TRUE, bg = "transparent",
	xlab = "Correlations (presence/absence data)", ylab = "Correlations (random forest)")

points(x = corrR0_pa[, correlR0_presAbs_0m], y = corrR0_rf[, correl_0m], pch = 17, # Triangle
	col = corrR0_rf[, colour])

points(x = corrR0_pa[, correlR0_presAbs_10m], y = corrR0_rf[, correl_10m], pch = 16, # Circle
	col = corrR0_rf[, colour])

# Add identity line
abline(a = 0, b = 1, lwd = 2)

tikzAnnotate(paste0("\\definecolor{shadeL}{RGB}{", paste(shadeL, collapse = ","), "}"))
tikzAnnotate(paste0("\\definecolor{shadeM}{RGB}{", paste(shadeM, collapse = ","), "}"))
tikzAnnotate(paste0("\\definecolor{shadeH}{RGB}{", paste(shadeH, collapse = ","), "}"))

tikzCoord(-1, 1, "posLegend")

tikzAnnotate("
	\\matrix [right] at (posLegend) {
		\\node [shape = rectangle, fill = shadeL, label = right:Low] {}; &
		\\node [shape = rectangle, fill = shadeM, label = right:Medium] {}; &
		\\node [shape = rectangle, fill = shadeH, label = right:High] {}; \\\\
	};
")
dev.off()
