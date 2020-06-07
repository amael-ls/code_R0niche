
#### Aim of prog: Compute the relative difference between correl_0m and correl_10m
## We compute the relative diffenrence for both correlations:
#		R0 <-> randomForest
#		R0 <-> presAbs data
#
## Remark:
# We set the reference to be correl_0m. We use absolute values on the denominator
# to keep the sign of the difference.

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)
library(xtable)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data
correlR0_presAbs = readRDS("./correlation_R0_presAbsData.rds")[, .(species, correlR0_presAbs_0m, correlR0_presAbs_10m)]
correlR0_randForest = readRDS("./correlation.rds")[, .(species, correl_0m, correl_10m)]
tsn_dt = readRDS("../growth/tsn.rds")[, .(species_id, tolLevel)]
setnames(tsn_dt, old = "species_id", new = "species")

#### Compute relative differences (percentage)
correlR0_presAbs[, relDiff_presAbs := (correlR0_presAbs_10m - correlR0_presAbs_0m)*100/abs(correlR0_presAbs_0m)]
correlR0_randForest[, relDiff_randF := (correl_10m - correl_0m)*100/abs(correl_0m)]

correlR0_randForest = correlR0_randForest[tsn_dt, on = "species"]

tikz(paste0("relDiff.tex"), width = 5.1, height = 5.1, standAlone = TRUE)
op = par(mar = c(2.5, 4, 0.2, 0.2), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(correlR0_randForest[, tolLevel], correlR0_randForest[, relDiff_randF],
	xlab = "Shade tolerance", ylab = "Relative difference")
dev.off()
