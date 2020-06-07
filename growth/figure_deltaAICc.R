
#### Aim of prog: Generate an extra figure to show the impact of some explanatory variables
## I compare only the folowing models
#		- formula2 (glmm_ReML=FALSE.R), T_a, P_a, dbh, cs, quadratic terms. (Selected model), model 2
#		- formula2 (glmm_fixef), T-a, P_a, quadratic terms, model 18
#		- formula3 (glmm_fixef), cs, model 19
#		- formula6 (glmm_fixef), dbh, quadratic terms, model 22
#

#### Load package and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(tikzDevice)

#### Common variables
## List folders
ls_folders = list.files(path = "./", pattern = "^array_[0-9]")
nbSpecies = length(ls_folders)
modelsToKeep = paste0("model", c(2, 18, 19, 22)) # cf comments introduction

mods = data.table(keptModels = modelsToKeep, variables = c("selected", "climate", "competition", "dbh"))

## Read the tables containing the Delta_AICc and R², 25 rows (i.e., models) per species
ls_delta_aic = vector(mode = "list", length = nbSpecies)

for (i in 1:nbSpecies)
{
	ls_delta_aic[[i]] = setDT(readRDS(paste0("./", ls_folders[i], "/aic_mega.rds")))
	ls_delta_aic[[i]][, Modnames := as.character(Modnames)]
	ls_delta_aic[[i]][, model := paste0(i, "_", Modnames)]
	ls_delta_aic[[i]] = ls_delta_aic[[i]][Modnames %in% modelsToKeep]
	ls_delta_aic[[i]][, rel_Delta_AICc := Delta_AICc/min(Delta_AICc)]
}

ls_delta_aic = rbindlist(ls_delta_aic)

## Select only the models abovementionned
ls_delta_aic = ls_delta_aic[, .(Modnames, AICc, Delta_AICc, rel_Delta_AICc, model)]
ls_delta_aic[, log10_Delta_AICc := log10(rel_Delta_AICc)]

ls_delta_aic[, average := mean(log10_Delta_AICc), by = Modnames]

x_pos = unique(ls_delta_aic[, .(Modnames, average)])
setorderv(x_pos, "average", -1)
setnames(x_pos, old = "Modnames", new = "keptModels")
x_pos[, rank := 1:.N]
x_pos[, x := 2*rank - 1]

x_pos = x_pos[mods, on = "keptModels", nomatch = 0]

ls_delta_aic[, x := x_pos[keptModels == Modnames, x], by = Modnames]

m_x = min(ls_delta_aic[, x]) - 0.5
M_x = max(ls_delta_aic[, x]) + 0.5

m_y = floor(min(ls_delta_aic[, log10_Delta_AICc])) - 0.05
M_y = ceiling(max(ls_delta_aic[, log10_Delta_AICc]))

tikz('figure_deltaAICc.tex', width = 5, height = 5)
plot(ls_delta_aic[, .(x, log10_Delta_AICc)], xlab = "", ylab = "", main = "",
    yaxt = "n", xaxt = "n", yaxs="i", xaxs="i", lwd = 1, xlim = c(m_x, M_x), ylim = c(m_y, M_y),
	type = "p", pch = 19, col = rgb(30, 30, 30, maxColorValue = 255))
points(unique(ls_delta_aic[, .(x, average)]), col = "#FFAD33", pch = 19)

axis(1, at = x_pos[, x], labels = x_pos[, variables])
axis(2, at = 0:M_y, labels = 0:M_y, las = 2)

title(xlab = "Model", ylab = "$ \\varphi $")
dev.off()
