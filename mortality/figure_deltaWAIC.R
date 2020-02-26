
#### Aim of prog: Generate an extra figure to show the impact of some explanatory variables
## I compare only the folowing models
#		- formula7 (selectClimaticVariables.R), T_m, P_m, dbh, cs, quadratic terms
#		- formula2 (submodels), T-a, P_a, quadratic terms, submodel 2
#		- formula3 (submodels), cs, submodel 3
#		- formula6 (submodels), dbh, quadratic terms, submodel 6
#

#### Load package and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(tikzDevice)
library(stringi)

#### Common variables
## List folders
ls_folders = list.files(path = "./", pattern = "^array_[0-9]")
nbSpecies = length(ls_folders)
modelsToKeep = c("model7", paste0("submodel", c(2, 3, 6))) # cf comments introduction

mods = data.table(keptModels = modelsToKeep, variables = c("best model", "climate", "competition", "dbh"))

## Read the tables containing the WAIC, 14 rows (i.e., models) per species
ls_delta_waic = vector(mode = "list", length = nbSpecies)

for (i in 1:nbSpecies)
{
	ls_delta_waic[[i]] = setDT(readRDS(paste0("./", ls_folders[i], "/waic_mega.rds")))
	ls_delta_waic[[i]][, model_id := paste0(i, "_", model)]
	ls_delta_waic[[i]] = ls_delta_waic[[i]][model %in% modelsToKeep]
	ls_delta_waic[[i]][, delta_WAIC := waic - min(waic)]
	ls_delta_waic[[i]][, sp_id := stri_sub(str = ls_folders[i],
		from = stri_locate_first(str = ls_folders[i], regex = "[0-9]")[1])]
}

ls_delta_waic = rbindlist(ls_delta_waic)

ls_delta_waic[, log10_delta_waic := log10(delta_WAIC + 1)]

ls_delta_waic[, average := mean(log10_delta_waic), by = model]

x_pos = unique(ls_delta_waic[, .(model, average)])
setorderv(x_pos, "average", -1)
setnames(x_pos, old = "model", new = "keptModels")
x_pos[, rank := 1:.N]
x_pos[, x := 2*rank - 1]

x_pos = x_pos[mods, on = "keptModels", nomatch = 0]
setorderv(x_pos, "rank", +1)

ls_delta_waic[, x := x_pos[keptModels == model, x], by = model]

m_x = min(ls_delta_waic[, x]) - 0.5
M_x = max(ls_delta_waic[, x]) + 0.5

m_y = floor(min(ls_delta_waic[, log10_delta_waic])) - 0.05
M_y = ceiling(max(ls_delta_waic[, log10_delta_waic]))

pdf("figure_deltaWAIC.pdf", height = 8, width = 8)
par(mar = c(5,6,2,2))
plot(ls_delta_waic[, .(x, log10_delta_waic)], xlab = "", ylab = "", main = "",
    yaxt = "n", xaxt = "n", yaxs="i", xaxs="i", lwd = 1, xlim = c(m_x, M_x), ylim = c(m_y, M_y),
	type = "p", pch = 19)
lines(x_pos[, .(x, average)], type = "l", lwd = 3,
	col = "#5CA7DD")
points(unique(ls_delta_waic[, .(x, average)]), col = "#FFAD33", pch = 19, cex = 1.2)

axis(1, at = x_pos[, x], labels = x_pos[, variables], cex.axis = 1.5)
axis(2, at = 0:M_y, labels = 0:M_y, cex.axis = 1.5, las = 2)

title(xlab = "Model", ylab = bquote(log[10](Delta~"AIC"[c])), cex.lab = 1.5)
dev.off()

tikz('figure_deltaWAIC.tex', width = 5, height = 5)
plot(ls_delta_waic[, .(x, log10_delta_waic)], xlab = "", ylab = "", main = "",
    yaxt = "n", xaxt = "n", yaxs="i", xaxs="i", lwd = 1, xlim = c(m_x, M_x), ylim = c(m_y, M_y),
	type = "p", pch = 19, col = rgb(30, 30, 30, maxColorValue = 255))
lines(x_pos[, .(x, average)], type = "l", lwd = 2,
	col = "#5CA7DD")
points(unique(ls_delta_waic[, .(x, average)]), col = "#FFAD33", pch = 19)

axis(1, at = x_pos[, x], labels = x_pos[, variables])
axis(2, at = 0:M_y, labels = 0:M_y, las = 2)

title(xlab = "Model", ylab = "$ \\log_{10}(\\Delta \\text{WAIC} + 1) $")
dev.off()
