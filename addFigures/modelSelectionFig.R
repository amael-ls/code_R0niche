
#### Aim of prog: plot on the same figure phi and psi defined in the parameterisation section
## For growth, I compare only 4 models
#		- formula2 (glmm_ReML=FALSE.R), T_a, P_a, dbh, cs, quadratic terms. (Selected model), model 2
#		- formula2 (glmm_fixef), T_a, P_a, quadratic terms, model 18
#		- formula3 (glmm_fixef), cs, model 19
#		- formula6 (glmm_fixef), dbh, quadratic terms, model 22
#
## For mortality, I compare only 4 models
#		- formula7 (selectClimaticVariables.R), T_m, P_m, dbh, cs, quadratic terms
#		- formula2 (submodels), T_m, P_m, quadratic terms, submodel 2
#		- formula3 (submodels), cs, submodel 3
#		- formula6 (submodels), dbh, quadratic terms, submodel 6
#

#### Load package and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(tikzDevice)
library(stringi)

######## Part I: Growth
#### Common variables
## List folders
growthPath = "../growth/"
ls_folders = list.files(path = growthPath, pattern = "^array_[0-9]")
nbSpecies = length(ls_folders)
modelsToKeep = paste0("model", c(2, 18, 19, 22)) # cf comments introduction

mods = data.table(keptModels = modelsToKeep, variables = c("best model", "climate", "competition", "dbh"))

## Read the tables containing the Delta_AICc and R², 25 rows (i.e., models) per species
ls_delta_aic = vector(mode = "list", length = nbSpecies)

for (i in 1:nbSpecies)
{
	ls_delta_aic[[i]] = setDT(readRDS(paste0(growthPath, ls_folders[i], "/aic_mega.rds")))
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

x_pos_growth = unique(ls_delta_aic[, .(Modnames, average)])
setorderv(x_pos_growth, "average", -1)
setnames(x_pos_growth, old = "Modnames", new = "keptModels")
x_pos_growth[, rank := 1:.N]
x_pos_growth[, x := 2*rank - 1]

x_pos_growth = x_pos_growth[mods, on = "keptModels", nomatch = 0]

ls_delta_aic[, x := x_pos_growth[keptModels == Modnames, x], by = Modnames]

# Shift x pos to the left, to avoid overlaping with mortality
ls_delta_aic[, x := x - 0.25]

# Axes lim
m_x_growth = min(ls_delta_aic[, x]) - 0.5
M_x_growth = max(ls_delta_aic[, x]) + 0.5

m_y_growth = floor(min(ls_delta_aic[, log10_Delta_AICc])) - 0.05
M_y_growth = ceiling(max(ls_delta_aic[, log10_Delta_AICc]))

######## Part II: Mortality
#### Common variables
## List folders
mortalityPath = "../mortality/"
ls_folders = list.files(path = mortalityPath, pattern = "^array_[0-9]")
nbSpecies = length(ls_folders)
modelsToKeep = c("model7", paste0("submodel", c(2, 3, 6))) # cf comments introduction

mods = data.table(keptModels = modelsToKeep, variables = c("best model", "climate", "competition", "dbh"))

## Read the tables containing the WAIC, 14 rows (i.e., models) per species
ls_delta_waic = vector(mode = "list", length = nbSpecies)

for (i in 1:nbSpecies)
{
	ls_delta_waic[[i]] = setDT(readRDS(paste0(mortalityPath, ls_folders[i], "/waic_mega.rds")))
	ls_delta_waic[[i]][, model_id := paste0(i, "_", model)]
	ls_delta_waic[[i]] = ls_delta_waic[[i]][model %in% modelsToKeep]
	ls_delta_waic[[i]][, delta_WAIC := waic - min(waic)]
	ls_delta_waic[[i]][, sp_id := stri_sub(str = ls_folders[i],
		from = stri_locate_first(str = ls_folders[i], regex = "[0-9]")[1])]
}

ls_delta_waic = rbindlist(ls_delta_waic)

ls_delta_waic[, log10_delta_waic := log10(delta_WAIC + 1)]

ls_delta_waic[, average := mean(log10_delta_waic), by = model]

x_pos_mortality = unique(ls_delta_waic[, .(model, average)])
setorderv(x_pos_mortality, "average", -1)
setnames(x_pos_mortality, old = "model", new = "keptModels")
x_pos_mortality[, rank := 1:.N]
x_pos_mortality[, x := 2*rank - 1]

x_pos_mortality = x_pos_mortality[mods, on = "keptModels", nomatch = 0]

namesToModify = names(x_pos_mortality)[!(names(x_pos_mortality) == "variables")]
setnames(x_pos_mortality,
	old = namesToModify,
	new = paste0(namesToModify, "_mortality"))

x_pos_mortality = x_pos_mortality[x_pos_growth, on = "variables"]
setorderv(x_pos_mortality, "rank", +1)

ls_delta_waic[, x := x_pos_mortality[keptModels_mortality == model, x], by = model]

# Shift x pos to the right, to avoid overlaping with growth
ls_delta_waic[, x := x + 0.25]

# Axes lim
m_x_mortality = min(ls_delta_waic[, x]) - 0.5
M_x_mortality = max(ls_delta_waic[, x]) + 0.5

m_y_mortality = floor(min(ls_delta_waic[, log10_delta_waic])) - 0.05
M_y_mortality = ceiling(max(ls_delta_waic[, log10_delta_waic]))

######## Part III: lot figure
#### Open figure
tikz('aic_waic.tex', width = 5, height = 5)

## add extra space to right margin of plot within frame
par(mar = c(5, 4, 4, 6) + 0.1)

## Growth
plot(ls_delta_aic[, .(x, log10_Delta_AICc)], xlab = "", ylab = "", main = "",
	yaxt = "n", xaxt = "n", yaxs="i", xaxs="i", lwd = 1, xlim = c(0, 8),
	ylim = c(m_y_growth, M_y_growth),
	type = "p", pch = 19, col = rgb(30, 30, 30, maxColorValue = 255))
points(unique(ls_delta_aic[, .(x, average)]), col = "#FFAD33", pch = 19)

axis(1, at = x_pos_growth[, x], labels = x_pos_growth[, variables])
axis(2, at = 0:M_y_growth, labels = 0:M_y_growth, las = 2)

mtext("Performance growth model ($ \\varphi $)", side = 2, line = 2.5)

## Allow a second plot on the same graph
par(new = TRUE)

plot(ls_delta_waic[, .(x, log10_delta_waic)], xlab = "", ylab = "", main = "",
	yaxt = "n", xaxt = "n", yaxs="i", xaxs="i", lwd = 1, xlim = c(0, 8),
	ylim = c(m_y_mortality, M_y_mortality),
	type = "p", pch = 17, col = rgb(30, 30, 30, maxColorValue = 255))
points(unique(ls_delta_waic[, .(x, average)]), col = "#00F9FF", pch = 17)

axis(4, at = 0:M_y_mortality, labels = 0:M_y_mortality, las = 2)
mtext("Performance mortality model ($ \\psi $)", side = 4, line = 2.5)

par(xpd = NA)
legend(0, 4.5, c("growth ($ G $)", "mortality ($ \\mu $)"), pch = c(19, 17),
	bty = "n", horiz = TRUE)

dev.off()
