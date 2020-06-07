
#### Aim of prog: plot MCMCs of mortality's parameters
#

#### Load packages
library(tikzDevice)
library(rstanarm)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool function
## Row bind a list of vectors to a single vector
rbindVec_ls = function(ls)
{
	# Allocate size of vector
	n = length(ls)
	size = 0
	for (i in 1:n)
		size = size + length(ls[[i]])

	results = numeric(size)
	first_index = 1

	for (i in 1:n)
	{
		last_index = first_index + length(ls[[i]]) - 1
		results[first_index:last_index] = ls[[i]]
		first_index = last_index + 1
	}
	return(results)
}

#### Common variables
## List array folders
ls_folders = dir(path = "./", pattern = "array_[0-9]{1,}")

## Chosen model
chosenModel = "model_7"

## List to store Rhat
rhat_ls = vector(mode = "list", length = length(ls_folders))
names(rhat_ls) = ls_folders

#### Plots
for (folder in ls_folders)
{
	savePath = paste0(folder, "/mcmc_plots/")
	if (!dir.exists(savePath))
		dir.create(savePath)

	species = dir(path = folder, pattern = "^[0-9]{4,}.*.txt$")
	species = stri_sub(str = species,
		to = stri_locate(str = species, regex = ".txt")[1] - 1)

	model = readRDS(paste0(folder, "/", chosenModel, ".rds"))

	rhat_ls[[folder]] = bayesplot::rhat(model)

	coeffToPlots = names(model$coefficients)
	coeffToPlots = coeffToPlots[!stri_detect(str = coeffToPlots,
		regex = "^b\\[\\(Intercept\\) plot_id")]

	jpeg(paste0(savePath, species, "_hist.jpg"), height = 1080, width = 1080, quality = 100)
	print(plot(x = model, plotfun = "mcmc_hist", pars = coeffToPlots))
	dev.off()

	jpeg(paste0(savePath, species, "_traces.jpg"), height = 1080, width = 1080, quality = 100)
	print(plot(x = model, plotfun = "trace", pars = coeffToPlots))
	dev.off()

	print(paste0("species ", species, " done"))
}

#### Plot the distribution of rhat all parameters
## rbind rhat
rhat = rbindVec_ls(rhat_ls)
(min(rhat))
(max(rhat))

## Plot
tikz("rhat_distrib.tex", width = 3.1, height = 3.1)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
hist(rhat, main = "", xlab = "$ \\hat{r} $")
dev.off()
