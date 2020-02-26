
#### Aim of prog: plot the confidence intervals for each species
## Each plot will be for one parameter. It is easy to do it per species

#### Load package and clear memory
library(data.table)
library(tikzDevice)
library(stringi)

#### Load data and set common variables
## Uncertainty file
confInt = fread("./confInt_allSpecies_mortality.csv")

## Clean the names in term to avoid problems with tikz
# Remove opening brackets
confInt[stri_detect(str = term, regex = "\\("),
	term := stri_replace_all(str = term, replacement = "", regex = "\\(")]

# Remove closing brackets
confInt[stri_detect(str = term, regex = "\\)"),
	term := stri_replace_all(str = term, replacement = "", regex = "\\)")]

# Remove underscore (reserved to mathmode latex)
confInt[stri_detect(str = term, regex = "\\_"),
	term := stri_replace_all(str = term, replacement = "", regex = "\\_")]

# Remove hat
confInt[stri_detect(str = term, regex = "\\^"),
	term := stri_replace_all(str = term, replacement = "", regex = "\\^")]

## Common variables
nbSpecies = length(unique(confInt[, species]))
nbParameters = confInt[, .N]/nbSpecies # Should be 12
nbPlotsPerPage = 4
nbPages = nbParameters %/% nbPlotsPerPage # Number of full pages
lastPlot = nbParameters %% nbPlotsPerPage

if (!dir.exists("./confIntPlot"))
	dir.create("./confIntPlot")

#### Plot
for (i in 1:nbPages)
{
	tikz(paste0("./confIntPlot/page_", i, "_mortality.tex"), width = 3.1, height = 3.1)
	for (j in 1:nbPlotsPerPage)
	{
		## Extract all the species-specific values for params (i - 1)*nbPlotsPerPage + j
		index = seq((i - 1)*nbPlotsPerPage + j, confInt[, .N], by = nbParameters)
		sub_dt = confInt[index]
		sub_dt[, species_id := stri_sub(str = species,
			from = stri_locate_first(str = species, regex = "-")[,1] + 1)]

		graphName = unique(sub_dt[, term])

		# Check there is only one parameter
		if (length(graphName) != 1)
		{
			print("*** Error (from sub_dt) ***: there should be only one parameter")
			print(paste0("page ", i, ", iteration ", j, " skipped; names = ", graphName))
			next;
		}

		# Decreasing alphabetical order
		setorderv(sub_dt, "species_id", -1)

		## Plot
		# Bounds of plot for climVar
		min_bound = min(sub_dt[, min_5p])
		max_bound = max(sub_dt[, max_95p])

		op = par(mar = c(2.5, 5, 0.1, 0.4), mgp = c(1.5, 0.2, 0),
			oma = c(0,0,0.9,0), tck = -0.015, las = 1)
		plot(x = NULL, y = NULL, xlim = c(min_bound, max_bound),
			ylim = c(1, nbSpecies), axes = FALSE, xlab = paste0(graphName, " estimation"),
			ylab = "")

		# if (stri_detect(str = climVar[i], regex = "precipitation"))
		# 	pos = round(seq(min_bound, max_bound, length.out = 5), 0)
		# if (stri_detect(str = climVar[i], regex = "temperature"))
		# 	pos = round(seq(min_bound, max_bound, length.out = 5), 1)

		pos = sort(c(round(min_bound, 2), round(max_bound, 2), 0))
		axis(side = 1, at = pos)

		if (j %% 2 == 1)
			axis(side = 2, at = seq(1, nbSpecies, length.out = nbSpecies),
				labels = sub_dt[, species_id], tck = 0, lwd = "")
		for (k in 1:nbSpecies)
		{
			m = unlist(sub_dt[k, min_5p])
			M = unlist(sub_dt[k, max_95p])
			a = unlist(sub_dt[k, median])

			segments(x0 = m, y0 = k, x1 = M, y1 = k, lwd = 2)
			points(x = a, y = k, pch = 20)
		}

		abline(v = 0, lwd = 2, lty = 5)
	}
	dev.off()
}

if (lastPlot != 0)
{
	tikz(paste0("./confIntPlot/page_", nbPages + 1, "_mortality.tex"), width = 3.1, height = 3.1)
	for (j in 1:lastPlot)
	{
		## Extract all the species-specific values for params (i - 1)*nbPlotsPerPage + j
		index = seq((i - 1)*nbPlotsPerPage + j, confInt[, .N], by = nbParameters)
		sub_dt = confInt[index]
		sub_dt[, species_id := stri_sub(str = species,
			from = stri_locate_first(str = species, regex = "-")[,1] + 1)]

		graphName = unique(sub_dt[, term])

		# Check there is only one parameter
		if (length(graphName) != 1)
		{
			print("*** Error (from sub_dt) ***: there should be only one parameter")
			print(paste0("page ", i, ", iteration ", j, " skipped; names = ", graphName))
			next;
		}

		# Decreasing alphabetical order
		setorderv(sub_dt, "species_id", -1)

		## Plot
		# Bounds of plot for climVar
		min_bound = min(sub_dt[, min_5p])
		max_bound = max(sub_dt[, max_95p])

		op = par(mar = c(2.5, 5, 0.1, 0.4), mgp = c(1.5, 0.2, 0),
			oma = c(0,0,0.9,0), tck = -0.015, las = 1)
		plot(x = NULL, y = NULL, xlim = c(min_bound, max_bound),
			ylim = c(1, nbSpecies), axes = FALSE, xlab = paste0(graphName, " estimation"),
			ylab = "")

		# if (stri_detect(str = climVar[i], regex = "precipitation"))
		# 	pos = round(seq(min_bound, max_bound, length.out = 5), 0)
		# if (stri_detect(str = climVar[i], regex = "temperature"))
		# 	pos = round(seq(min_bound, max_bound, length.out = 5), 1)

		pos = sort(c(round(min_bound, 2), round(max_bound, 2), 0))
		axis(side = 1, at = pos)

		if (j %% 2 == 1)
			axis(side = 2, at = seq(1, nbSpecies, length.out = nbSpecies),
				labels = sub_dt[, species_id], tck = 0, lwd = "")
		for (k in 1:nbSpecies)
		{
			m = unlist(sub_dt[k, min_5p])
			M = unlist(sub_dt[k, max_95p])
			a = unlist(sub_dt[k, median])

			segments(x0 = m, y0 = k, x1 = M, y1 = k, lwd = 2)
			points(x = a, y = k, pch = 20)
		}

		abline(v = 0, lwd = 2, lty = 5)
	}
	dev.off()
}
