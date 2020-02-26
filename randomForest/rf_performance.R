
#### Aims of prog: Check the performance of the random forests
#

#### Load packages
library(randomForest)
library(data.table)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool function
getSpecies = function(str)
{
	# Clean str if ends by "_cal.rds"
	str = stri_sub(str[stri_detect(str, regex = "_cal.rds$")], to = stri_locate_last(str, regex = "_")[1] - 1)
	# Replace underscores by scores
	species = stri_replace_all(str = str, replacement = "-", regex = "_")

	# Switch the Taxonomic Serial Number (tsn) with species code
	tsn_pos = stri_locate_last(str = species, regex = "-")[1]
	tsn = stri_sub(str = species, from = tsn_pos + 1)
	sp_code = stri_sub(str = species, to = tsn_pos - 1)
	species = paste0(tsn, "-", sp_code)
	return (species)
}

#### Load the results
## List files finishing by _cal.rds from ./results/calibration/
(loadPath = "./results/calibration/")
ls_results = list.files(path = loadPath, pattern = "_cal.rds$")

n = length(ls_results)
performance_rf = data.table(species = getSpecies(ls_results),
	false_error = numeric(length = n),
	true_error = numeric(length = n))

for (i in 1:n)
{
	results = readRDS(paste0(loadPath, ls_results[i]))
	performance_rf[i, false_error := results["confusion"][[1]]["FALSE", "class.error"]]
	performance_rf[i, true_error := results["confusion"][[1]]["TRUE", "class.error"]]
	print(paste0("species: ", performance_rf[i, species], " done"))
}

saveRDS(performance_rf, "./performance.rds")

# all.equal(4910/(42814 + 4910), results["confusion"][[1]][1, "class.error"])
