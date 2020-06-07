
library(data.table)

ls_species = list.files("./results/")
n = length(ls_species)

results = data.table(species_id = ls_species, min_time = numeric(n), max_time = numeric(n),
	mean_time = numeric(n), median_time = numeric(n))

for (folder in ls_species)
{
	time_integ = fread(paste0("./results/", folder, "/time_integ.csv"))
	results[species_id == folder, min_time := min(time_integ)]
	results[species_id == folder, max_time := max(time_integ)]
	results[species_id == folder, mean_time := mean(time_integ$V1)]
	results[species_id == folder, median_time := median(time_integ$V1)]
}

print(results)
