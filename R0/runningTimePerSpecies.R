
#### Aim of prog: Count how many pixels per species
## Comment:
# This prog helps me to estimate how much time do I need to compute R0

library(data.table)

ls_folder = dir("./Matlab_data/")
n = length(ls_folder)

results = data.table(species = ls_folder, dim = numeric(n),
	timeBeluga = numeric(n), timeCedar = numeric(n))

for (i in 1:n)
{
	dt = fread(paste0("./Matlab_data/", ls_folder[i], "/matlabGrowth_above.csv"))
	results[i, dim := dt[, .N]]
}

results[, timeBeluga := dim*2935/(201119*3600)] # Time it took for Picea rubens, with 32 CPUs on Beluga
results[, timeCedar := dim*2345/(201119*3600)] # Time it took for Picea rubens, with 48 CPUs on Cedar

results[, timeCedar_ACRU := dim*15634/(1641791*3600)] # Time it took for Acer Rubrum, with 48 CPUs on Cedar
results[, timeCedar_FAGR := dim*11089/(1162143*3600)] # Time it took for Fagus grandifolia, with 48 CPUs on Cedar

results[, timeCedar_ACRU - timeCedar] # Difference in hours
results[, (timeCedar_ACRU - timeCedar_FAGR)*3600] # To get the difference in seconds
