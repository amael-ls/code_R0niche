#### Aim of prog: correlation between R0 and the distance to closest edge
## Load data
#	- Little's range (1971)
#	- R0 calculated by Matlab
#
## Distance closest edge
#	- Distance between a point and its orthogonal projection on Little's edge
#
## Plot
#	- tikzDevice for latex

#### Load package and clear memory
library(data.table)
library(tikzDevice)
library(doParallel)
library(stringi)
library(raster)
library(sf)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Create the cluster
## Cluster variables
array_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste0("array id = ", array_id))

nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

print("nodeslist")
print(nodeslist)
print("end nodeslist")

## Make cluster
cl = makeCluster(nodeslist, type = "PSOCK")
registerDoParallel(cl)

print("cluster done")

#### Common variables
ls_14species = dir(path = "./results", pattern = "[0-9]{4,}")
(species = ls_14species[array_id])

#### Compute correlation R0 -- dist
## Load data
# Shapefile
shpPath = "../little1971/"
little1971 = st_read(dsn = paste0(shpPath, species))
little1971 = st_transform(little1971,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

little1971 = st_geometry(st_union(little1971))
little1971 = st_cast(little1971, "MULTILINESTRING")

# R0 and coordinates
clim_2010 = readRDS(paste0("./results/", species, "/lonLatR0cropped.rds"))

## Compute distance and correlation
# Distance
st_crs(clim_2010) == st_crs(little1971)
dist = st_distance(clim_2010, little1971, by_element = TRUE)
saveRDS(dist, paste0("./results/", species, "/distance.rds"))

# Coerce to data table and correlation
st_geometry(clim_2010) = NULL
correl_R0edge = cor(clim_2010[, R0], dist)
saveRDS(correl_R0edge, paste0("./results/", species, "/correl_R0_distEdge.rds"))

print(paste0(species, " done. Correlation = ", correl_R0edge))
