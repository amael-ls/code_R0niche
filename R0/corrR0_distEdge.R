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
#
## Remarks
# 1 - Note that geosphere::dist2line, is computationally intensive
#
# 2 - I had to modify the function alongTrackDistance from the geosphere package:
# see the closed issue https://github.com/rspatial/geosphere/issues/3
# The subtlety comes from the fact that acos (which is the recirpocal of cos) is
# defined on [-1, 1]. Depending on how the computer represent numbers, the function
# geosphere::alongTrackDistance threw an error because it creates numbers slightly
# bigger than 1 (or smaller than -1). By slightly I mean a difference of ~2e-16 for
# instance...

#### Load package and clear memory
library(data.table)
library(doParallel)
library(geosphere)
library(stringi)
library(raster)
library(sf)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Tool functions
source("./dist_geosphere.R")

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

#### Compute distance
## Load data
# Shapefile, do not crop it otherwise it would create artificial boundaries!
shpPath = "../little1971/"
little1971 = st_read(dsn = paste0(shpPath, species))
little1971 = st_transform(little1971,
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

little1971 = st_geometry(st_union(little1971))
little1971 = st_cast(little1971, "MULTILINESTRING")
little1971 = as(little1971, "Spatial")

# R0 and coordinates with a competition, canopy height s* = 10 m
clim_2010 = stack(paste0("./results/", species, "/10m/lonLatR0cropped.grd"))
clim_2010 = projectRaster(clim_2010,
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

## Compute distance
identical(st_crs(clim_2010), st_crs(little1971))

clim_2010 = data.table(rasterToPoints(clim_2010), keep.rownames = TRUE)
coords = as.matrix(clim_2010[!is.na(R0), .(x, y)])
lim = nrow(coords)

print(paste0("Number of coordinates to treate: ", lim))

start = Sys.time()
dist_results = foreach(i = 1:lim, .packages = c("geosphere", "sp")) %dopar%
	my_dist2Line(coords[i,], little1971, distfun = geosphere::distGeo) # It is using my version of dist2line, cf remarks
end = Sys.time()
end - start

dist_results = do.call(rbind, dist_results)
dist_results = data.table(dist_results)
saveRDS(dist_results, paste0("./results/", species, "/distance.rds"))

#### Compute correlation between R0_10m and dist to closest edge (s* = 10m)
## Total
correl_R0edge = cor(clim_2010[, R0], dist_results[, distance])
saveRDS(correl_R0edge, paste0("./results/", species, "/correl_R0_distEdge.rds"))

centroid = clim_2010[, lapply(.SD, mean), .SDcols = c("x", "y")] # Average lon/lat of the cropped data

## Northern part
# Northern region, i.e., latitudes northern to centroid
ind_north = which(clim_2010[, y] >= centroid[, y])
correl_R0edge_north = cor(clim_2010[ind_north, R0], dist_results[ind_north, distance])
saveRDS(correl_R0edge_north, paste0("./results/", species, "/correl_R0_distEdge_north.rds"))

# Points with a projection in the north (i.e., closer to north than south)
ind_north = which(dist_results[, lat] >= centroid[, y])
correl_R0edge_north = cor(clim_2010[ind_north, R0], dist_results[ind_north, distance])
saveRDS(correl_R0edge_north, paste0("./results/", species, "/correl_R0_distEdge_north_proj.rds"))

## Southern part
# Southern region, i.e., latitudes southern to centroid
ind_south = which(clim_2010[, y] < centroid[, y])
correl_R0edge_south = cor(clim_2010[ind_south, R0], dist_results[ind_south, distance])
saveRDS(correl_R0edge_south, paste0("./results/", species, "/correl_R0_distEdge_south.rds"))

# Points with a projection in the south (i.e., closer to south than north)
ind_south = which(dist_results[, lat] < centroid[, y])
correl_R0edge_south = cor(clim_2010[ind_south, R0], dist_results[ind_south, distance])
saveRDS(correl_R0edge_south, paste0("./results/", species, "/correl_R0_distEdge_south_proj.rds"))

################################################################################
######## PART II, without competition
#### Load data
## R0 without competition, s* = 0m
clim_2010 = stack(paste0("./results/", species, "/0m/lonLatR0cropped.grd"))
clim_2010 = projectRaster(clim_2010,
	crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

clim_2010 = data.table(rasterToPoints(clim_2010), keep.rownames = TRUE)

#### Compute correlation between R0_10m and dist to closest edge (s* = 10m)
## Total
correl_R0edge = cor(clim_2010[, R0], dist_results[, distance])
saveRDS(correl_R0edge, paste0("./results/", species, "/correl_R0_distEdge_0m.rds"))

centroid = clim_2010[, lapply(.SD, mean), .SDcols = c("x", "y")] # Average lon/lat of the cropped data

## Northern part
# Northern region, i.e., latitudes northern to centroid
ind_north = which(clim_2010[, y] >= centroid[, y])
correl_R0edge_north = cor(clim_2010[ind_north, R0], dist_results[ind_north, distance])
saveRDS(correl_R0edge_north, paste0("./results/", species, "/correl_R0_distEdge_north_0m.rds"))

# Points with a projection in the north (i.e., closer to north than south)
ind_north = which(dist_results[, lat] >= centroid[, y])
correl_R0edge_north = cor(clim_2010[ind_north, R0], dist_results[ind_north, distance])
saveRDS(correl_R0edge_north, paste0("./results/", species, "/correl_R0_distEdge_north_proj_0m.rds"))

## Southern part
# Southern region, i.e., latitudes southern to centroid
ind_south = which(clim_2010[, y] < centroid[, y])
correl_R0edge_south = cor(clim_2010[ind_south, R0], dist_results[ind_south, distance])
saveRDS(correl_R0edge_south, paste0("./results/", species, "/correl_R0_distEdge_south_0m.rds"))

# Points with a projection in the south (i.e., closer to south than north)
ind_south = which(dist_results[, lat] < centroid[, y])
correl_R0edge_south = cor(clim_2010[ind_south, R0], dist_results[ind_south, distance])
saveRDS(correl_R0edge_south, paste0("./results/", species, "/correl_R0_distEdge_south_proj_0m.rds"))
