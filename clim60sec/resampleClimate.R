
#### Aim of prog: Resample the rasters that are misaligned with respect to the first (reference)
## Description
# This is done within a year for bioclim and pcp
# The program checkRefRaster.R verified that the references are the same among years within
# bioclim and pcp respectively (TRUE), and that a reference of bioclim is the same as pcp (FALSE)
# Given the pcp reference is not the same as bioclim reference, I decided to resample all the pcp
# on the bioclim reference (not that the maximum difference of coordinates is 0.0085 between the
# two references, which is almost 1 km difference)
#
## WARNING: I overrite the data, if you want to keep the original climate data,
# you need to copy them first

#### Load library and clear memory
library(doParallel)
library(raster)

options(max.print = 500)
rm(list = ls())

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

#### Read & resample data
## Common variables
allYear = 1956:2013
(year = allYear[array_id])
checkProj = TRUE

## Load path bioclim
(loadPath_bio = paste0("./climateData/bioclim/", year, "/"))
ls_files_bio = list.files(path = loadPath_bio, pattern = ".asc")
nb_raster_bio = length(ls_files_bio)

if (nb_raster_bio != 19)
	print(paste0("*** ERROR: bioclimatic variables missing, year = ", year))

## Load path pcp
(loadPath_pcp = paste0("./climateData/pcp/", year, "/"))
ls_files_pcp = list.files(path = loadPath_pcp, pattern = ".asc")
nb_raster_pcp = length(ls_files_pcp)

if (nb_raster_pcp != 12)
	print(paste0("*** ERROR: pcp month missing, year = ", year))

## Reference raster /!\ (bio & pcp) /!\
ref_raster = raster(paste0(loadPath_bio, "bio60_01.asc"))
ref_coordinates = coordinates(ref_raster)
ref_proj = crs(ref_raster)

## Resample bioclim if required
for (i in 2:nb_raster_bio) # Starts at 2, the reference is aligned by definition
{
	# Load current raster
	current_raster = raster(paste0(loadPath_bio, ls_files_bio[i]))
	coord = coordinates(current_raster)
	proj = crs(current_raster)

	# Check coordinates and projections are the same
	booleanCoord = isTRUE(all.equal(ref_coordinates, coord))
	booleanProj = isTRUE(all.equal(ref_proj, proj))

	# If not the same coordinates, resample
	if (!booleanCoord & booleanProj)
	{
		print(paste0("biovar ", i, " has different coordinates from biovar 1"))
		print(paste0("InfNorm(diff coordinates) = ", max(abs(coord - ref_coordinates))))

		raster::resample(x = current_raster, y = ref_raster, method = "ngb",
			filename = paste0(loadPath_bio, ls_files_bio[i]), overwrite = TRUE)
	}

	if (!booleanProj)
	{
		print(paste0("biovar ", i, " has a different projection from biovar 1"))
		checkProj = FALSE
	}
}

if(!checkProj)
{
	print(paste0("WARNING: projection (bioclim) for year = ", year, " is different from reference"))
	checkProj = TRUE
}

## Resample pcp if required
for (i in 1:nb_raster_pcp) # Starts at 1, we know the first pcp is misaligned
{
	# Load current raster
	current_raster = raster(paste0(loadPath_pcp, ls_files_pcp[i]))
	coord = coordinates(current_raster)
	proj = crs(current_raster)

	# Check coordinates and projections are the same
	booleanCoord = isTRUE(all.equal(ref_coordinates, coord))
	booleanProj = isTRUE(all.equal(ref_proj, proj))

	# If not the same coordinates, resample
	if (!booleanCoord & booleanProj)
	{
		print(paste0("pcpvar ", i, " has different coordinates from pcp 1"))
		print(paste0("InfNorm(diff coordinates) = ", max(abs(coord - ref_coordinates))))

		raster::resample(x = current_raster, y = ref_raster, method = "ngb",
			filename = paste0(loadPath_pcp, ls_files_pcp[i]), overwrite = TRUE)
	}

	if (!booleanProj)
	{
		print(paste0("pcpvar ", i, " has a different projection from biovar 1"))
		checkProj = FALSE
	}
}

if(!checkProj)
	print(paste0("WARNING: projection (pcp) for year = ", year, " is different from reference"))
