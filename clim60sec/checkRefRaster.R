
#### Aim of prog: Check all the first bioclim have the same coordinates along years
# This is important because first bioclim will then be the reference for resampling
# when required. Because I do it in parallel (for years), I first must be sure the
# chosen reference (i.e., 1st bioclim) is independent of the chosen year.
# Same with pcp

#### Load library and clear memory
library(raster)

options(max.print = 500)
rm(list = ls())

#### Read data
## Common variables
allYear = 1956:2013
nb_years = length(allYear)
ref_raster_bio = raster(paste0("./climateData/bioclim/", allYear[1], "/bio60_01.asc"))
ref_coordinates_bio = coordinates(ref_raster_bio)

ref_raster_pcp = raster(paste0("./climateData/pcp/", allYear[1], "/pcp60_01.asc"))
ref_coordinates_pcp = coordinates(ref_raster_pcp)

## Check the two references are equal
comparison2Ref = all.equal(ref_coordinates_bio, ref_coordinates_pcp)
if (!isTRUE(comparison2Ref))
{
	print(comparison2Ref)
	print("Need to choose the bioclim reference for pcp")
	print(paste0("InfNorm(diff coordinates) = ", max(abs(ref_coordinates_bio - ref_coordinates_pcp))))
}

for (i in 2:nb_years)
{
	print(paste0("year: ", allYear[i]))

	# Check within bioclim
	current_raster = raster(paste0("./climateData/bioclim/", allYear[i], "/bio60_01.asc"))
	coords = coordinates(current_raster)
	print(all.equal(ref_coordinates_bio, coords))

	# Check within pcp
	current_raster = raster(paste0("./climateData/pcp/", allYear[i], "/pcp60_01.asc"))
	coords = coordinates(current_raster)
	print(all.equal(ref_coordinates_pcp, coords))
}
