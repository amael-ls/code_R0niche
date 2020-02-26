# Read me of clim60sec

## Order of the programs

The programs must be run in this order:

1. downloadBioClimData.R
2. downloadPrecipData.R
3. checkRefRaster.R
4. resampleClimate.R
5. subsetClimForTrees.R
6. createClimDatabase.R
7. clim2010.R

## Remarks

There are 19 bioclimatic variables, named bio60_01, ..., bio60_19.
Further information about them can be found at:

0. https://pubs.usgs.gov/ds/691/ds691.pdf
1. https://fennerschool.anu.edu.au/files/anuclim61.pdf

*** Section 6.1, page 61
The  descriptions  below  assume that BIOCLIM is using the weekly time step (the default).
If the monthly time step is selected, monthly values are used when calculating these parameters.
The quarterly parameters are not aligned to any calendar quarters. BIOCLIM's definition of a
quarter is any 13 consecutive weeks, (or any consecutive 3 months if running with a monthly
time step). For example, the driest quarter will be the 13 consecutive weeks that are drier than
any other set of 13 consecutive weeks.

*** /!\ Section 1.4, page 10 /!\
When outputting results from ESOCLIM, BIOCLIMand GROCLIMof AN UCLIMVersion 5.1
to Arc/InfoUNGENERATE files (point data), Arc/InfoASCIIGRID or IDRISI ASCII image files,
the output values were multiplied by 10 or 100 and then rounded to the nearest integer.
This was done to preserve appropriate precision in the output data with smaller output files.
However this process led to some confusion and inconvenience, particularly when these outputs
were compared with data from other sources.

## Names of BIO variables (http://worldclim.org/bioclim)

- BIO1  = Annual Mean Temperature
- BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))
- BIO3  = Isothermality (BIO2/BIO7)
- BIO4  = Temperature Seasonality (Coefficient of Variation)
- BIO5  = Max Temperature of Warmest Month
- BIO6  = Min Temperature of Coldest Month
- BIO7  = Temperature Annual Range (BIO5-BIO6)
- BIO8  = Mean Temperature of Wettest Quarter
- BIO9  = Mean Temperature of Driest Quarter
- BIO10 = Mean Temperature of Warmest Quarter
- BIO11 = Mean Temperature of Coldest Quarter
- BIO12 = Annual Precipitation
- BIO13 = Precipitation of Wettest Month
- BIO14 = Precipitation of Driest Month
- BIO15 = Precipitation Seasonality (Coefficient of Variation)
- BIO16 = Precipitation of Wettest Quarter
- BIO17 = Precipitation of Driest Quarter
- BIO18 = Precipitation of Warmest Quarter
- BIO19 = Precipitation of Coldest Quarter
