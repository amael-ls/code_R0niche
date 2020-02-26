# Execution order of the programs

## Mandatory programs

1. Clean the raw tree data, using **treeCleaner.R**
2. Slice the tree dataset using **slicer.R**
3. Calculate the competition for each slice (the dead trees are kept, but their height is 0, hence dead trees do not play any role in the competition for light)) using **calculateCompetition.R**
4. Check all the slice have been done using **listFiles.R**
5. Reconstruct the tree database from the slices using **rebuild.R**
6. Create growth and mortality data, using **createGrowthMortality.R**:
   - Merge tree data with climate (already formatted)
   - Create growth data
   - Create mortality data
7. Create the data for the sdm (require the raw climate data) using **sdmData.R**

## Optional programs (must be run after mandatory)

1. Describe the $ G $ and $ \mu $ databases. Used to create the tables in supplementary information using **describeDatabase.R**
4. To read Purves 2007's file automatically. This program can be used to add new species (as long as they are in Purves 2007's file) by modifying the data table tsn_dt ligne 47. The program **readS3_purves2007.R** was used to create the purves2007_allometries.rds
5. To create maps: **mapData.R**, and **mapDataConvexHull.R**

## Remarks

- The script **parametersAllometries.R** is called by other scripts. This program was done manually, and then automatised in **readS3_purves2007.R**
