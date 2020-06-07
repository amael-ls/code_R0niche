# Execution order of the programs

Due to the stiffness of the function $\exp \left[-\int \frac{1}{G}\right]$ , I used Matlab instead of R for this part. Matlab is a non-open software, but is available in most universities.

## Mandatory programs

1. **main.m** to calculate the time it takes to grow up to a given height
2. **calc_time.R** to extract min, max, mean, and median time with the 2006-2010 climate data

## Dependent scripts

All the following files are **functions** called by main.m:

- **growthFct.m**, the individual growth function
