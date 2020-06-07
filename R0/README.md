# Execution order of the programs

Due to the stiffness of the function $\exp \left[-\int \frac{\mu}{G}\right]$ , I used Matlab instead of R for this part. Matlab is a non-open software, but is available in most universities.

## Mandatory programs

1. **main.m** to calculate $R_0$ with competition 10m
2. **main_fundNiche.m** to calculate $R_0$ without competition
3. **grad_R0.R** to compute the gradients
4. **polar_grad.R** to plot the gradient in polar coordinates
5. **maps/map.R** to plot the maps
6. **corrR0_distEdge.R** to compute the correlation between $R_0$ and the closest edge (defined by Little 1971)
7. **correlTable.R** to create the .tex table file
8. **plot3correl.R** to create the plots
9. **table_R0.R** to create a table for non-scaled values of $R_{0}$

## Dependent scripts

All the following files are **functions** called by main.m:

- **dbhToCrownArea.m** to calculate the crown and height are from dbh
- **integrand.m**, the integrand of the $R_0$ equation
- **survivorship.m**, the number of individuals that survive up to the canopy

The following file is called by **corrR0_distEdge.R**

- **dist_geosphere.R**

## Optional program

- **plotCentroid.R**
- **runningTimePerSpecies.R**
