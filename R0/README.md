# Execution order of the programs

Due to the stiffness of the function $\exp \left[-\int \frac{\mu}{G}\right]$ , I used Matlab instead of R for this part. Matlab is a non-open software, but is available in most universities.

## Mandatory programs

1. **main.m** to calculate $R_0$
2. **corrR0_distEdge.R** to compute the correlation between $R_0$ and the closest edge (defined by Little 1971)
3. **correlTable.R** to create the .tex table file

## Dependent scripts

All the following files are **functions** called by main.m:

- **dbhToCrownArea.m** to calculate the crown and height are from dbh
- **integrand.m**, the integrand of the $R_0$ equation
- **survivorship.m**, the number of individuals that survive up to the canopy
