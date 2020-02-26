# Execution order of the programs

Due to the stiffness of the function $\exp \left[-\int \frac{\mu}{G}\right]$ , I used Matlab instead of R for this part. Matlab is a non-open software, but is available in most universities.

## Mandatory programs

1. **createValidationData.R**, to create the sdm data for matlab
2. **main.m** to calculate $R_0$
3. **validation.R**, to assess the correlation beween $R_0$ and the sdm, both calculated at the same location (provided by the data used in growth and mortality estimates)

## Dependent scripts

All the following files are **functions** called by main.m:

- **dbhToCrownArea.m** to calculate the crown and height are from dbh
- **integrand.m**, the integrand of the $R_0$ equation
- **survivorship.m**, the number of individuals that survive up to the canopy