# Execution order of the programs

## Mandatory programs

To select the best model for mortality and reproduce the results, use:

1. **rstanarm.R**, designed to run on the supercomputer mammouth, which use [SLURM workload manager](https://slurm.schedmd.com/documentation.html "SLURM doc").
2. **waic_mortality.R** to select the best model base on WAIC
3. **waic_analyses.R**
4. **rstanarm_submodels.R** to run alternative submodels
5. **waic_mortality_submodels.R** 
6. **waic_analyses_submodels.R**
7. **extract_fixef.R**, to extract the estimated slopes of the fixed effect
8. **rstanarm_rsq.R** to estimate the $ R^2 $ of each model. However, they cannot be compared
9. **mega_rsq_WAIC.R**, to merge $R^2$ and WAIC
10. **rebuildParameters.R**
11. **figure_deltaWAIC.R**, for the model selection figure
12. **plotMortalityDBH.R**, to show the extrapolation to the juveniles

# Scripts related to the WAIC calculus

The following scripts are from **rstanarm package**, I modified some functions to run the calculus on the hard drive rather than on the RAM. These scripts are called from other scripts, and cannot be run individually.

- as.matrix.stanreg.R
- lArgsStanReg_submodels.R
- lArgsStanReg_waic.R
- log_lik.R
- loo.R
- misc.R
- stan_glm.fit.R

# Dependent programs

The following scripts cannot be called independently.

- **selectClimaticVariables.R**
- **submodels.R**