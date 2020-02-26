# Execution order of the programs

## Mandatory programs

To select the best model for growth, we follow Zuur et. al (2009, section 5.7 p. 121-122).

- Create the beyond optimal model (OM), ReML = TRUE, to find best random structure
- Then search for best fixed structure, ReML = FALSE
- Rerun the selected model with ReML = TRUE

1. To select the best model and reproduce the results, use the program **growth_selection.R**. This program is designed to run on the supercomputer mammouth, which use [SLURM workload manager](https://slurm.schedmd.com/documentation.html "SLURM doc"). This program calls:
   - **glmm_beyondOM.R**, to find the beyond optimal model (*i.e.,* find the random structure)
   - **glmm_ReML=FALSE.R**, to find the best climatic variables
   - **glmm_fixef.R**, to find the best model and validate the importance of climate, competition, and dbh
2. **finalRun.R**
3. Use the program **mega_rsq_AICc.R** to calculate the AICc for the 14 models (5 from glmm_ReML=FALSE.R, and 9 from glmm_fixef.R). The tex file generated is almost the one used in the article.
4. Use **mega_vif.R** to calculate the variance inflation factor
5. **final_rsq_AICc.R**
6. **AICc_analyses.R**
7. **figure_deltaAICc.R** to make the figure used in the article to select the growth model
8. **plotGrowthDBH.R** to show the extrapolation for the juveniles

## Optional programs (must be run after mandatory)

1. **residuals.R**, to check the residuals of the selected model
