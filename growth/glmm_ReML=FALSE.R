
#### Aims of prog:
## Select the best climate variables:
# 		- Mean temperature growing season (K) and precip per year (PP)
# 		- Min temperature year (k) and precip per year (PP)
# 		- Mean temperature growing season (K) and precip per year (pp)
# 		- Min temperature year (k) and precip growing season (pp)
#
# Diverse remarks:
# Rk1: I expect mean temperature growing season (K) to be in the model.
# Rk2: For precip, maybe per year (PP) could compensate for drier growing season (pp)

print("starting glmm_ReML=FALSE.R")

## K & PP (second expected climate variables to be selected)
formula1 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(annual_mean_temperature + I(annual_mean_temperature^2) +
		mean_temperature_of_wettest_quarter + I(mean_temperature_of_wettest_quarter^2) +
		annual_precipitation + I(annual_precipitation^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2))"

## k & PP
formula2 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(annual_mean_temperature + I(annual_mean_temperature^2) +
		annual_precipitation + I(annual_precipitation^2))"

## K & pp (third expected climate variables to be chosen)
formula3 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(annual_mean_temperature + I(annual_mean_temperature^2) +
		mean_temperature_of_wettest_quarter + I(mean_temperature_of_wettest_quarter^2) +
		annual_precipitation + I(annual_precipitation^2) +
		precipitation_of_coldest_quarter + I(precipitation_of_coldest_quarter^2))"

## k & pp
formula4 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(annual_mean_temperature + I(annual_mean_temperature^2) +
		mean_temperature_of_driest_quarter + I(mean_temperature_of_driest_quarter^2) +
		annual_precipitation + I(annual_precipitation^2) +
		precipitation_of_coldest_quarter + I(precipitation_of_coldest_quarter^2))"

formula5 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_driest_quarter +
		I(mean_temperature_of_driest_quarter^2) +
		mean_temperature_of_wettest_quarter + I(mean_temperature_of_wettest_quarter^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2) +
		precipitation_of_coldest_quarter + I(precipitation_of_coldest_quarter^2))"

formula6 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_driest_quarter +
		I(mean_temperature_of_driest_quarter^2) +
		annual_precipitation + I(annual_precipitation^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2) +
		precipitation_of_coldest_quarter + I(precipitation_of_coldest_quarter^2))"

## All variables (first expected model to be selected, but interpretability + variance inflation factor?)
formula7 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(annual_mean_temperature + I(annual_mean_temperature^2) +
		temperature_annual_range + I(temperature_annual_range^2) +
		mean_temperature_of_wettest_quarter + I(mean_temperature_of_wettest_quarter^2) +
		mean_temperature_of_driest_quarter + I(mean_temperature_of_driest_quarter^2) +
		mean_temperature_of_warmest_quarter + I(mean_temperature_of_warmest_quarter^2) +
		mean_temperature_of_coldest_quarter + I(mean_temperature_of_coldest_quarter^2) +
		annual_precipitation + I(annual_precipitation^2) +
		precipitation_of_wettest_quarter + I(precipitation_of_wettest_quarter^2) +
		precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2) +
		precipitation_of_coldest_quarter + I(precipitation_of_coldest_quarter^2))"

## Formula8
formula8 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_warmest_quarter +
		I(mean_temperature_of_warmest_quarter^2) +
		precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2))"

##
formula9 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_coldest_quarter +
		I(mean_temperature_of_coldest_quarter^2) +
		precipitation_of_wettest_quarter + I(precipitation_of_wettest_quarter^2))"

##
formula10 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_warmest_quarter +
		I(mean_temperature_of_warmest_quarter^2) +
		mean_temperature_of_coldest_quarter + I(mean_temperature_of_coldest_quarter^2) +
		precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2) +
		precipitation_of_wettest_quarter + I(precipitation_of_wettest_quarter^2))"

##
formula11 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_warmest_quarter +
		I(mean_temperature_of_warmest_quarter^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2))"

##
formula12 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_wettest_quarter +
		I(mean_temperature_of_wettest_quarter^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2))"

##
formula13 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_warmest_quarter +
		I(mean_temperature_of_warmest_quarter^2) +
		precipitation_of_wettest_quarter + I(precipitation_of_wettest_quarter^2))"

##
formula14 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_warmest_quarter +
		I(mean_temperature_of_warmest_quarter^2) +
		mean_diurnal_range + I(mean_diurnal_range^2) +
		precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2))"

##
formula15 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_wettest_quarter +
		I(mean_temperature_of_wettest_quarter^2))*
		(precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2))"

##
formula16 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	(canopy_status + dbh + I(dbh^2))*(mean_temperature_of_warmest_quarter +
		I(mean_temperature_of_warmest_quarter^2))*
		(precipitation_of_wettest_quarter + I(precipitation_of_wettest_quarter^2))"

## List all formulae
formulae_dt = data.table(name = paste0("model", 1:(limReML - limBeyond_OM), "_ReML=FALSE"),
	formulae = c(formula1, formula2, formula3, formula4, formula5, formula6, formula7,formula8,
		formula9, formula10, formula11, formula12, formula13, formula14, formula15, formula16))

## Run, /!\ REML = FALSE /!\, third step of Zuur 2009
current_modname = formulae_dt[formula_id, name]
model = lmer(formula = as.formula(formulae_dt[formula_id, formulae]), data = growth_data_norm,
	REML = FALSE)

saveRDS(model, paste0(savePath, current_modname, ".rds"))

r_squared = r.squaredGLMM(model)
print(r_squared)

## Variance inflation factor (VIF) among climatic variables
varInfFactor = car::vif(model)
saveRDS(varInfFactor, paste0(savePath, current_modname, "_vif.rds"))
