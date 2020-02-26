
#### Aims of prog:
## Estimate the growth parameters using the dataset created in the folder ../createData
#		For each demographic parameter, I keep the best model using Δaic comparison
#		For growth, I expect the temperature variable "mean_temp_period_3" to provide the best fit, (TO FIND ref physio plant). For precipitation, no idea
#
## Description of variables. For climatic variables, cf
# 	1) http://cfs.nrcan.gc.ca/projects/3/8
# 	2) http://cfs.nrcan.gc.ca/projects/3/9 for more details
# "plot_id"
# "tree_id"
# "species_id"
# "growth"
# "dbh"
# "min_temp_coldest_period": The lowest temperature of any monthly minimum temperature
# "mean_temp_period_3": mean temperature for period 3 (growing season)
# "tot_pp_for_period_3": total precipitation for period 3 (growing season)
# "tot_annual_pp": The sum of all the monthly precipitation estimates
# "canopy_status": competition. If above canopy = 1; else 0 (cf create data folder and progs for more details)
#
## Diverse remarks:
# Rk1: "mean_temp_period_3", "tot_annual_pp", "tot_pp_period3", "min_temp_coldest_period"
# Since mean(x) >= min(x), I use 'capital K' for mean_temp et 'k' for min_temp.
# T (which was smart for temperature) is reserved for TRUE in R, hence I use K
# K is dedicated to (lord) Kelvin (although the temperature is in celsius)
# Same idea, precip_year >= precip_period3, so I use PP and pp respectively

#### List of optimal models
## Description of the fixed effect:
#		- Species specific intercept
#		- Canopy status is species specific and modifies intercept + slopes climate
#		- Climate effect is species specific
#		- The dbh is species specific
#		- The dbh effect interacts with climate
## Description of the group effect structure
#		- plot_id is interacting with species
#		- plot_id is interacting with species and year_measured is interacting with plot_id.
#			The second interaction is there to represent acid rain for instance (North-East)
#		- I did not tried nested group effects, as I think plot is not nested in species
#		- The three-way interaction is: plot_id group effect is species specific and change
#			given the year

print("starting glmm_beyondOM.R")

##
formula1 = "growth ~ 1 + (1 | plot_id) +
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

##
formula2 = "growth ~ 1 + (1 | plot_id) +
	(1 | plot_id:year_measured) +
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

##
formula3 = "growth ~ 1 + (1 | plot_id) +
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

##
formula4 = "growth ~ 1 + (1 | plot_id:year_measured) +
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

## List all formulae
formulae_dt = data.table(name = paste0("model", 1:4, "_OM"),
	formulae = c(formula1, formula2, formula3, formula4))

## Run
current_modname = formulae_dt[formula_id, name]
model = lmer(formula = as.formula(formulae_dt[formula_id, formulae]), data = growth_data_norm,
	REML = TRUE)

saveRDS(model, paste0(savePath, current_modname, ".rds"))

r_squared = r.squaredGLMM(model)
print(r_squared)
