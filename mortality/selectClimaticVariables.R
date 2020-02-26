
#### List of the formulae to select the best climatic combination
# Due to variance inflation factors (see growth estimation), I cannot combine
# all the variables. Thus, I limit the model to single temperature-precipitation
# combinations. Annual variables are particularly correlated to their derived variables
# such as x_wettest, x_driest, x_warmest, ..., where x is precipitation or temperature
#
# Due to computation time (and convergence troubles), I also use a simpler random structure.
# For growth, most of the variance was caught by plot_id, hence I use plot_id only
#
# There was a curious bug with atom (my text editor) and mac, so my way of programing from growth could not be applied here:
# For the growth, I did a data table, and then did as.formula(formulae_dt[formula_id, formulae])
# However, here, for whatever reason, it ran only the first line of the formula!
# Hence, I copy-paste my code and used if conditions (quick and dirty way, sorry!)
#
# The climate--dbh interactions were also creating a lot of convergence problems, and made some
# models probably nearly unidentifiable. Hence, I also deleted this relationship from the growth model.

## Annual variables (Lines 2010). These two variables are correlated to many others, but interpretability?
if (formula_id == 1)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(annual_mean_temperature + I(annual_mean_temperature^2) +
			annual_precipitation + I(annual_precipitation^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_1.rds"))
}

## Driest variables (water deficiency)
if (formula_id == 2)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(mean_temperature_of_driest_quarter +
			I(mean_temperature_of_driest_quarter^2) +
			precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_2.rds"))
}

## Warmest variables (hottest period, potentially combined with water deficiency)
if (formula_id == 3)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(mean_temperature_of_warmest_quarter +
			I(mean_temperature_of_warmest_quarter^2) +
			precipitation_of_warmest_quarter + I(precipitation_of_warmest_quarter^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_3.rds"))
}

## Driest month pp, driest quarter temp (extrem drought)
if (formula_id == 4)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(mean_temperature_of_driest_quarter +
			I(mean_temperature_of_driest_quarter^2) +
			precipitation_of_driest_month + I(precipitation_of_driest_month^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_4.rds"))
}

## Annual pp, driest quarter temp (in case of there is a water storage effect)
if (formula_id == 5)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(mean_temperature_of_driest_quarter +
			I(mean_temperature_of_driest_quarter^2) +
			annual_precipitation + I(annual_precipitation^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_5.rds"))
}

## Min temperature, driest precipitation month (Xylem froze + extrem drought)
if (formula_id == 6)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(min_temperature_of_coldest_month +
			I(min_temperature_of_coldest_month^2) +
			precipitation_of_driest_month + I(precipitation_of_driest_month^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_6.rds"))
}

## Min temperature, driest precipitation 3 months (Xylem froze + longer period drought)
if (formula_id == 7)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(min_temperature_of_coldest_month +
			I(min_temperature_of_coldest_month^2) +
			precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2)) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "model_7.rds"))
}

#### Bridge sampling to calculate ratios, but it does not work
## Note that it requires diagnostic_file, that I commented since bs does not work
# bs = bridge_sampler(model, cores = length(nodeslist))
# saveRDS(bs, paste0(savePath, current_modname, "_bs.rds"))
