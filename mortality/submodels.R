
# Species specific linear response to climate
if (formula_id == 1)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		min_temperature_of_coldest_month + precipitation_of_driest_quarter,
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_1.rds"))
}

# Species specific quadratic response to climate
if (formula_id == 2)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		min_temperature_of_coldest_month + I(min_temperature_of_coldest_month^2) +
			precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_2.rds"))
}

# No climate response, growth change only with canopy status (correspond to Strigul 2008)
if (formula_id == 3)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status,
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_3.rds"))
}

# Sp-specific climate and canopy, climate effect might be stronger/weaker on canopy trees:
# 		- stronger negative temperature effect on canopy trees
# 		- stronger negative drougth effect on understorey trees (weak root system)
if (formula_id == 4)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status*(min_temperature_of_coldest_month +
			I(min_temperature_of_coldest_month^2) +
			precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2)),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_4.rds"))
}

# Sp-specific linear response to individual size only
if (formula_id == 5)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		dbh,
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_5.rds"))
}

# Sp-specific quadratic response to individual size only
if (formula_id == 6)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_6.rds"))
}

# Sp-specific quadratic response to individual size, climate, and canopy status
if (formula_id == 7)
{
	model = stan_glmer(formula = deltaState ~ 1 + (1 | plot_id) +
		canopy_status + min_temperature_of_coldest_month +
			I(min_temperature_of_coldest_month^2) +
			precipitation_of_driest_quarter + I(precipitation_of_driest_quarter^2) + dbh + I(dbh^2),
		data = mortality_data_norm,
		family = binomial, iter = 2000, cores = 4, chains = 4,
		prior = normal(),
		prior_intercept = normal(),
		prior_aux = exponential(),
		prior_covariance = decov(),
		prior_PD = FALSE,
		# diagnostic_file = paste0(savePath, "diagnostic", formula_id, ".csv"),
		algorithm = "sampling")
		saveRDS(model, paste0(savePath, "submodel_7.rds"))
}
