print("Starting glmm_fixef.R")

# Species specific linear response to climate
formula1 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	annual_mean_temperature + annual_precipitation"

# Species specific quadratic response to climate
formula2 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	annual_mean_temperature + I(annual_mean_temperature^2) +
	annual_precipitation + I(annual_precipitation^2)"

# No climate response, growth change only with canopy status (correspond to Strigul 2008)
formula3 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	canopy_status"

# Sp-specific climate and canopy, climate effect might be stronger/weaker on canopy trees:
# 		- stronger negative temperature effect on canopy trees
# 		- stronger negative drougth effect on understorey trees (weak root system)
formula4 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	canopy_status * (annual_mean_temperature +
		I(annual_mean_temperature^2) +
		annual_precipitation + I(annual_precipitation^2))"

# Sp-specific linear response to individual size only
formula5 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) + dbh"

# Sp-specific quadratic response to individual size only
formula6 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) + dbh + I(dbh^2)"

# Sp-specific quadratic response to individual size, climate, and canopy status
formula7 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) + canopy_status +
	annual_mean_temperature + I(annual_mean_temperature^2) +
	annual_precipitation + I(annual_precipitation^2) +
	dbh + I(dbh^2)"

# Beyond Optimal Model without dbh-climate interaction; clim-canopy is interacting
formula8 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	canopy_status*(annual_mean_temperature +
		I(annual_mean_temperature^2) +
		annual_precipitation + I(annual_precipitation^2)) +
	dbh + I(dbh^2)"

# Beyond Optimal Model, only dbh-PP interaction; clim-canopy is interacting
formula9 = "growth ~ 1 + (1 | plot_id) +
	(1 | year_measured) + (1 | plot_id:year_measured) +
	canopy_status*(annual_mean_temperature +
		I(annual_mean_temperature^2) +
		annual_precipitation + I(annual_precipitation^2)) +
	(annual_precipitation + I(annual_precipitation^2))*(dbh + I(dbh^2))"

# For formula 9, bigger trees required more water and are more exposed to sun

# List all formulae
formulae_dt = data.table(name = paste0("model", 1:9, "_fixef"),
	formulae = c(formula1, formula2, formula3, formula4, formula5, formula6,
		formula7, formula8, formula9))

## Run, /!\ REML = FALSE /!\, third step of Zuur 2009
(current_modname = formulae_dt[formula_id, name])
model = lmer(formula = as.formula(formulae_dt[formula_id, formulae]), data = growth_data_norm,
	REML = FALSE)

saveRDS(model, paste0(savePath, current_modname, ".rds"))

r_squared = r.squaredGLMM(model)
print(r_squared)

## Variance inflation factor (VIF) among climatic variables
varInfFactor = car::vif(model)
saveRDS(varInfFactor, paste0(savePath, current_modname, "_vif.rds"))
