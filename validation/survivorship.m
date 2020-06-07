function [ output_args ] = survivorship( s_star, params_g, params_m, scaling_g, ...
	dbh_scaling_g, dbh_scaling_m)
%% Variable access
% --- Growth
g0 = params_g.beta0;
g1 = params_g.beta1;
g2 = params_g.beta2;

% --- Mortality
m0 = params_m.beta0;
m1 = params_m.beta1;
m2 = params_m.beta2;

% --- Scaling
mu_g = scaling_g.mu;
sd_g = scaling_g.sd;
mu_dbh_g = dbh_scaling_g.mu;
sd_dbh_g = dbh_scaling_g.sd;

mu_dbh_m = dbh_scaling_m.mu;
sd_dbh_m = dbh_scaling_m.sd;

% Growth and mortality functions (annual mortality => offset = log(1) = 0)
growth = @(dbh, c0, c1, c2, scmu_g, scsd_g, scmu_dbh, scsd_dbh) ...
	exp(scsd_g * (c0 + ...
    c1*(dbh - scmu_dbh)/scsd_dbh + ...
    c2*((dbh - scmu_dbh)/scsd_dbh)^2) + scmu_g);

mortality = @(dbh, c0, c1, c2, scmu_dbh, scsd_dbh) 1 - ...
	exp(-exp(c0 + c1*(dbh - scmu_dbh)/scsd_dbh + c2*((dbh - scmu_dbh)/scsd_dbh)^2));

%% Result of the integrand 'expected reproduction of an individual at fixed size x'
output_args = exp(-integral(@(y) mortality(y, m0, m1, m2, mu_dbh_m, sd_dbh_m) ./ ...
    growth(y, g0, g1, g2, mu_g, sd_g, mu_dbh_g, sd_dbh_g), ...
	0, s_star, 'ArrayValued', true));

end
