function [ output_args ] = integrand( x, s_star, params_g, params_m, scaling_g, ...
	dbh_scaling_g, dbh_scaling_m, allometries, C0_C1, mm )
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

% --- Allometries
a = allometries.a;
b = allometries.b;
T = allometries.T;

%% Growth and mortality functions
growth = @(dbh, c0, c1, c2, scmu_g, scsd_g, scmu_dbh, scsd_dbh) ...
	exp(scsd_g * (c0 + ...
    c1*(dbh - scmu_dbh)/scsd_dbh + ...
    c2*((dbh - scmu_dbh)/scsd_dbh)^2) + scmu_g);

mortality = @(dbh, c0, c1, c2, scmu_dbh, scsd_dbh) 1 ./ ...
    (1 + exp(-(c0 + c1*(dbh - scmu_dbh)/scsd_dbh + c2*((dbh - scmu_dbh)/scsd_dbh)^2)));

%% Result of the integrand 'expected reproduction of an individual at fixed size x'
output_args = dbhToCrownArea( x, s_star, a, b, T, C0_C1, mm ) .* ...
	1./growth(x, g0, g1, g2, mu_g, sd_g, mu_dbh_g, sd_dbh_g) .* ...
    exp(-integral(@(y) mortality(y, m0, m1, m2, mu_dbh_m, sd_dbh_m) ./ ...
    growth(y, g0, g1, g2, mu_g, sd_g, mu_dbh_g, sd_dbh_g), ...
	s_star, x, 'ArrayValued', true));

end
