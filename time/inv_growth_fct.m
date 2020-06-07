function [ inv_growth ] = inv_growth_fct( x, params_g, scaling_g, dbh_scaling_g )
%% Variable access
% --- Growth
g0 = params_g.beta0;
g1 = params_g.beta1;
g2 = params_g.beta2;

% --- Scaling
mu_g = scaling_g.mu;
sd_g = scaling_g.sd;
mu_dbh_g = dbh_scaling_g.mu;
sd_dbh_g = dbh_scaling_g.sd;

%% Inverse growth function (inverse, not reciprocal!)
inv_growth = 1/exp(sd_g * (g0 + ...
    g1*(x - mu_dbh_g)/sd_dbh_g + ...
    g2*((x - mu_dbh_g)/sd_dbh_g)^2) + mu_g);

end
