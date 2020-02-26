function [ crownArea ] = dbhToCrownArea( dbh, phi_star, a, b, T_param, C0_C1, mm )
% Due to precision, I get sometimes negative values close to zero
if dbh < 0
	dbh = 0;
elseif phi_star < 0
	phi_star = 0;
end

% Table S2
R0_C0 = C0_C1('R0', 'C0').C0;
R0_C1 = C0_C1('R0', 'C1').C1;

R40_C0 = C0_C1('R40', 'C0').C0;
R40_C1 = C0_C1('R40', 'C1').C1;

M_C0 = C0_C1('M', 'C0').C0;
M_C1 = C0_C1('M', 'C1').C1;

B_C0 = C0_C1('B', 'C0').C0;
B_C1 = C0_C1('B', 'C1').C1;

% Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
R0 = (1 - T_param)*R0_C0 + T_param*R0_C1;
R40 = (1 - T_param)*R40_C0 + T_param*R40_C1;

M = (1 - T_param)*M_C0 + T_param*M_C1;
B = (1 - T_param)*B_C0 + T_param*B_C1;

% Calculate potential max radius Eq S1.6, /!\ dbh in cm /!\
Rp_max = R0 + (R40 - R0)*dbh/40;
if mm
	Rp_max = R0 + (R40 - R0)*dbh/400;
end

% All the following part is vectorised
% Convert dbh to height. /!\ If dbh is in mm, then height is in m /!\
dbhToHeight = @(dbh_loc, a_loc, b_loc) 10^(a_loc - b_loc + b_loc*log10(dbh_loc));

height = dbhToHeight(dbh, a, b);
height_star = dbhToHeight(phi_star, a, b);

% Calculate distance to the top
distToTop = height - height_star;

% Vector crown radius
crownRadius = Rp_max * ...
	( min(distToTop, height*M) / (height*M) )^B; % R_{i, y}^p

crownArea = pi * crownRadius * crownRadius;

end

