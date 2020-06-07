%% Clear console and variables
% clear all
% clc

disp('Matlab main.m is starting')

%% Parallel pool cluster
% create a local cluster object
pc = parcluster('local');

% explicitly set the JobStorageLocation to the temp directory that was
% created in the sbatch script submission.sh
pc.JobStorageLocation = strcat('/scratch/amael/', getenv('SLURM_JOB_ID'));

array_id = str2num(getenv('SLURM_ARRAY_TASK_ID'));
disp(['array id: ', num2str(array_id)])

slurm_cpus = str2num(getenv('SLURM_CPUS_PER_TASK'));

parpool(pc, str2num(getenv('SLURM_CPUS_PER_TASK')))

%% Load data
% --- Allometries
C0_C1 = readtable('../createMatlabData/C0_C1.csv');
allometries = readtable('../createMatlabData/purves2007_allometries.csv');

C0_C1.Properties.RowNames = C0_C1.parameter;
allometries.Properties.RowNames = allometries.species;

% --- Species names
ls_species = readtable('../createMatlabData/ls_species.csv');

% --- Species-specific integral bounds (dbh corresponding to 45m height)
integral_bounds = readtable('../createMatlabData/dbh_params.csv');
integral_bounds.Properties.RowNames = integral_bounds.species_id;

%% Run
% --- Species-specific parameters and data
currentSpecies = ls_species.x{array_id};
disp(['species id: ', currentSpecies])

scalingGrowth = readtable('../createMatlabData/growthScaling.csv');
scalingGrowth.Properties.RowNames = scalingGrowth.species_id;
scalingGrowth = scalingGrowth(currentSpecies, {'mu', 'sd'});

dbh_scalingGrowth = readtable('../createMatlabData/growthDbhScaling.csv');
dbh_scalingGrowth.Properties.RowNames = dbh_scalingGrowth.species_id;
dbh_scalingGrowth = dbh_scalingGrowth(currentSpecies, {'mu', 'sd'});

mu_g = scalingGrowth.mu;
sd_g = scalingGrowth.sd;
mu_dbh_g = dbh_scalingGrowth.mu;
sd_dbh_g = dbh_scalingGrowth.sd;

dbh_scalingMortality = readtable('../createMatlabData/mortalityDbhScaling.csv');
dbh_scalingMortality.Properties.RowNames = dbh_scalingMortality.species_id;
dbh_scalingMortality = dbh_scalingMortality(currentSpecies, {'mu', 'sd'});

climate_over_g = readtable(char(strcat('./Matlab_data/', currentSpecies, '/matlabGrowth_above.csv')));
climate_over_g = climate_over_g(:, {'beta0', 'beta1', 'beta2'});
climate_under_g = readtable(char(strcat('./Matlab_data/', currentSpecies, '/matlabGrowth_below.csv')));
climate_under_g = climate_under_g(:, {'beta0', 'beta1', 'beta2'});

climate_over_m = readtable(char(strcat('./Matlab_data/', currentSpecies, '/matlabMortality_above.csv')));
climate_over_m = climate_over_m(:, {'beta0', 'beta1', 'beta2'});
climate_under_m = readtable(char(strcat('./Matlab_data/', currentSpecies, '/matlabMortality_below.csv')));
climate_under_m = climate_under_m(:, {'beta0', 'beta1', 'beta2'});

s_inf = integral_bounds(currentSpecies, 'dbh_inf').dbh_inf;
age_max = integral_bounds(currentSpecies, 'age_max').age_max;
local_s = readtable(char(strcat('./results/', currentSpecies, '/local_s_inf.csv')), 'ReadVariableNames', false);
local_s = local_s.Var1;

allometries = allometries(currentSpecies, {'a', 'b', 'T'});

n = height(climate_over_g);
fec = 0.0071;

% --- Create folder to store the results
if ~exist('./results', 'dir')
	mkdir('./results')
end

if ~exist(strcat('./results/', currentSpecies), 'dir')
	mkdir(strcat('./results/', currentSpecies))
end

% --- Define vector of results to save
R0_10m = zeros(n, 1);

disp('parfor loop starting')

tic
parfor (j = 1:n)
	fixef_growth_over = climate_over_g(j, :);
	fixef_mortality_over = climate_over_m(j, :);

	fixef_growth_under = climate_under_g(j, :);
	fixef_mortality_under = climate_under_m(j, :);

	local_s_inf = local_s(j);

	current_s_inf = min(local_s_inf, s_inf);

	R0_10m(j) = pi*fec * ...
	integral( @(x) integrand( x, 0, fixef_growth_over, fixef_mortality_over, scalingGrowth, dbh_scalingGrowth, dbh_scalingMortality, allometries, C0_C1, 'true' ), 0, current_s_inf, 'ArrayValued', true);
end
toc
csvwrite(char(strcat('./results/', currentSpecies, '/R0_0m.csv')), R0_10m)

delete(gcp);
exit;

% https://www.mathworks.com/help/matlab/ref/rowfun.html
% climate.TestAvg = mean(climate{:,4:end}, 2);
