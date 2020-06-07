%% Clear console and variables
% clear all
% clc

disp('Matlab main.m is starting')

%% Parallel pool cluster
% create a local cluster object
pc = parcluster('local')

% explicitly set the JobStorageLocation to the temp directory that was
% created in the sbatch script submission.sh
pc.JobStorageLocation = strcat('/scratch/amael/', getenv('SLURM_JOB_ID'))

array_id = str2num(getenv('SLURM_ARRAY_TASK_ID'));
disp(['array id: ', num2str(array_id)])

slurm_cpus = str2num(getenv('SLURM_CPUS_PER_TASK'));
disp(['number of CPUs: ', num2str(slurm_cpus)])

parpool(pc, str2num(getenv('SLURM_CPUS_PER_TASK')))

%% Load data
% --- Species names
ls_species = readtable('../createMatlabData/ls_species.csv');

% --- Species-specific integral bounds (dbh corresponding to 45m height)
integral_bounds = readtable('../createMatlabData/dbh_params.csv');
integral_bounds.Properties.RowNames = integral_bounds.species_id;

%% Run
% --- Species-specific parameters and data
currentSpecies = ls_species.x{array_id};
disp(['species id: ', currentSpecies])

scalingGrowth = readtable(char(strcat('../createMatlabData/growthScaling.csv')));
scalingGrowth.Properties.RowNames = scalingGrowth.species_id;
scalingGrowth = scalingGrowth(currentSpecies, {'mu', 'sd'});

dbh_scalingGrowth = readtable(char(strcat('../createMatlabData/growthDbhScaling.csv')));
dbh_scalingGrowth.Properties.RowNames = dbh_scalingGrowth.species_id;
dbh_scalingGrowth = dbh_scalingGrowth(currentSpecies, {'mu', 'sd'});

climate_under_g = readtable(char(strcat('../R0/Matlab_data/', currentSpecies, '/matlabGrowth_below.csv')));
climate_under_g = climate_under_g(:, {'beta0', 'beta1', 'beta2'});

s_inf = integral_bounds(currentSpecies, 'dbh_infinity').dbh_infinity;

n = height(climate_under_g);
fec = 0.0071;

% --- Create folder to store the results
if ~exist('./results', 'dir')
	mkdir('./results')
end

if ~exist(strcat('./results/', currentSpecies), 'dir')
	mkdir(strcat('./results/', currentSpecies))
end

% --- Define vector of results to save
time_integ = zeros(n, 1);

disp('parfor loop starting')

tic
parfor (j = 1:n)
	fixef_growth_under = climate_under_g(j, :);
	time_integ(j) = integral( @(x) inv_growth_fct( x, fixef_growth_under, scalingGrowth, dbh_scalingGrowth ), 0, s_inf, 'ArrayValued', true);
end
toc
csvwrite(char(strcat('./results/', currentSpecies, '/time_integ.csv')), time_integ)

delete(gcp);
exit;

% https://www.mathworks.com/help/matlab/ref/rowfun.html
% climate.TestAvg = mean(climate{:,4:end}, 2);
