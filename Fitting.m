%% Fitting.m - Main System Identification Script
% This script loads experimental plasma data, configures a physical grey-box 
% chamber model, and performs iterative parameter estimation (system identification).
% It outputs the fitted model, visualizes the performance against validation data,
% and saves the results.

%% Start-up
clear all; close all;

% Add all subfolders to the path (to include all helpers and classes). https://g.co/gemini/share/ecd4827b3d77
current_folder = fileparts(which(mfilename));
addpath(genpath(current_folder));

rmpath(genpath("FitResults"))

% Start/wipe diary
diary off
delete("Temp\temp_diary.txt")
diary("Temp\temp_diary.txt")

%% Settings
% ---- Data source ----
% Define a cell array of files to stitch together and define time windows for each file.
settings.data_filenames = {"..\Data\47118measuredTest.mat"};
settings.start_time = [0.44];
settings.end_time   = [0.89];

%settings.data_filenames = {"..\Data\47119measuredTest.mat"};
%settings.start_time = [0.44];
%settings.end_time   = [0.8];

%settings.data_filenames = {"..\Data\47118measuredTest.mat","..\Data\47119measuredTest.mat"};
%settings.start_time = [0.44, 0.44];
%settings.end_time   = [0.89, 0.80];

%settings.data_filenames = {"..\Data\50654measuredTest.mat"};
%settings.start_time = [0.4];
%settings.end_time   = [0.85];

%settings.data_filenames = {"..\Data\50654measuredTest.mat", "..\Data\50648measuredTest.mat", "..\Data\50643measuredTest.mat"};
%settings.start_time = [0.40, 0.35, 0.35];%t=0.305
%settings.end_time   = [0.85, 0.85, 0.9];

%settings.data_filenames = {"..\Data\TCV69147.mat"};
%settings.start_time = [1.15];
%settings.end_time   = [1.68];


% --- Model Selection ---
% Available models: FourChamberModelUniversal(), ThreeChamberModelUniversal(), TwoChamberModelIons(), TwoChamberModelDivertor()
settings.chamber_model = FourChamberModelUniversal();

% --- Data Processing Settings ---
settings.save_data = true; % Set to 'true' to save all results
settings.save_data_folder_prefix = '1061aRobust'; % Prefix for results folder
settings.inspect_data = 3; % 0:None, 1:Time domain, 2: frequency response, 3:Corr/SPA, 4:Full LPM, 5:Hankel
                           
settings.validation_split = 1 - 1; % Fraction of data used for VALIDATION.

% IMPORTANT: For multi-experiment data, auto-normalization (1 or 2) can be dangerous
% if the physical scaling differs between shots. It is safer to use 0 (No norm) or provide a fixed Matrix [offset, scale].
settings.normalize_data = 1; % Normalization method:
                             % 0/false: No normalization
                             % 1/true: Scale data so the maximum is 1
                             % 2: Scale data so min is 0 and max is 1

settings.detrend_data.setting = true; % 'true' to enable data detrending
settings.filter_data = true; % 'true' to enable data filtering
settings.smooth_data = 0; % Smoothing window (seconds)
settings.sample_time = 2.5e-3;% The new sample time for the data. In general should be the worst of all you diagnostics. For MSTU that is 2.5e-3, for TCV that is 1.25e-3

% --- Fitting settings ---
settings.verbose = 4;%2: MSE/Loss plot, 3: +Parameter evolution plot, 4: +Fitting time plot

settings.search_method = 'auto'; % Optimization algorithm
settings.step_size = 5;
settings.output_weight = settings.chamber_model.output_weight;%diag([1, 1, 5, 1]);% How much weight each output has. Can be "noise", can also be diag([1, 1, 5, 1]); (third channel has 20 times more weight). Or an empty array [] if you want it all equally.

settings.batch_cycles = 30;% Number of batch/coordinate fiting/decent.
settings.max_iterations = 1*10^2;% Maximum number of itterations per cycle

% Robustness Analysis Settings
settings.robustness.num_runs = 16;% Set to 1 if you don't want to do robustness analysis.
settings.robustness.perturbation = 0.2;

settings % Display settings

%% Getting data and pre-proccecing
% ---- Settings and initial guesses of chamber model ----
settings.parameters = settings.chamber_model.default_parameters;
settings.margin_factors = settings.chamber_model.default_margin_factors;
settings.fcn_type = 'c';

% Load, preprocess, merge and split data using helper function
[data, training_data, validation_data, shot, settings] = loadAndPreprocessData(settings);

% --- Define Batch Strategy ---
ode_function = @(varargin) settings.chamber_model.getMatrices(varargin{:});
sys_temp = idgrey(ode_function, settings.parameters, settings.fcn_type);
sys_temp = setParameterBounds(sys_temp, settings.parameters, settings.margin_factors);
    
    % Get all parameter names from the model
    all_param_names = {sys_temp.Structure.Parameters.Name};
    
    % Batch 1: Fit only the 'C' (output scaling) parameters.
    batch1.data = training_data;
    params_b1 = {'C1', 'C2', 'C3', 'C4'};
    batch1.params_to_fit = params_b1(ismember(params_b1, all_param_names)); 
    batch1.iterations = settings.max_iterations;
    batch1.max_cycle = 10;
    
    % Batch 2: Fit 'time constant' and 'transport' parameters.
    batch2.data = training_data;
    params_b2 = {'Ionization time', 'Confinement time', 'Recombination time', 'Div Ionization time', 'k_diff_core_div', 'k_diff_div_neu'};
    batch2.params_to_fit = params_b2(ismember(params_b2, all_param_names));
    batch2.iterations = settings.max_iterations;

    % Batch 2a: Core Physics & Core->Div Transport
    batch2a.data = training_data;
    params_b2a = {'Ionization time', 'Confinement time', 'Recombination time', 'k_diff_core_div'};% 'k_diff_div_neu'};
    batch2a.params_to_fit = params_b2a(ismember(params_b2a, all_param_names));
    batch2a.iterations = settings.max_iterations;
    
    % Batch 2b: Divertor Physics & Div->Neu Transport
    batch2b.data = training_data;
    params_b2b = {'Div Ionization time'};%, 'k_diff_div_neu'};
    batch2b.params_to_fit = params_b2b(ismember(params_b2b, all_param_names));
    batch2b.iterations = settings.max_iterations;

    % Batch 2c: Divertor Physics & Div->Neu Transport
    batch2c.data = training_data;
    params_b2c = {'Confinement time', 'Recombination time', 'k_diff_core_div'};
    batch2c.params_to_fit = params_b2c(ismember(params_b2c, all_param_names));
    batch2c.iterations = settings.max_iterations;

    % Batch 2d: Divertor Physics & Div->Neu Transport
    batch2d.data = training_data;
    params_b2d = {'Recombination time', 'k_diff_neu_ion', 'k_diff_ion_neu'};
    batch2d.params_to_fit = params_b2d(ismember(params_b2d, all_param_names));
    batch2d.iterations = settings.max_iterations;
    
    % Batch 3: Fine-tune pumping.
    batch3.data = training_data;
    params_b3 = {'Pumping time 4', 'Pumping time 1'};%, 'Leaking time'};
    batch3.params_to_fit = params_b3(ismember(params_b3, all_param_names));
    batch3.iterations = settings.max_iterations;
    
    % Batch 4: Coupling/Splitting fractions.
    batch4.data = training_data;
    params_b4 = {'Divertor fraction', 'Ionization fraction'};%, 'Recycling fraction'};
    batch4.params_to_fit = params_b4(ismember(params_b4, all_param_names)); 
    batch4.iterations = settings.max_iterations;
    
    % Batch 5: Polish All (Optional)
    batch_all_parameters.data = training_data;
    batch_all_parameters.params_to_fit = {'ALL'};
    batch_all_parameters.iterations = 1;
    
    
    %batch_strategy = {batch1, batch2, batch3, batch4};
    % batch_strategy = {batch1, batch2, batch3, batch4, batch_all_parameters};
    %batch_strategy = {batch1, batch2a, batch2b, batch3, batch4};
    batch_strategy = {batch1, batch2a, batch3, batch4};
    %batch_strategy = {batch1, batch2a, batch3};
    %batch_strategy = {batch1, batch2c, batch4};
    %batch_strategy = {batch1, batch2b, batch2d, batch3, batch4};
    % batch_strategy = {batch1, batch2a, batch2b, batch3, batch4, batch_all_parameters};

%% Inspect the data
% Loop through each experiment individually for inspection
num_experiments = size(data, 4);
for k = 1:num_experiments
    fprintf('\n--- Inspecting Experiment %d / %d ---\n', k, num_experiments);
    
    % Extract the single experiment data
    exp_data = getexp(data, k);
    
    % Create a temporary settings struct for this specific iteration
    temp_settings = settings;
    
    % Handle array vs scalar settings safely
    if length(settings.start_time) >= k
        temp_settings.start_time = settings.start_time(k);
    elseif ~isempty(settings.start_time)
         temp_settings.start_time = settings.start_time(1); % Fallback to scalar
    end
    
    if length(settings.end_time) >= k
        temp_settings.end_time = settings.end_time(k);
    elseif ~isempty(settings.end_time)
         temp_settings.end_time = settings.end_time(1); % Fallback to scalar
    end

    if isfield(settings, 'LPM_data') && length(settings.LPM_data) >= k
        temp_settings.LPM_data = settings.LPM_data(k);
    end

    % Run inspection on the single experiment
    inspectData(exp_data, temp_settings);
end
inspectFrequencyResponse(data, settings.LPM_data, [], "Combined");

%% Estimate delay
disp('Estimating initial delay for each input-output pair...');
%num_outputs = size(data.y, 2);
%num_inputs = size(data.u, 2);

num_experiments = size(data, 4);
num_outputs = size(data, 2); % For multi id data
num_inputs = size(data, 3);  % For multi id data
initial_delays = zeros(num_outputs, num_inputs); % This will store delays in samples

for i = 1:num_outputs % Loop through all combinations of inputs and outputs
    for j = 1:num_inputs
        for k = 1:num_experiments
            % Create a SISO iddata object for the current I/O pair
            siso_data = getexp(data(:, i, j), k);% Do it for each experiment
            
            initial_delays(i, j) = delayest(siso_data);
            
            % Display the result in seconds for clarity
            delay_in_seconds = initial_delays(i, j) * siso_data.Ts;
            disp(sprintf('Initial delay from input ''%s'' to output ''%s'': %.4f seconds (%d samples).', data.InputName{j}, data.OutputName{i}, delay_in_seconds, initial_delays(i, j)));
        end
    end
end

%keyboard
%% Use idgrey
tic
disp(' ')
disp('Fitting grey-box system')

% Determine number of loops based on settings
if settings.robustness.num_runs>1
    total_runs = settings.robustness.num_runs;
    fprintf('Robustness analysis ENABLED. Running %d fits with random initial perturbations.\n', total_runs);
else
    total_runs = 1;
end

% Initialize storage for the ensemble of models
dummy_struct = struct('sys', [], 'history', [], 'mse', [], 'initial_params', []);
all_fits = repmat(dummy_struct, 1, total_runs);

% Make parallel pool, needs toolbox. You can also disable this and use a normal for loop.
current_pool = gcp('nocreate');
if isempty(current_pool)
    disp("Starting parallel proccecing")
    parpool; % Start default pool
end

% --- START PARALLEL LOOP ---
%parfor% Use when many runs. Since this enables runs to be done in parallel.
%for% Use for few runs since this enables the actuall fitting to be multi-core. Or when you don't have the parallel proccecing package.
parfor run_idx = 1:total_runs%Using parfor to multithread
    
    % Note: fprintf inside parfor may not appear in order in the command window
    % but we print the start to indicate activity.
    fprintf('Starting Fit Run %d / %d on worker...\n', run_idx, total_runs);
    
    % Create a local copy of settings for this worker to modify safely
    local_settings = settings;
    
    % Reduce verbosity for workers to prevent console spamming
    % (Only keep it high if you really need to debug a specific worker crash)
    local_settings.verbose = settings.verbose;
    %local_settings.verbose = 0;
    
    % --- 1. Prepare Parameters for this run ---    
    % Logic: New = Old * (1 + perturbation * (random number between -1 and 1))
    current_parameters = settings.parameters; % Start with nominal
    
    if run_idx > 1
        n_params = size(current_parameters, 1);
        for p = 1:n_params
            val = current_parameters{p, 2};
            % Only perturb if the value is non-zero
            if val ~= 0
                % rand() is safe in parfor (streams are independent)
                noise = (rand() * 2 - 1) * local_settings.robustness.perturbation; 
                current_parameters{p, 2} = val * (1 + noise);
            end
        end
    end
    disp(current_parameters)
    
    % --- 2. Initialize System ---
    % We must redefine the function handle locally if it relies on workspace variables,
    % though usually passing 'local_settings' is enough.
    ode_function = @(varargin) local_settings.chamber_model.getMatrices(varargin{:});
    
    sys_init = idgrey(ode_function, current_parameters, local_settings.fcn_type);
    
    % Apply bounds based on the *current* randomized parameters
    sys_init = setParameterBounds(sys_init, current_parameters, local_settings.margin_factors);

    % --- 3. Perform Fitting ---
    % Use the local settings and sys_init
    [sys_run, history_run] = fitGreyboxInBatches(sys_init, batch_strategy, local_settings, validation_data);
    sys_run.Name = local_settings.chamber_model.name;
    
    % --- 4. Pack Results for Output ---
    % In parfor, we must construct the struct entirely, then assign it to the slice.
    fit_result = struct();
    fit_result.sys = sys_run;
    fit_result.history = history_run;
    fit_result.mse = sys_run.Report.Fit.MSE;
    fit_result.initial_params = current_parameters;
    
    all_fits(run_idx) = fit_result;
end

% --- Shutdown Parallel Pool ---
delete(gcp('nocreate')); % Close the pool to free resources

% --- 5. Select Best Fit ---
% Find the run with the lowest MSE to use for the rest of the script
[min_mse, best_idx] = min([all_fits.mse]);

fprintf('\n=== Robustness Analysis Complete ===\n');
fprintf('Best Fit: Run %d with MSE: %f\n', best_idx, min_mse);
fprintf('Worst Fit: Run %d with MSE: %f\n', find([all_fits.mse] == max([all_fits.mse]), 1), max([all_fits.mse]));

% Assign the first result to the variables expected by the rest of the script
sys_grey = all_fits(1).sys;
history = all_fits(1).history;

% --- Plot History and print fit ---
% Reconstruct the initial system used for this specific run (for plots)
ode_function = @(varargin) settings.chamber_model.getMatrices(varargin{:});
sys_init_plot = idgrey(ode_function, all_fits(1).initial_params, settings.fcn_type);

% Call the plotting function
plotBatchHistory(history, sys_init_plot, settings.verbose);

% Restore nominal parameters to settings so we don't save the randomized ones as "settings"
%settings.parameters = nominal_parameters; 
fit_report = sys_grey.Report;
temp = history(end);
temp = temp{1,1}.history;
printGreyFitInfo(sys_grey, 3, temp.up_to_date_pvec_sd)% Print a summary of the fit (parameter values, uncertainties)

% ---- Inspect the estimation ----
toc

if settings.robustness.num_runs>1
    %fprintf('\n--- Analyzing Robustness ---\n');

    analyzeRobustness(all_fits, sys_grey, batch_strategy, settings);
end

% --- Inspect Estimation Results (Per Experiment) ---
inspectFrequencyResponse(data, [], sys_grey, "Model only");

% Loop through each experiment to generate inspection plots. If you want combined that can be done in ModelTesting.m
num_experiments = size(data, 4);

for k = 1:num_experiments
    fprintf('\n--- Inspecting Estimation Result for Experiment %d / %d ---\n', k, num_experiments);
    
    % Extract single experiment data
    exp_data = getexp(data, k);
    
    % Prepare settings with scalar start/end times for this specific experiment
    % This ensures functions like frf_cutsignal receive the correct scalar time t1/t2
    temp_settings = settings;
    if length(settings.start_time) >= k
        temp_settings.start_time = settings.start_time(k);
    elseif ~isempty(settings.start_time)
        temp_settings.start_time = settings.start_time(1); % Fallback
    end
    
    if length(settings.end_time) >= k
        temp_settings.end_time = settings.end_time(k);
    elseif ~isempty(settings.end_time)
         temp_settings.end_time = settings.end_time(1); % Fallback
    end
    
    % Create a descriptive title including Experiment Number
    fit_name = append("Grey-box (Exp ", string(k), "), MSE=", string(fit_report.Fit.MSE));
    
    % Run inspection on the single experiment
    inspectEstimation(sys_grey, exp_data, temp_settings, fit_name);
end

%% Save data
% --- Initialize data saving ---
tic
disp("Saving data")
% Use the first filename in the list as the base for the output filename
[filepath, filename, fileextension] = fileparts(settings.data_filenames{1});

timestamp = datestr(now, 'yyyy-mm-dd_HH.MM'); % Create timestamp string

% Create a scalar string for the shot identifier to prevent directory_name becoming an array
if length(shot) > 1
    shot_str = string(shot(1)) + "_plus_" + string(length(shot)-1);
else
    shot_str = string(shot);
end

% Construct the folder name using the scalar shot_str and settings.normalize_data
folder_name = append(settings.save_data_folder_prefix, '_', timestamp, '_', ...
    settings.chamber_model.name, '_shot_', shot_str, ...
    '_norm_', string(settings.normalize_data_original_setting), ...
    '_valsplit_', string(settings.validation_split), ...
    '_detrend_', string(settings.detrend_data.setting), ...
    '_batches_', string(settings.batch_cycles), '\');
    %'_maxiter_', sprintf('%.1e', settings.max_iterations), '\');

if settings.save_data
    directory_name = append('..\FitResults\', folder_name);
else
    directory_name = append('..\FitResults\NonRelevant\', folder_name);
end

mkdir(directory_name); % Create the directories to save data/plots in.
mkdir(append(directory_name, '\PDF'));
mkdir(append(directory_name, '\FIG'));
mkdir(append(directory_name, '\SVG'));

% --- Copy files used for fitting --- 
% Save a snapshot of the *exact* models/code that was run.
copyfile('fitting.m', append(directory_name, 'fitting.m'))
copyfile('ChamberModelsClasses', append(directory_name, 'ChamberModelsClasses'))
%copyfile('Helpers', append(directory_name, 'Helpers'))

% --- Save workspace ---
save(append(directory_name, filename, ".mat"))

% --- Adjust figures and save ---
figHandles = findall(groot, 'Type', 'figure'); % Find all currently open figure windows
for i = 1:length(figHandles)
    fig = figHandles(i);
    
    styleAndResizeFigure(fig); % Change figures to white mode

    if ~isempty(fig.Name)
        % Use the figure's name. We append the figure number to ensure
        % uniqueness, especially for plots like 'Bode' that are
        % generated in a loop with the same name.
        figname = sprintf('figure_%d_%s', fig.Number, fig.Name);
    else
        % Fallback for any figures that might not have a name set
        figname = sprintf('figure_%d', fig.Number);
    end

    %exportgraphics(fig, append(directory_name, figname, ".png"), 'Resolution', 600);% Save as (high resolution) PNG for easy viewing.
    print(fig, append(directory_name, figname, ".png"), '-dpng', '-r600');

    try
        exportgraphics(fig, append(directory_name, "\PDF\", figname, ".pdf"), 'ContentType', 'vector');% Save as PDF for use in report.
    catch ME
        warning('Could not save PDF for %s: %s', figname, ME.message);
    end

    try
        savefig(fig, append(directory_name, "\FIG\", figname, ".fig"));% Export as Matlab Fig for editibility.
    catch ME
        warning('Could not save FIG for %s: %s', figname, ME.message);
    end

    try
        exportgraphics(fig, append(directory_name, "\SVG\", figname, ".svg"), 'ContentType', 'vector');% Export as SVG also for editibility.
    catch ME
        warning('Could not save SVG for %s: %s', figname, ME.message);
    end
end

toc
diary off

pause(2.5)
movefile("Temp\temp_diary.txt", append(directory_name, filename, '.txt')) % Move the temporary diary log file into the final results folder so that things can be read back.

%close all;