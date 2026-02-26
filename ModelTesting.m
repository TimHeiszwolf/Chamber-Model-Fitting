%% ModelTesting.m - Post-Fit Analysis Script
% This script loads a previously fitted grey-box model and its settings.
% It recalculates the preprocessed data and allows for extensive post-fit 
% analysis, including time/frequency domain inspection, Monte Carlo 
% simulations, and Bode sensitivity analysis.

%% Start-up
clear all; close all;

% Add all subfolders to the path (to include all helpers and classes).
current_folder = fileparts(which(mfilename));
addpath(genpath(current_folder));
% rmpath(genpath("FitResults")) % Not needed if FitResults is not a subfolder
rmpath(genpath("ChamberModelsClasses"))
%rmpath(genpath("Helpers"))

%% Settings
perturbation_fraction = 0.5;% The pertubation factor for the bode-plot. This is how much each paramter is perturbed.

%%% --- Data/model to load ---
fit_results = '..\FitResults\1061a_2026-01-18_01.16_4 Chamber_shot_47118_plus_1_norm_1_valsplit_0_detrend_true_batches_3.0e+01\47118measuredTest.mat';%Semi old stuff.
%fit_results = '..\FitResults\1061dlong_2026-01-18_04.46_3 Chamber_shot_47118_plus_1_norm_1_valsplit_0_detrend_true_batches_1.0e+02\47118measuredTest.mat';%Still short version
%fit_results = '..\FitResults\1081a_2026-01-18_04.44_4 Chamber_shot_47118_plus_1_norm_1_valsplit_0_detrend_true_batches_30\47118measuredTest.mat';%Semi old stuff.


%fit_results = '..\FitResults\1052c_2026-01-19_07.21_3 Chamber_shot_50654_norm_1_valsplit_0_detrend_true_maxiter_1.0e+02\TCV69147.mat';
%fit_results = '..\FitResults\1034_2026-01-18_19.34_3 Chamber_shot_69147_norm_1_valsplit_0_detrend_true_maxiter_1.0e+02\50654measuredTest.mat';

% Load files
[folder_path, ~, ~] = fileparts(fit_results);
addpath(genpath(folder_path));
load(fit_results)

%%% --- Data to compare the model to/with ---
data_to_compare = "..\Data\47118measuredTest.mat"; settings.start_time=0.44;settings.end_time=0.89;%Good? maybe go to 0.48 instead of 0.44?
%data_to_compare = "..\Data\47119measuredTest.mat";settings.start_time=0.44;settings.end_time=0.8;%Good? maybe go to 0.48 instead of 0.44?% 50 hz vessel
%data_to_compare = "..\Data\50648measuredTest.mat";settings.start_time=0.4; settings.end_time=0.85;

%data_to_compare = "..\Data\TCV68861.mat";settings.start_time=1.15; settings.end_time=1.68;


% Other loading
load(data_to_compare, "data", "shot", "f_exc", "P", "f0");% The original data used for training still exists as "training_data"
close all;

sys_grey.Name = "Chamber";%Change name if desired. Is used in the plot.

settings.detrend_data.f0 = f0;
settings.detrend_data.f_exc = f_exc;
settings.detrend_data.P = P;
raw_data = data;
[data, settings] = settings.chamber_model.preProccesData(raw_data, settings);
%settings.normalize_data = settings.normalize_data_original_setting;

%% Inspect the data
inspectData(data, settings, 2);% 0:None, 1:Time domain, 2: frequency response, 3:Corr/SPA, 4:Full LPM, 5:Hankel

%% Do the comparision
%try
    %inspectEstimation(sys_ssest, data, settings, "SSEST")
    %disp("SSEST:")
    %fit_report = sys_ssest.Report;
    %printSSESTFitInfo(sys_ssest, 2)
    %printSSESTFitInfo(sys_ssest, 1)
%end

%try
    %inspectFrequencyResponse(data, settings.LPM_data, sys_grey)% TODO FIND BETTER PLACE
    inspectFrequencyResponse(data, [], sys_grey, "Model only");
    inspectEstimation(sys_grey, data, settings, "Grey")
    disp("Grey-box:")
    fit_report = sys_grey.Report;
    printGreyFitInfo(sys_grey, 3)
    
%end

%{
C = sys_grey.C;
C(2,2) = C(2,2)*2;
C(3,3) = C(3,3)*5;
sys_ss = idss(sys_grey);
set(sys_ss, 'C', C);
inspectEstimation(sys_ss, data, settings, "Grey")
%}


%{
all_param_names = {sys_grey.Structure.Parameters.Name};
batch1.data = data;
params_b1 = {'C1', 'C2', 'C3', 'C4'};
batch1.params_to_fit = params_b1(ismember(params_b1, all_param_names)); % Only use params that exist
batch1.iterations = settings.max_iterations;
batch1.max_cycle = 50;
batch_strategy = {batch1};

settings.verbose = 3;
settings.batch_cycles = 20;        
settings.max_iterations = 1*10^2

sys_init = idgrey(sys_grey)
[sys_grey_test, history] = fitGreyboxInBatches(sys_init, batch_strategy, settings);
inspectEstimation(sys_grey_test, data, settings, "Grey")
%}


%inspectEstimation(sys_ssest_free, data, settings, "SSEST fully free")

%inspectEstimation(sys_N4SID, data, settings, "N4SID")

%% --- Benchmark Models ---
% We train these on the training_data.
if exist('training_data', 'var')
    disp('--- Training Benchmark Models on original training data ---');
    
    % 1. N4SID (Subspace State-Space)
    opt_n4sid = n4sidOptions('Focus', 'simulation');
    sys_n4sid = n4sid(training_data, 4, opt_n4sid);% Order 4 to match the 4-chamber physical model structure.
    sys_n4sid.Name = 'N4SID';
    
    % 2. ARX (Linear Difference Equation)
    % Heuristic: Use memory depth of 4 (similar to order 4)
    %ny = size(training_data.OutputData, 2);
    %nu = size(training_data.InputData, 2);
    ny = size(training_data, 2); % For multi input stuff
    nu = size(training_data, 3); % For multi input stuff
    na = 4 * ones(ny, ny);
    nb = 4 * ones(ny, nu);
    nk = 1 * ones(ny, nu); % 1 sample delay
    sys_arx = arx(training_data, [na nb nk]);
    sys_arx.Name = 'ARX';
    
    % 3. Non-Linear ARX (Sigmoid Network)
    % Checks if non-linearities (saturation, recombination limits) are dominant.
    %sys_nlarx = nlarx(training_data, [na nb nk], 'sigmoidnet');
    %sys_nlarx.Name = 'NLARX';

    % 4. Transfer Function (Simple Continuous Black-box)
    % Estimating a transfer function with 4 poles to match the 4-chamber structure
    sys_tf = tfest(training_data, 4);
    sys_tf.Name = 'Transfer Fcn';
    
    % --- Compare on Test Data ---
    disp('--- Comparing Models on Test Data ---');
    figure;
    opt = compareOptions('InitialCondition', 'e');
    training_data.Name = "Data";% This ensures the legend says "Data" instead of long strings like "Validation Data: ...".
    data.Name = "Data";% This ensures the legend says "Data" instead of long strings like "Validation Data: ...".

    %compare(data, sys_grey, sys_n4sid, sys_arx, opt);
    %compare(data, sys_n4sid, sys_arx, sys_grey, opt);
    compare(data, sys_n4sid, sys_tf, sys_grey, opt);
    % Fix the look and add the line
    styleTimeDomainComparison(data, settings.validation_split, "Benchmark Models");
    
    % Quantitative Metrics (NRMSE)
    [yhat_grey, fit_grey] = compare(data, sys_grey, compareOptions);
    [yhat_n4sid, fit_n4sid] = compare(data, sys_n4sid, compareOptions);
    [yhat_tf, fit_tf] = compare(data, sys_tf, compareOptions);
    %[~, fit_arx] = compare(training_data, sys_arx, compareOptions);
    %[~, fit_nlarx] = compare(training_data, sys_nlarx, compareOptions);
    
    % Handle Multi-Experiment Results (Convert Cell Arrays to Matrix)
    if iscell(fit_grey), fit_grey = cell2mat(fit_grey); end
    if iscell(fit_n4sid), fit_n4sid = cell2mat(fit_n4sid); end
    if iscell(fit_tf), fit_tf = cell2mat(fit_tf); end
    %if iscell(fit_arx), fit_arx = cell2mat(fit_arx); end
    %if iscell(fit_nlarx), fit_nlarx = cell2mat(fit_nlarx); end

    % --- Calculate Mean Absolute Error (MAE) ---
    % Extract arrays (OutputData).
    % Handle cases where yhat is a cell array of iddata objects (multi-experiment).
    y_true_data = data.OutputData;
    
    % Helper to flatten data: If cell (multi-exp), vertcat them; else use directly.
    if iscell(y_true_data)
        y_true_data = vertcat(y_true_data{:});
    end

    if iscell(yhat_grey)
        % yhat_grey is a cell array of iddata objects -> Extract OutputData from each, then stack.
        temp_grey = cellfun(@(x) x.OutputData, yhat_grey, 'UniformOutput', false);
        y_grey_data = vertcat(temp_grey{:});
    else
        y_grey_data = yhat_grey.OutputData;
    end

    if iscell(yhat_n4sid)
        temp_n4sid = cellfun(@(x) x.OutputData, yhat_n4sid, 'UniformOutput', false);
        y_n4sid_data = vertcat(temp_n4sid{:});
    else
        y_n4sid_data = yhat_n4sid.OutputData;
    end

    if iscell(yhat_tf)
        temp_tf = cellfun(@(x) x.OutputData, yhat_tf, 'UniformOutput', false);
        y_tf_data = vertcat(temp_tf{:});
    else
        y_tf_data = yhat_tf.OutputData;
    end

    % Calculate MAE per channel (dim 1 is time, dim 2 is channel)
    mae_grey_per_ch = mean(abs(y_true_data - y_grey_data), 1);
    mae_n4sid_per_ch = mean(abs(y_true_data - y_n4sid_data), 1);
    mae_tf_per_ch = mean(abs(y_true_data - y_tf_data), 1);
    y_smooth = smoothdata(y_true_data, 1, 'gaussian', 4); 
    mae_noise_floor = mean(abs(y_true_data - y_smooth), 1);

    fprintf('\n=== Test Set Results (Mean NRMSE) ===\n');
    fprintf('Grey-Box (Physical): %.2f%%\n', mean(fit_grey(:)));
    fprintf('N4SID (Black-box):   %.2f%%\n', mean(fit_n4sid(:)));
    fprintf('Transfer Fcn:        %.2f%%\n', mean(fit_tf(:)));
    %fprintf('ARX (Linear):        %.2f%%\n', mean(fit_arx(:)));
    %fprintf('NLARX (Non-linear):  %.2f%%\n', mean(fit_nlarx(:)));

    fprintf('\n=== Test Set Results (MAE) ===\n');
    fprintf('1. Noise Floor (Smooth): Avg: %.4f | Per Channel: [%s]\n', mean(mae_noise_floor), num2str(mae_noise_floor, '%.4f '));%Likely complete crap
    fprintf('Grey-Box (Physical): Avg: %.4f | Per Channel: [%s]\n', mean(mae_grey_per_ch), num2str(mae_grey_per_ch, '%.4f '));
    fprintf('N4SID (Black-box):   Avg: %.4f | Per Channel: [%s]\n', mean(mae_n4sid_per_ch), num2str(mae_n4sid_per_ch, '%.4f '));
    fprintf('Transfer Fcn:        Avg: %.4f | Per Channel: [%s]\n', mean(mae_tf_per_ch), num2str(mae_tf_per_ch, '%.4f '));

else
    warning('training_data not found in fit_results. Benchmarks skipped.');
end


%% Perform Sensitivity Analysis on the Grey-Box Model

if exist('sys_grey', 'var')
    disp('--- Starting Sensitivity Analysis ---')
    tic
    performSensitivityBode(sys_grey, settings, perturbation_fraction);
    toc

    tic
    %performMonteCarloAnalysis(sys_grey, data, settings, 6)
    toc
    %performSensitivityMonteCarlo(sys_grey, 100)
else
    disp('Variable "sys_grey" not found, skipping sensitivity analysis.');
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
    '_maxiter_', sprintf('%.1e', settings.max_iterations), '\');

directory_name = append('..\ModelTestingResults\', folder_name);

mkdir(directory_name); % Create the directories to save data/plots in.
mkdir(append(directory_name, '\PDF'));
mkdir(append(directory_name, '\FIG'));
mkdir(append(directory_name, '\SVG'));

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