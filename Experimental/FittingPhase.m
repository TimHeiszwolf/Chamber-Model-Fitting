%% Fitting_FrequencyDomain_PhaseOnly.m
% This script fits the parameters of the Chamber Model by minimizing the
% error between the Model Phase and the LPM Phase (Frequency Domain).
%
% It deliberately IGNORES the Magnitude response (and therefore fixes the C/Gain parameters)
% to avoid issues with signal scaling/calibration.

%% 1. Start-up
clear; close all; clc;
current_folder = fileparts(which(mfilename));
addpath(genpath(current_folder));

% Initialize Diary
diary off; if exist('Temp', 'dir'), delete("Temp\temp_diary_freq.txt"); end
if ~exist('Temp', 'dir'), mkdir('Temp'); end
diary("Temp\temp_diary_freq.txt");

%% 2. Settings
% ---- Data Files ----
settings.data_filenames = { ...
    "..\Data\47118measuredTest.mat", ...
    "..\Data\47119measuredTest.mat" ...
    };

% ---- Time Windows (per file) ----
settings.start_time = [0.44, 0.44]; 
settings.end_time   = [0.89, 0.80];

% ---- Data Processing Settings ----
settings.save_data = false; 
settings.inspect_data = 0; % We will inspect manually later
settings.validation_split = 0; % Use all data for fitting
settings.normalize_data = 0; % Normalization matters less for phase, but 0 is safest
settings.detrend_data.setting = true; 
settings.filter_data = true; 
settings.smooth_data = 0; 
settings.data_frequency = 2.5e-3; 

% ---- Model Initialization ----
settings.chamber_model = FourChamberModel(); % Or FourChamberModelTCV()
settings.parameters = settings.chamber_model.default_parameters;
settings.margin_factors = settings.chamber_model.default_margin_factors;
settings.fcn_type = 'c';

%% 3. Load Data & Compute LPM
% This function calculates the LPM FRF (including Phase) and stores it in settings.LPM_data
disp('--- Loading Data and Computing LPM ---');
[data, ~, ~, ~, settings] = loadAndPreprocessData(settings);

%% 4. Prepare Optimization
% We need to separate parameters into "Gain" (C-matrix, ignored) and "Kinetic" (A/B-matrix, fitted).
% Scalar gains (C1, C2...) do not affect phase. Fitting them to phase is impossible.

disp('--- preparing Optimization Parameters ---');

% Create dummy model to parse parameters
ode_function = @(varargin) settings.chamber_model.getMatrices(varargin{:});
sys_init = idgrey(ode_function, settings.parameters, settings.fcn_type);
sys_init = setParameterBounds(sys_init, settings.parameters, settings.margin_factors);

all_params = sys_init.Structure.Parameters;
num_params = length(all_params);
param_names = {all_params.Name};
param_values = [all_params.Value];
param_min = [all_params.Minimum];
param_max = [all_params.Maximum];

% Filter: Exclude 'C' parameters (Gains) from optimization
% Heuristic: If name starts with 'C' followed by a number, exclude it.
is_gain_param = false(1, num_params);
for i = 1:num_params
    if startsWith(param_names{i}, 'C') && ~isempty(regexp(param_names{i}, '^C\d+$', 'once'))
        is_gain_param(i) = true;
    end
end

% Indices map
idx_fit = find(~is_gain_param);
idx_fixed = find(is_gain_param);

p_guess = param_values(idx_fit);
lb = param_min(idx_fit);
ub = param_max(idx_fit);

fprintf('Optimizing the following parameters (Phase Only):\n');
disp(param_names(idx_fit)');
fprintf('Fixing the following parameters (Gain):\n');
disp(param_names(idx_fixed)');

%% 5. Define Objective Function
% The objective function calculates the weighted phase difference vector.

objective_fun = @(p_subset) costFunctionPhase(p_subset, idx_fit, param_values, settings, ode_function);

%% 6. Run Optimization
disp('--- Starting lsqnonlin Optimization ---');
options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxIterations', 50, ...
    'FunctionTolerance', 1e-6, ...
    'StepTolerance', 1e-6);

tic;
% Request 'jacobian' as the 7th output argument
[p_opt_subset, resnorm, residual, exitflag, output, lambda, jacobian] = ...
    lsqnonlin(objective_fun, p_guess, lb, ub, options);
fit_time = toc;

fprintf('\nOptimization finished in %.2f seconds.\n', fit_time);

% --- Calculate Uncertainties ---
% 1. Degrees of Freedom
num_residuals = length(residual);
num_fitted_params = length(p_opt_subset);
dof = num_residuals - num_fitted_params;

% 2. Estimate Variance of the Residuals (MSE)
% Note: If your weights were perfect (1/sigma), MSE should be close to 1.
% If not, this scales the uncertainty to match the actual fit quality.
mse = resnorm / dof;

% 3. Calculate Covariance Matrix: Cov = inv(J'*J) * MSE
% Use pseudo-inverse or backslash for numerical stability
covariance_matrix = (jacobian' * jacobian) \ eye(num_fitted_params) * mse;

% 4. Standard Deviations (square root of diagonal elements)
p_std_subset = sqrt(diag(covariance_matrix))';

%% 7. Update Model with Results
p_final = param_values;
p_final(idx_fit) = p_opt_subset;

% Update the settings struct for future use
settings.parameters(:, 2) = num2cell(p_final);

% Create the final system object
sys_final = idgrey(ode_function, settings.parameters, settings.fcn_type);
sys_final = setParameterBounds(sys_final, settings.parameters, settings.margin_factors);

%% 8. Report and Visualize

% --- Print Parameter Changes with Uncertainty ---
fprintf('\n%-20s | %-12s | %-12s | %-12s\n', 'Parameter', 'Initial', 'Final', 'Std. Dev.');
fprintf('%s\n', repmat('-', 1, 65));

% Create a full vector for SDs (initially zero)
p_std_full = zeros(size(param_values));
p_std_full(idx_fit) = p_std_subset;

for i = 1:num_params
    if ismember(i, idx_fit)
        % Fitted Parameter
        fprintf('%-20s | %-12.4e | %-12.4e | %-12.4e (*)\n', ...
            param_names{i}, param_values(i), p_final(i), p_std_full(i));
    else
        % Fixed Parameter
        fprintf('%-20s | %-12.4e | %-12.4e | %-12s\n', ...
            param_names{i}, param_values(i), p_final(i), 'Fixed');
    end
end
fprintf('(*) = Fitted parameter\n');

% --- Plotting ---
disp('--- Generating Bode Plots ---');
inspectFrequencyResponse(data, settings.LPM_data, sys_final);

diary off;

%% ---------------------------------------------------------
%  Local Functions
% ---------------------------------------------------------

function residuals = costFunctionPhase(p_subset, idx_fit, p_full_initial, settings, ode_fcn)
    % Reconstruct full parameter vector
    p_current = p_full_initial;
    p_current(idx_fit) = p_subset;
    
    % Update settings with current parameters (needed for getMatrices format)
    % Note: ode_function expects the specific parameter args, but idgrey wrapping handles
    % vectors differently. Here we manually construct the system to get Bode.
    
    % We need to convert p_current vector back to the cell array format if strictly needed,
    % but calculating matrices A,B,C,D directly is faster.
    
    % However, settings.chamber_model.getMatrices takes individual args. 
    % We can use the cell expansion of the values.
    param_cell = num2cell(p_current);
    
    try
        [A, B, C, D] = settings.chamber_model.getMatrices(param_cell{:});
        sys = ss(A, B, C, D);
    catch
        % If unstable parameters cause crash, return high residual
        residuals = 1e6 * ones(1000, 1); 
        return;
    end
    
    residuals = [];
    
    % Loop over experiments stored in settings.LPM_data
    num_experiments = length(settings.LPM_data);
    
    for k = 1:num_experiments
        lpm = settings.LPM_data(k);
        freqs = lpm.freq; % [1 x Nf]
        
        if isempty(freqs), continue; end
        
        % Calculate Model Response at these frequencies
        % bode returns mag: [Ny, Nu, Nf], phase: [Ny, Nu, Nf]
        w = 2 * pi * freqs;
        [~, phase_model_all] = bode(sys, w);
        
        [Ny, Nu, Nf] = size(phase_model_all);
        
        % Iterate over Inputs (j) and Outputs (i)
        for j = 1:Nu
            for i = 1:Ny
                % Data
                phase_data_rad = squeeze(lpm.phase(i, j, :)); % [Nf x 1] radians
                uncertainty_rad = squeeze(lpm.s2p(i, j, :));  % [Nf x 1] radians (2-sigma)
                
                % Model
                phase_model_deg = squeeze(phase_model_all(i, j, :)); % [Nf x 1] degrees
                phase_model_rad = deg2rad(phase_model_deg);
                
                % --- Residual Calculation ---
                % 1. Difference
                diff_rad = phase_model_rad - phase_data_rad;
                
                % 2. Handle Phase Wrapping (-pi to pi)
                % The robust way is to take the angle of the complex phasor
                diff_rad_wrapped = angle(exp(1i * diff_rad));
                
                % 3. Weighting
                % Weight by inverse of uncertainty (standard deviation)
                % s2p is 2*sigma, so sigma = s2p/2.
                weights = 1 ./ (uncertainty_rad / 2);
                
                % Handle cases where uncertainty is 0 or NaN
                weights(isnan(weights) | isinf(weights)) = 1; 
                
                % Weighted Residual
                res = diff_rad_wrapped .* weights;
                
                % Flatten and append
                residuals = [residuals; res(:)];
            end
        end
    end
    
    % Remove NaNs from residuals just in case
    residuals(isnan(residuals)) = [];
end



%% Change figure color
figHandles = findall(groot, 'Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);
    
    styleAndResizeFigure(fig)
end