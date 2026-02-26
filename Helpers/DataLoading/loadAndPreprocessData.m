function [data, training_data, validation_data, shot, settings] = loadAndPreprocessData(settings)
% loadAndPreprocessData Loads, preprocesses, merges, and splits experiment data.
%
%   [data, training_data, validation_data, shot, settings] = loadAndPreprocessData(settings)
%
%   Description:
%       Iterates through a list of experiment filenames defined in 'settings'.
%       For each file, it:
%       1. Loads the raw data and experiment parameters (f_exc, P, etc.).
%       2. Preprocesses the time-domain data (detrending, normalization) using 
%          the chamber model's specific `preProccesData` method.
%       3. Computes the Frequency Response Function (FRF) using the Local 
%          Polynomial Method (LPM) to establish a "ground truth" Bode plot 
%          with uncertainty bounds.
%       4. Merges all experiments into a single multi-experiment `iddata` object.
%       5. Splits the merged data into training and validation sets.
%
%   Inputs:
%       settings - Struct containing configuration options:
%           .data_filenames   - Cell array of file paths (strings).
%           .start_time       - Scalar or array of start times [s].
%           .end_time         - Scalar or array of end times [s].
%           .chamber_model    - Model object with a .preProccesData() method.
%           .validation_split - Fraction (0-1) of data to reserve for validation.
%           .normalize_data   - (Optional) Normalization flag.
%
%   Outputs:
%       data            - Merged `iddata` object containing all experiments.
%       training_data   - `iddata` subset used for model fitting.
%       validation_data - `iddata` subset used for validation.
%       shot            - Array of shot numbers loaded.
%       settings        - Updated settings struct, now including:
%           .LPM_data   - Struct array of frequency domain results per experiment.

    %% Initialization
    settings.normalize_data_original_setting = settings.normalize_data;
    merged_data = [];
    shot_list = [];

    % Initialize LPM (Frequency Domain) storage structure
    settings.LPM_data = struct('freq', {}, 'magnitude', {}, 'phase', {}, 's2g', {}, 's2p', {}, 'f0', {}, 'P', {});

    config_start_times = settings.start_time;
    config_end_times   = settings.end_time;

    %% Loop through all experiment files
    for i = 1:length(settings.data_filenames)
        current_file = settings.data_filenames{i};
        fprintf('Processing file (%d/%d): %s \n', i, length(settings.data_filenames), current_file);
        
        % --- Handle Time Window Settings ---
        % Allow global settings (scalar) or per-experiment settings (array)
        if length(config_start_times) > 1
            settings.start_time = config_start_times(i);
        else
            settings.start_time = config_start_times;
        end
        
        if length(config_end_times) > 1
            settings.end_time = config_end_times(i);
        else
            settings.end_time = config_end_times;
        end
        
        % --- Load Raw Data ---
        try
            loaded = load(current_file, "data", "shot", "f_exc", "P", "f0");
        catch ME
            warning('Could not load file %s: %s', current_file, ME.message);
            continue;
        end
        
        % --- Pre-processing ---
        % Pass parameters needed for detrending to the settings struct
        settings.detrend_data.f0 = loaded.f0;
        settings.detrend_data.f_exc = loaded.f_exc;
        settings.detrend_data.P = loaded.P;
        
        % Clean, detrend, and normalize the time-domain data
        [temp_data, settings] = settings.chamber_model.preProccesData(loaded.data, settings);
        
        % --- Calculate LPM Frequency Response (Ground Truth) ---
        % We calculate the FRF now while we have the exact raw data timing 
        % and excitation parameters available.
        fprintf('   Calculating LPM Frequency Response...\n');
        num_inputs_shot = size(temp_data, 3);
        num_outputs_shot = size(temp_data, 2);
        num_freqs = length(loaded.f_exc);
        
        % Initialize 3D storage matrices: [Outputs x Inputs x Frequencies]
        lpm_mag   = zeros(num_outputs_shot, num_inputs_shot, num_freqs);
        lpm_phase = zeros(num_outputs_shot, num_inputs_shot, num_freqs);
        lpm_s2g   = zeros(num_outputs_shot, num_inputs_shot, num_freqs); % Gain Uncertainty
        lpm_s2p   = zeros(num_outputs_shot, num_inputs_shot, num_freqs); % Phase Uncertainty
        
        fs = 1/temp_data.Ts;
        t = temp_data.SamplingInstants;
        
        % Iterate through every Input-Output pair to compute the transfer function
        for u_idx = 1:num_inputs_shot
            for y_idx = 1:num_outputs_shot
                u_sig = temp_data.u(:, u_idx);
                y_sig = temp_data.y(:, y_idx);
                
                % Cut signals to exact periods (P) to minimize leakage
                uc = frf_cutsignal(u_sig, t, settings.start_time, loaded.f_exc(1), fs, loaded.P, settings.end_time);
                yc = frf_cutsignal(y_sig, t, settings.start_time, loaded.f_exc(1), fs, loaded.P, settings.end_time);
                
                % Run Local Polynomial Method (LPM)
                [FRF_res, ~, ~, ~, ~, ~] = frf_lpm_main(uc, yc, loaded.f0, loaded.f_exc, loaded.P, 0);
                
                % Store results
                lpm_mag(y_idx, u_idx, :)   = abs(FRF_res.TF);       % Magnitude (Linear)
                lpm_phase(y_idx, u_idx, :) = angle(FRF_res.TF);     % Phase (Radians)
                lpm_s2g(y_idx, u_idx, :)   = FRF_res.s2g;           % 2-sigma Gain Error
                lpm_s2p(y_idx, u_idx, :)   = FRF_res.s2p;           % 2-sigma Phase Error
            end
        end
        
        % Store frequency data in the settings struct array (one entry per experiment)
        settings.LPM_data(i).freq      = loaded.f_exc;
        settings.LPM_data(i).magnitude = lpm_mag;   
        settings.LPM_data(i).phase     = lpm_phase; 
        settings.LPM_data(i).s2g       = lpm_s2g;
        settings.LPM_data(i).s2p       = lpm_s2p;
        settings.LPM_data(i).f0        = loaded.f0;
        settings.LPM_data(i).P         = loaded.P;
        
        % --- Merge Data ---
        if isempty(merged_data)
            merged_data = temp_data;
        else
            % combine experiments into a multi-experiment iddata object
            merged_data = merge(merged_data, temp_data);
        end
        
        shot_list = [shot_list, loaded.shot];
    end

    if isempty(merged_data)
        error('No data could be loaded or merged.');
    end

    data = merged_data;
    shot = shot_list; 

    % Restore the full configuration arrays to settings for record-keeping
    settings.start_time = config_start_times;
    settings.end_time = config_end_times;

    %% Split into Training and Validation sets
    % Uses helper function to split based on time or experiment index
    [training_data, validation_data] = splitData(data, settings.validation_split);

end