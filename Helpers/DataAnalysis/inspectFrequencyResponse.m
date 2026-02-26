function inspectFrequencyResponse(data, LPM_data, sys_fitted, plot_title_suffix)
% inspectFrequencyResponse Generates a Bode plot comparing experimental data 
% (using Local Polynomial Method results) with a fitted model's response.
%
%   Inputs:
%       data              - iddata object containing input/output names and dimensions.
%       LPM_data          - (Optional) Struct array containing pre-calculated frequency response data.
%                           Can be empty [] if only plotting model.
%                           Expected fields per element (experiment):
%                             .freq      (Frequency vector [Hz])
%                             .magnitude (Magnitude array [Ny, Nu, Nf])
%                             .phase     (Phase array [Ny, Nu, Nf] in radians)
%                             .s2g       (2-sigma gain uncertainty)
%                             .s2p       (2-sigma phase uncertainty)
%       sys_fitted        - (Optional) Identify system model to plot against the data.
%       plot_title_suffix - (Optional) String appended to plot titles/filenames.

    arguments
        data
        LPM_data
        sys_fitted = [] % Optional
        plot_title_suffix = " " % Optional suffix for titles/filenames
    end

    % --- Hard coded settings --- %
    desired_line_width_model = 2.5;
    desired_line_width_LPM = 2;
    lpm_marker_size = 11; 
    lpm_stagger_amount = 0.035; % Fraction each channel is shifted for visibility.
    
    axis_label_font_size = 18; % Font size for X/Y-axis labels
    axis_tick_font_size = 13; % Font size for axis ticks
    
    %xlim_freq = [1e-0, 1e2]; % Frequency range in Hz
    xlim_freq = [5*1e-0, 2*1e2]; % Frequency range in Hz [5*1e-0, 2*1e2] is good.
    ylim_phase = [-360, 20]; % Phase range in degrees
    
    ylim_mag = false;%[-100, 25]; % ylim_mag = [-60, 10] or false if auto-scaling% Currently seeems to be broken. TODO FIX

    % --- Dimensions & Setup ---
    num_experiments = size(data, 4);
    num_outputs = size(data, 2);
    num_inputs = size(data, 3);

    % Colors for Outputs: Get standard 7, remove 3rd (Yellow), and cycle
    base_colors = lines(7);
    base_colors(3, :) = []; % Remove Yellow
    
    % Create enough colors for all outputs by repeating the base_colors
    output_colors = repmat(base_colors, ceil(num_outputs/size(base_colors,1)), 1);
    output_colors = output_colors(1:num_outputs, :);
    
    % --- Pre-calculate Model Response (if provided) ---
    if ~isempty(sys_fitted)
        % bode returns [Ny, Nu, Nw]
        [mag_all, phase_all, w_all] = bode(sys_fitted); 
        freq_model_hz = w_all / (2*pi);
    end

    % --- Loop over Inputs (Create one Figure per Input Channel) ---
    for j = 1:num_inputs
        
        fig_name = sprintf('Bode %s %s', data.InputName{j}, plot_title_suffix);
        figure('Name', fig_name);
        
        % --- Setup Axes (2 Subplots total: Mag & Phase) ---
        ax_mag = subplot(2, 1, 1);
        % MOVED UP: Set Axes properties first so they don't overwrite labels later
        set(ax_mag, 'XScale', 'log', 'FontSize', axis_tick_font_size);
        hold on; grid on;
        ylabel('Magnitude (dB)', 'FontWeight', 'bold', 'FontSize', axis_label_font_size);
        
        ax_phase = subplot(2, 1, 2);
        % MOVED UP: Set Axes properties first
        set(ax_phase, 'XScale', 'log', 'FontSize', axis_tick_font_size);
        hold on; grid on;
        ylabel('Phase (deg)', 'FontWeight', 'bold', 'FontSize', axis_label_font_size);
        xlabel('Frequency (Hz)', 'FontSize', axis_label_font_size);
        
        % --- Loop over Outputs (Plot all on the same axes) ---
        for i = 1:num_outputs
            current_color = output_colors(i, :);
            
            % --- 1. Plot Model (Reference) ---
            if ~isempty(sys_fitted)
                mag_model = squeeze(mag_all(i, j, :));
                phase_model = squeeze(phase_all(i, j, :));
                
                mag_model_db = 20*log10(mag_model);
                
                % Plot Model Lines
                plot(ax_mag, freq_model_hz, mag_model_db, '-', 'Color', current_color, 'LineWidth', desired_line_width_model, 'DisplayName', data.OutputName{i});
                plot(ax_phase, freq_model_hz, phase_model, '-', 'Color', current_color, 'LineWidth', desired_line_width_model, 'HandleVisibility', 'off');
            end
            
            % --- 2. Plot Experiments (From Stored LPM Data) ---
            if ~isempty(LPM_data)
                for k = 1:num_experiments
                    try
                        % Extract Stored Data
                        freqs = LPM_data(k).freq * (1 + (i-1)*lpm_stagger_amount);
                        
                        mag_linear = squeeze(LPM_data(k).magnitude(i, j, :));
                        phase_rad  = squeeze(LPM_data(k).phase(i, j, :));
                        
                        % Extract Uncertainty
                        if isfield(LPM_data(k), 's2g')
                            s2g = squeeze(LPM_data(k).s2g(i, j, :));
                            s2p = squeeze(LPM_data(k).s2p(i, j, :));
                        else
                            s2g = zeros(size(mag_linear));
                            s2p = zeros(size(phase_rad));
                        end
                        
                        % -- Process Magnitude (dB) --
                        mag_db = 20*log10(mag_linear);
                        
                        % Calc Error Bars (Linear -> dB)
                        upper_db = 20*log10(mag_linear + s2g);
                        lower_linear = mag_linear - s2g;
                        lower_linear(lower_linear <= 0) = 1e-9; 
                        lower_db = 20*log10(lower_linear);
                        
                        err_up = upper_db - mag_db;
                        err_down = mag_db - lower_db;
                        
                        % -- Process Phase (Deg) --
                        phase_deg = rad2deg(phase_rad);
                        phase_err = rad2deg(s2p);
                        
                        % -- Phase Matching --
                        for phase_index = 1:length(phase_deg)
                            window_center = mean(ylim_phase);
                            N_shifts = round((window_center - phase_deg(phase_index)) / 360);
                            phase_deg(phase_index) = phase_deg(phase_index) + N_shifts * 360;
                        end
                        
                        % Handle Legend
                        if isempty(sys_fitted) && k == 1
                             disp_name = data.OutputName{i};
                             handle_vis = 'on';
                        else
                             disp_name = '';
                             handle_vis = 'off';
                        end
    
                        % -- Plotting --
                        errorbar(ax_mag, freqs, mag_db, err_down, err_up, 'x', ...
                            'Color', current_color, 'LineWidth', desired_line_width_LPM, ...
                            'MarkerSize', lpm_marker_size, 'DisplayName', disp_name, 'HandleVisibility', handle_vis);
                        
                        errorbar(ax_phase, freqs, phase_deg, phase_err, phase_err, 'x', ...
                            'Color', current_color, 'LineWidth', desired_line_width_LPM, ...
                            'MarkerSize', lpm_marker_size, 'HandleVisibility', 'off');
                        
                    catch ME
                        warning('Error plotting LPM for Input %d, Output %d, Exp %d: %s', j, i, k, ME.message);
                    end
                end
            end
        end
        
        % Format Axes (applied once per figure)
        % Note: XScale and FontSize were moved to the top
        %legend(ax_mag, 'show', 'Location', 'best');
        legend(ax_mag, 'show', 'Location', 'northeast');% Best, SouthWest, SouthEast, NorthWest, NorthEast, BestOutside, EastOutside, NorthEastOutside, NorthWestOutside, SouthEastOutside, SouthWestOutside

        xlim(ax_mag, xlim_freq);
        xlim(ax_phase, xlim_freq);
        
        % Apply Magnitude limits if not set to false
        if ~isequal(ylim_mag, false)
            ylim(ax_mag, ylim_mag);
        end
        
        ylim(ax_phase, ylim_phase);
        linkaxes([ax_mag, ax_phase], 'x');
    end
end