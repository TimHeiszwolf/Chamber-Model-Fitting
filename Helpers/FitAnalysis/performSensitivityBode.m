function performSensitivityBode(sys_grey, settings, perturbation_fraction, mag_threshold, phase_threshold)
    % performSensitivityAnalysisWithInterval Analyzes model sensitivity and visualizes
    % it with both bar charts and shaded interval Bode plots for +/- perturbations.
    %
    %   Inputs:
    %       sys_grey              - The fitted idgrey model from greyest.
    %       settings              - The settings struct, containing the chamber_model object.
    %       perturbation_fraction - The fraction to perturb each parameter by.
    %       mag_threshold         - Threshold in dB to display parameter in filtered plot.
    %       phase_threshold       - Threshold in degrees to display parameter in filtered plot.
    
    arguments
        sys_grey
        settings
        perturbation_fraction = 0.05;
        mag_threshold = 0.5;
        phase_threshold = 1.00;
    end

    % --- Hard coded settings (Matched to inspectFrequencyResponse) ---
    desired_line_width_baseline = 2.5;
    axis_label_font_size = 18; % Font size for X/Y-axis labels
    axis_tick_font_size = 13;  % Font size for axis ticks
    xlim_freq = [5*1e-0, 2*1e2]; % Frequency range in Hz
    % -----------------------------------------------------------------

    %disp(sprintf('Starting combined sensitivity analysis (+/- %.2f%% perturbation)...', perturbation_percentage));

    % --- Step 1: Get baseline information ---
    [p_vec_baseline, ~] = getpvec(sys_grey);
    param_names = {sys_grey.Structure.Parameters.Name};
    num_params = length(param_names);
    
    % --- Generate Distinct Colors for All Parameters ---
    param_colors = lines(num_params); 
    
    tf_baseline = tf(sys_grey);
    [num_outputs, num_inputs] = size(tf_baseline);
    
    % Use the baseline model's response to define a common frequency vector
    [mag_base, phase_base, freqs] = bode(tf_baseline);
    
    freqs_hz = freqs / (2*pi); % Convert to Hz for plotting
    
    % --- Step 2: Loop through I/O pairs and parameters to gather all perturbed responses ---
    for in_idx = 1:num_inputs
        for out_idx = 1:num_outputs
            
            % Initialize matrices to store the Bode responses and metrics for this I/O pair
            all_perturbed_mags_pos = zeros(length(freqs), num_params);
            all_perturbed_phases_pos = zeros(length(freqs), num_params);
            all_perturbed_mags_neg = zeros(length(freqs), num_params);
            all_perturbed_phases_neg = zeros(length(freqs), num_params);
            rms_mag_sensitivity = zeros(1, num_params);
            rms_phase_sensitivity = zeros(1, num_params);
            
            % For this specific I/O pair, loop through all parameters
            for i = 1:num_params
                % --- Positive Perturbation ---
                p_vec_pos = p_vec_baseline;
                p_vec_pos(i) = p_vec_pos(i) * (1 + perturbation_fraction);
                param_cell_pos = settings.chamber_model.default_parameters;
                param_cell_pos(:, 2) = num2cell(p_vec_pos);
                [A_p, B_p, C_p, D_p] = settings.chamber_model.getMatrices(param_cell_pos{:,2});
                sys_pos = ss(A_p, B_p, C_p, D_p);
                [mag_p_pos, phase_p_pos] = bode(sys_pos, freqs);
                
                % --- Negative Perturbation ---
                p_vec_neg = p_vec_baseline;
                %p_vec_neg(i) = p_vec_neg(i) / (1 + perturbation_fraction);
                p_vec_neg(i) = p_vec_neg(i) * (1 - perturbation_fraction);
                param_cell_neg = settings.chamber_model.default_parameters;
                param_cell_neg(:, 2) = num2cell(p_vec_neg);
                [A_p, B_p, C_p, D_p] = settings.chamber_model.getMatrices(param_cell_neg{:,2});
                sys_neg = ss(A_p, B_p, C_p, D_p);
                [mag_p_neg, phase_p_neg] = bode(sys_neg, freqs);

                % Store full response curves for interval plots
                all_perturbed_mags_pos(:, i) = squeeze(mag_p_pos(out_idx, in_idx, :));
                all_perturbed_phases_pos(:, i) = squeeze(phase_p_pos(out_idx, in_idx, :));
                all_perturbed_mags_neg(:, i) = squeeze(mag_p_neg(out_idx, in_idx, :));
                all_perturbed_phases_neg(:, i) = squeeze(phase_p_neg(out_idx, in_idx, :));

                % --- Calculate RMS metric for Bar Charts (using positive perturbation) ---
                mag_base_vec = squeeze(mag_base(out_idx, in_idx, :));
                phase_base_vec = squeeze(phase_base(out_idx, in_idx, :));
                mag_diff_dB = 20*log10(all_perturbed_mags_pos(:, i)) - 20*log10(mag_base_vec);
                phase_diff_deg = all_perturbed_phases_pos(:, i) - phase_base_vec;
                rms_mag_sensitivity(i) = rms(mag_diff_dB);
                rms_phase_sensitivity(i) = rms(phase_diff_deg);
            end
            
            % --- Step 3a: Visualize sensitivity with Bar Charts ---
            figure;
            bar(rms_mag_sensitivity);
            title(sprintf('Magnitude Sensitivity: Input ''%s'' to Output ''%s''', tf_baseline.InputName{in_idx}, tf_baseline.OutputName{out_idx}), 'FontSize', axis_tick_font_size);
            ylabel('RMS Change in Magnitude (dB)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            xlabel('Model Parameters', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            set(gca, 'XTick', 1:num_params, 'XTickLabel', param_names, 'FontSize', axis_tick_font_size);
            xtickangle(45);
            grid on;

            figure;
            bar(rms_phase_sensitivity);
            title(sprintf('Phase Sensitivity: Input ''%s'' to Output ''%s''', tf_baseline.InputName{in_idx}, tf_baseline.OutputName{out_idx}), 'FontSize', axis_tick_font_size);
            ylabel('RMS Change in Phase (degrees)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            xlabel('Model Parameters', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            set(gca, 'XTick', 1:num_params, 'XTickLabel', param_names, 'FontSize', axis_tick_font_size);
            xtickangle(45);
            grid on;

            % --- Step 3b: Visualize sensitivity with Interval Plots (ALL PARAMETERS) ---
            figure('Name', sprintf('Sensitivity %s to %s', tf_baseline.InputName{in_idx}, tf_baseline.OutputName{out_idx}));
            
            % Get the baseline data for plotting
            mag_base_data = squeeze(20*log10(mag_base(out_idx, in_idx, :)));
            phase_base_data = squeeze(phase_base(out_idx, in_idx, :));
            
            % Prepare for legend
            legend_handles = gobjects(num_params + 1, 1);
            legend_labels = cell(num_params + 1, 1);
            legend_labels{1} = 'Baseline System';

            % Magnitude Plot
            subplot(2, 1, 1);
            hold on;
            legend_handles(1) = plot(freqs_hz, mag_base_data, 'k-', 'LineWidth', desired_line_width_baseline);

            for i = 1:num_params
                mag_pos_data = 20*log10(all_perturbed_mags_pos(:, i));
                mag_neg_data = 20*log10(all_perturbed_mags_neg(:, i));
                mag_lower_bound = min([mag_base_data, mag_pos_data, mag_neg_data], [], 2);
                mag_upper_bound = max([mag_base_data, mag_pos_data, mag_neg_data], [], 2);
                x_poly = [freqs_hz; flipud(freqs_hz)];
                y_poly = [mag_lower_bound; flipud(mag_upper_bound)];
                
                h = fill(x_poly, y_poly, param_colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
                legend_handles(i+1) = h;
                legend_labels{i+1} = sprintf('%s', param_names{i});
            end

            set(gca, 'XScale', 'log', 'FontSize', axis_tick_font_size);
            grid on;
            %title(sprintf('Bode Sensitivity Intervals (%s %%): Input ''%s'' to Output ''%s''', string(100 * perturbation_fraction), tf_baseline.InputName{in_idx}, tf_baseline.OutputName{out_idx}), 'FontSize', axis_tick_font_size);
            ylabel('Magnitude (dB)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            legend(legend_handles, legend_labels, 'Location', 'northeast');
            xlim(xlim_freq);
            
            % Phase Plot
            subplot(2, 1, 2);
            hold on;
            plot(freqs_hz, phase_base_data, 'k-', 'LineWidth', desired_line_width_baseline);
            for i = 1:num_params
                phase_pos_data = all_perturbed_phases_pos(:, i);
                phase_neg_data = all_perturbed_phases_neg(:, i);
                phase_lower_bound = min([phase_base_data, phase_pos_data, phase_neg_data], [], 2);
                phase_upper_bound = max([phase_base_data, phase_pos_data, phase_neg_data], [], 2);
                x_poly = [freqs_hz; flipud(freqs_hz)];
                y_poly = [phase_lower_bound; flipud(phase_upper_bound)];
                
                fill(x_poly, y_poly, param_colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            end
            
            set(gca, 'XScale', 'log', 'FontSize', axis_tick_font_size);
            grid on;
            ylabel('Phase (deg)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            xlabel('Frequency (Hz)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
            xlim(xlim_freq);

            % --- Step 3c: Visualize filtered sensitivity (Thresholded) ---
            % Find indices of parameters exceeding thresholds
            sig_mag_idx = find(rms_mag_sensitivity > mag_threshold);
            sig_phase_idx = find(rms_phase_sensitivity > phase_threshold);

            if ~isempty(sig_mag_idx) || ~isempty(sig_phase_idx)
                figure('Name', sprintf('Significant Sensitivity %s to %s', tf_baseline.InputName{in_idx}, tf_baseline.OutputName{out_idx}));
                
                % --- Filtered Magnitude Plot ---
                subplot(2, 1, 1);
                hold on;
                % Plot Baseline
                plot(freqs_hz, mag_base_data, 'k-', 'LineWidth', desired_line_width_baseline, 'DisplayName', 'Model');
                
                % Loop only through significant magnitude parameters
                for k = 1:length(sig_mag_idx)
                    idx = sig_mag_idx(k);
                    mag_pos_data = 20*log10(all_perturbed_mags_pos(:, idx));
                    mag_neg_data = 20*log10(all_perturbed_mags_neg(:, idx));
                    mag_lower = min([mag_base_data, mag_pos_data, mag_neg_data], [], 2);
                    mag_upper = max([mag_base_data, mag_pos_data, mag_neg_data], [], 2);
                    
                    x_poly = [freqs_hz; flipud(freqs_hz)];
                    y_poly = [mag_lower; flipud(mag_upper)];
                    
                    fill(x_poly, y_poly, param_colors(idx, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
                         'DisplayName', sprintf('%s', param_names{idx}));
                end
                set(gca, 'XScale', 'log', 'FontSize', axis_tick_font_size); grid on; xlim(xlim_freq);
                %title(sprintf('Significant Parameters (>%.1fdB, >%.1fdeg): %s to %s', mag_threshold, phase_threshold, tf_baseline.InputName{in_idx}, tf_baseline.OutputName{out_idx}), 'FontSize', axis_tick_font_size);
                ylabel('Magnitude (dB)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
                legend('show', 'Location', 'northeast');

                % --- Filtered Phase Plot ---
                subplot(2, 1, 2);
                hold on;
                plot(freqs_hz, phase_base_data, 'k-', 'LineWidth', desired_line_width_baseline, 'DisplayName', 'Baseline System');
                
                for k = 1:length(sig_phase_idx)
                    idx = sig_phase_idx(k);
                    phase_pos_data = all_perturbed_phases_pos(:, idx);
                    phase_neg_data = all_perturbed_phases_neg(:, idx);
                    phase_lower = min([phase_base_data, phase_pos_data, phase_neg_data], [], 2);
                    phase_upper = max([phase_base_data, phase_pos_data, phase_neg_data], [], 2);
                    
                    x_poly = [freqs_hz; flipud(freqs_hz)];
                    y_poly = [phase_lower; flipud(phase_upper)];
                    
                    fill(x_poly, y_poly, param_colors(idx, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
                         'DisplayName', sprintf('%s', param_names{idx}));
                end
                set(gca, 'XScale', 'log', 'FontSize', axis_tick_font_size); grid on; xlim(xlim_freq);
                ylabel('Phase (deg)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold'); xlabel('Frequency (Hz)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
                legend('show', 'Location', 'northeast');
            end
        end
    end
    %disp('Combined sensitivity analysis complete. All plots have been generated.');
end