function inspectData(data, settings, inspect_level)
% inspectData Visualizes input-output data for system identification with varying levels of detail.
%
% Syntax:
%   inspectData(data, settings)               % Uses settings.inspect_data
%   inspectData(data, settings, inspect_level) % Overrides settings
%
% Description:
%   Generates plots to inspect data quality, frequency response, correlations, 
%   and model order requirements.
%
% Inputs:
%   data          - An iddata object containing the input-output signals.
%   settings      - Struct containing settings. Must contain 'LPM_data' for Level 2.
%   inspect_level - (Optional) Integer (0-5) specifying detail level:
%                     0: None
%                     1: Time Domain (Input/Output plots)
%                     2: Frequency Domain (Bode plot via LPM)
%                     3: Correlations & Spectral Analysis (SPA)
%                     4: Full LPM Analysis (Detailed per-channel plots)
%                     5: Hankel Singular Values (Model order check)
%                     (Higher levels include all lower levels)

%% Determine inspection level
if nargin < 3
    inspect_level = settings.inspect_data;
end

% --- Hardcoded settings ---
lineWidth = 4;
axis_label_font_size = 18; % Font size for X/Y-axis labels
axis_tick_font_size = 13;  % Font size for axis numbers
legend_font_size = 11;     % For the text boxes and subplot titles
axis_line_width = 1.5;      % Width for Axis lines and Grid lines (default of MATlab is 0.5)

%% Plot data in the time domain
if inspect_level >= 1
    % --- Plot Output Channels ---
    figure('Name', 'Output Channels');
    num_outputs = size(data.y, 2);
    %sgtitle('Output Channels');%Title could be commented out for if you want it easier use it in the report.
    ax_list = [];
    for i = 1:num_outputs
        ax = subplot(num_outputs, 1, i);
        ax_list = [ax_list, ax];
        plot(data.SamplingInstants, data.y(:, i), 'LineWidth', lineWidth);

        % --- Styling ---
        % Add styled text box for channel name
        text(ax, 0.97, 0.02, data.OutputName{i}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', ...
            'Color', 'k', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'k', ...
            'Margin', axis_line_width, ...
            'FontSize', legend_font_size);

        set(ax, 'FontSize', axis_tick_font_size); % Set Ticks FIRST

        % X-Axis (Only show on bottom plot)
        if i == num_outputs
            xlabel('Time (s)', 'FontSize', axis_label_font_size);
        else
            xlabel('');
            set(ax, 'XTickLabel', []); % Hide numbers
        end

        grid on;
        title(''); % Remove default title
    end
    linkaxes(ax_list, 'x');
    
    % --- Add Global Y-Label for Output Channels ---
    han = axes(gcf, 'visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Amplitude (au.)', 'FontSize', axis_label_font_size);
    
    % --- Plot Input Channels ---
    figure('Name', 'Input Channels');
    num_inputs = size(data.u, 2);
    %sgtitle('Input Channels');%Title could be commented out for if you want it easier use it in the report.
    ax_list = [];
    for i = 1:num_inputs
        ax = subplot(num_inputs, 1, i);
        ax_list = [ax_list, ax];
        plot(data.SamplingInstants, data.u(:, i), 'LineWidth', lineWidth);
        
        % --- Styling ---
        text(ax, 0.97, 0.02, data.InputName{i}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', ...
            'Color', 'k', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'k', ...
            'LineWidth', 1.5, ...
            'Margin', 1, ...
            'FontSize', legend_font_size);

        set(ax, 'FontSize', axis_tick_font_size); % Set Ticks FIRST

        if i == num_inputs
            xlabel('Time (s)', 'FontSize', axis_label_font_size);
        else
            xlabel('');
            set(ax, 'XTickLabel', []);
        end

        grid on;
        title('');
    end
    linkaxes(ax_list, 'x');

    % --- Add Global Y-Label for Input Channels ---
    han = axes(gcf, 'visible', 'off');
    han.YLabel.Visible = 'on';
    ylabel(han, 'Amplitude (au.)', 'FontSize', axis_label_font_size);
end

%% Plot Frequency Response (Bode) using LPM Data
if inspect_level >= 2
    inspectFrequencyResponse(data, settings.LPM_data, [], "Inspection");% Use the separate function for clean Bode plots
    % This does have a bug. If multiple data sets are used then LPM_data is an array. The software will then only show the first one. Can be circumvented.
end

%% Plot data with Spa and check correlations
if inspect_level >= 3
    num_outputs = size(data.y, 2);
    num_inputs = size(data.u, 2);

    % --- SPA Analysis ---
    gs = spa(data);
    figure('Name', 'Spa');
    bodeplot(gs);
    grid on;
    %title('Bode Plot from iddata using spa');%Title could be commented out for if you want it easier use it in the report.

    % --- Correlation Analysis ---
    total_plots_per_output = num_inputs + num_outputs - 1;
    
    for i = 1:num_outputs
        fig = figure;
        %sgtitle(sprintf('Correlation Analysis for Output: %s', data.OutputName{i}), 'Interpreter', 'none');%Title could be commented out for if you want it easier use it in the report.
        sgtitle(sprintf('Output: %s', data.OutputName{i}), 'Interpreter', 'none', ...
            'Color', 'k', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');

        plot_idx = 1;
        ax_list = [];
    
        % Correlation with each input signal
        for j = 1:num_inputs
            ax = subplot(total_plots_per_output, 1, plot_idx);
            ax_list = [ax_list, ax];
            crosscorr(data.y(:, i), data.u(:, j));
            legend(ax, 'off'); % Remove default legend
            hStem = findobj(gca, 'Type', 'stem');% Find markers
            set(hStem, 'MarkerSize', 4);% Decrease marker size
            ylim([-1 1]);
            
            % Snap x-axis to data range
            x_data = hStem.XData;
            xlim([min(x_data), max(x_data)]);

            title(sprintf('vs. Input: %s', data.InputName{j}), 'Interpreter', 'none', 'FontSize', legend_font_size);
            
            % Clean axes & Fonts
            set(ax, 'FontSize', axis_tick_font_size);
            ylabel(''); % Remove individual Y label

            if plot_idx == total_plots_per_output
                xlabel('Lag (samples)', 'FontSize', axis_label_font_size);
            else
                xlabel('');
                set(ax, 'XTickLabel', []);
            end
            
            plot_idx = plot_idx + 1;
        end
    
        % Correlation with other output signals
        for k = 1:num_outputs
            if i == k
                continue; % Skip auto-correlation
            end
            ax = subplot(total_plots_per_output, 1, plot_idx);
            ax_list = [ax_list, ax];
            crosscorr(data.y(:, i), data.y(:, k));
            legend(ax, 'off'); % Remove default legend
            hStem = findobj(gca, 'Type', 'stem');% Find markers
            set(hStem, 'MarkerSize', 4);% Decrease marker size
            ylim([-1 1]);
            
            % Snap x-axis to data range
            x_data = hStem.XData;
            xlim([min(x_data), max(x_data)]);

            title(sprintf('vs. Output: %s', data.OutputName{k}), 'Interpreter', 'none', 'FontSize', legend_font_size);
            
            % Clean axes & Fonts
            set(ax, 'FontSize', axis_tick_font_size);
            ylabel(''); % Remove individual Y label
            
            if plot_idx == total_plots_per_output
                xlabel('Lag (samples)', 'FontSize', axis_label_font_size);
            else
                xlabel('');
                set(ax, 'XTickLabel', []);
            end
            
            plot_idx = plot_idx + 1;
        end
        
        % Add global Y label using an invisible axes
        han = axes(fig, 'visible', 'off');
        han.YLabel.Visible = 'on';
        ylabel(han, 'Sample Cross-Correlation', 'FontSize', axis_label_font_size);
    end
end

%% Do LPM
if inspect_level >= 4
    disp('Starting frequency domain and correlation analysis...');
    
    % --- LPM Analysis ---
    plotsettings.fig_on = 1;
    plotsettings.phaserange = [-300 0];
    plotsettings.freqrange_FRF = [0.5 200];
    plotsettings.freqrange_lines = [0.5 100];
    plotsettings.correct_phase = 0;
    plotsettings.xscale = 'linear';
    
    num_inputs = size(data.u, 2);
    num_outputs = size(data.y, 2);

    for j = 1:num_inputs % Loop over each input
        for i = 1:num_outputs % Loop over each output
            fs = 1 / data.Ts;
            u_j = data.u(:, j);
            y_i = data.y(:, i);
            
            uc = frf_cutsignal(u_j, data.SamplingInstants, settings.start_time, settings.detrend_data.f_exc(1), fs, settings.detrend_data.P, settings.end_time);
            yc = frf_cutsignal(y_i, data.SamplingInstants, settings.start_time, settings.detrend_data.f_exc(1), fs, settings.detrend_data.P, settings.end_time);
    
            [~, ~, ~, ~, ~, ~] = frf_lpm_main(uc, yc, settings.detrend_data.f0, settings.detrend_data.f_exc, settings.detrend_data.P, plotsettings.fig_on);
            sgtitle(sprintf('LPM FRF: From Input ''%s'' to Output ''%s''', data.InputName{j}, data.OutputName{i}));
        end
    end
end

%% Check Hankel matrices for model order
if inspect_level >= 5
    order_range = 1:10;
    figure;
    n4sid(data, order_range);
    title('Hankel Singular Values');
    xlabel('Model Order');
    ylabel('Singular Value Magnitude');
    grid on;
end

end