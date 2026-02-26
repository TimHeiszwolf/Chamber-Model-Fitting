function styleTimeDomainComparison(data, validation_split, plot_title_suffix)
%% Post-processes a MATLAB compare() plot for reports.
% This function applies consistent styling, thickens lines, manages axis labels, 
% and draws a vertical line at the training/validation data split point.
%
%   USAGE:
%       compare(data, sys1, sys2, ...);
%       styleTimeDomainComparison(data, 0.7, "Experiment A");
%
%   INPUTS:
%       data              - The iddata object used in the original compare() call.
%       validation_split  - (optional)Double (0-1) representing the split ratio (from settings.validation_split. Default is 1.
%       plot_title_suffix - (Optional) String appended to the figure window name.

    arguments
        data
        validation_split double = 1
        plot_title_suffix string = " "
    end

    % --- Configuration ---
    desired_line_width = 4;
    split_line_color = [0.85, 0.325, 0.098]; % Orange-red
    axis_label_font_size = 18; % Font size for X/Y-axis labels
    axis_tick_font_size = 13;   % Font size for axis numbers
    legend_font_size = 11;
    axis_line_width = 1.5;      % Width for Axis lines and Grid lines (default of MATlab is 0.5)

    % Outdates configuration
    yaxis_label_rotation = 90;% Rotation of the axis label, 0 is horizontal, 90 is vertical.
    
    % --- Determine Split Time ---
    % We calculate where the training data ends based on the split ratio
    [training_data, ~] = splitData(data, validation_split);
    
    if ~isempty(training_data)
        split_time = training_data.SamplingInstants(end);
    else
        split_time = 0;
    end

    % --- Apply Styling to Current Figure ---
    fig_handle = gcf;

    % Remove default "Simulated Response Comparison" title added by compare()
    all_axes = findall(fig_handle, 'Type', 'axes');
    for i = 1:numel(all_axes)
         if contains(string(all_axes(i).Title.String), "Response Comparison", "IgnoreCase", true)
             all_axes(i).Title.String = ""; 
         end
    end
    
    new_title = sprintf('Time-Domain %s', plot_title_suffix);
    set(fig_handle, 'Name', new_title);% Set the window name
    %sgtitle(fig_handle, new_title, 'Interpreter', 'none'); %Title could be commented out for if you want it easier use it in the report.
    
    % 2. Thick Lines
    % Find all data lines (excluding the split line we might have just added)
    all_lines = findall(fig_handle, 'Type', 'line');
    set(all_lines, 'LineWidth', desired_line_width);

    % 3. Style Axes and Add Split Line
    % Find all plot axes in the figure (excluding legends)
    all_plot_axes = findall(fig_handle, 'Type', 'axes', '-not', 'Tag', 'legend');

    for k = 1:numel(all_plot_axes)
        ax = all_plot_axes(k);
        
        % -- Font and Label Styling --
        set(ax, 'FontSize', axis_tick_font_size);
        set(ax.XLabel, 'FontSize', axis_label_font_size);
        %ax.XLabel.String = "Time (s)";%REMOVE OR  ADJUST IF YOU WANT DIFFERENT X-AXIS LABEL, HARD-CODED HERE BECAUSE I COULDN'T FIGURE OUT WERE IT WAS OTHERWISE SET.
        set(ax.YLabel, 'FontSize', axis_label_font_size);
        set(ax.Title, 'FontSize', axis_label_font_size);
        
        % Move channel name from Y-Label to a styled Text Box inside the plot
        channel_name = string(ax.YLabel.String);
        if channel_name ~= "" && ~contains(channel_name, "Amplitude")
            % Clear standard title
            ax.Title.String = "";
            
            % Create text box in bottom right
            % Units 'normalized': 0,0 is bottom-left, 1,1 is top-right.
            text(ax, 0.97, 0.02, channel_name, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'bottom', ...
                'Color', 'k', ...
                'BackgroundColor', 'w', ...
                'EdgeColor', 'k', ...
                'LineWidth', axis_line_width, ...
                'Margin', 1, ...
                'FontSize', legend_font_size); % Use axis size for channel labels
        end

        % Set standard Y-axis label
        if ~contains(string(ax.YLabel.String), "Amplitude", "IgnoreCase", true)%Remove everything but the standard amplitude label
            ax.YLabel.String = "";
            %set(ax.YLabel, 'Rotation', yaxis_label_rotation, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        else
            ax.YLabel.String = "Amplitude (au.)";% If it does contain the amplitude it must be the "global" axis label and we fix it.
            set(ax.YLabel, 'Rotation', yaxis_label_rotation, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end

        % -- Split Line --
        if validation_split > 0 && validation_split < 1
            hold(ax, 'on');
            xline(ax, split_time, '--', 'LineWidth', desired_line_width,'Color', split_line_color, 'Alpha', 0.9, 'HandleVisibility', 'off'); % Hide from legend
            hold(ax, 'off');
        end
        grid(ax, 'on'); % Ensure grid is on
    end

    % 4. Clean up Legend Text by updating the underlying objects' DisplayNames
    all_plot_objects = findall(fig_handle, '-property', 'DisplayName');
    
    for k = 1:numel(all_plot_objects)
        obj = all_plot_objects(k);
        current_str = obj.DisplayName;
        
        if isempty(current_str)
            continue;
        end

        % Regex magic to fix legend labels
        current_str = regexprep(current_str, ':\s*-?[\d\.]+%?\s*$', '');% Remove fit percentage
        current_str = regexprep(current_str, '\s*\([^)]+\)\s*$', '');% Remove channel name (basically everything in parathesis)
        
        if contains(current_str, 'Validation data', 'IgnoreCase', true) || startsWith(current_str, 'Data ', 'IgnoreCase', true)
            current_str = 'Data';
        end
        
        obj.DisplayName = current_str;
    end

    % --- 5. Single Legend Logic ---
    % Sort axes by position (Top to Bottom) to identify the first channel.
    % We reuse 'all_plot_axes' found in step 3.
    if ~isempty(all_plot_axes)
        % Extract positions. Position is [left, bottom, width, height]
        % We sort by 'bottom' (index 2) in descending order to find the top plot.
        positions = vertcat(all_plot_axes.Position);
        [~, sortIdx] = sort(positions(:, 2), 'descend');
        sorted_axes = all_plot_axes(sortIdx);

        % Keep legend on the first (top) axis, remove others
        legend(sorted_axes(1), 'show'); 
        for k = 2:numel(sorted_axes)
             legend(sorted_axes(k), 'off');
        end
    end
end