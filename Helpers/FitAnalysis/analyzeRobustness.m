function analyzeRobustness(all_fits, sys_baseline, batch_strategy, settings)
% analyzeRobustness Visualizes parameter stability across multiple runs.
%
%   analyzeRobustness(all_fits, sys_baseline, batch_strategy, settings)
%
%   Inputs:
%       all_fits       - Struct array containing results from multiple runs.
%       sys_baseline   - The reference system (usually the best fit).
%       batch_strategy - The cell array defining which params were in which batch.
%       settings       - Struct containing .margin_factors.

    fprintf('\n--- Analyzing Robustness ---\n');

    all_params = sys_baseline.Structure.Parameters;
    all_param_names = {all_params.Name};
    num_params_total = length(all_params);
    
    % --- 1. Identify Fitted Parameters ---
    % A parameter is considered "fitted" if:
    % A) It was included in at least one batch (or 'ALL' was used).
    % B) Its margin factor allowed it to be Free (margin != 1).

    % A. Check Batch Strategy
    included_in_batch = false(1, num_params_total);
    for i = 1:length(batch_strategy)
        batch_params = batch_strategy{i}.params_to_fit;
        
        if any(strcmpi(batch_params, 'ALL'))
            included_in_batch(:) = true;
            break; % 'ALL' covers everything
        else
            % Find indices of parameters named in this batch
            % strcmpi compares strings case-insensitively
            is_in_this_batch = ismember(lower(all_param_names), lower(batch_params));
            included_in_batch = included_in_batch | is_in_this_batch;
        end
    end

    % B. Check Margin Factors
    % If margin is 1, setParameterBounds forces Free = false.
    allowed_by_margins = true(1, num_params_total);
    mf = settings.margin_factors;
    
    if isscalar(mf)
        if mf == 1, allowed_by_margins(:) = false; end
    elseif size(mf, 2) == 1 % Vector [N x 1]
        allowed_by_margins = (mf ~= 1)';
    elseif size(mf, 2) == 2 % Matrix [N x 2]
        % Fixed if BOTH min and max factors are 1
        is_fixed = (mf(:,1) == 1) & (mf(:,2) == 1);
        allowed_by_margins = ~is_fixed';
    end
    
    % Combine conditions
    fitted_indices = find(included_in_batch & allowed_by_margins);
    fitted_names = all_param_names(fitted_indices);

    if isempty(fitted_indices)
        warning('No fitted parameters found based on batch strategy and margins. Skipping robustness plot.');
        return;
    end
    
    num_params_fitted = length(fitted_indices);
    num_runs = length(all_fits);
    
    % --- 2. Extract and Filter Values ---
    % Matrix dimensions: [Number of Runs] x [Number of Fitted Parameters]
    param_values_matrix = zeros(num_runs, num_params_fitted);
    
    for k = 1:num_runs
        vals = [all_fits(k).sys.Structure.Parameters.Value];
        param_values_matrix(k, :) = vals(fitted_indices);
    end
    
    % --- 3. Calculate Stats ---
    robustness_sd_fitted = std(param_values_matrix, 0, 1);
    
    % Map back to full parameter vector for printing (fill non-fitted with 0)
    full_robustness_sd = zeros(1, num_params_total);
    full_robustness_sd(fitted_indices) = robustness_sd_fitted;
    
    % --- 4. Normalization ---
    baseline_vals = [all_params(fitted_indices).Value];
    baseline_vals(baseline_vals == 0) = 1; % Prevent divide by 0
    
    normalized_matrix = param_values_matrix ./ baseline_vals;
    
    % --- Visualization ---
    figure('Name', 'Robustness Analysis');
    hold on;
    
    % A. Plot "Spaghetti" lines (The Ensemble)
    %plot(1:num_params_fitted, normalized_matrix', 'Color', [0.7 0.7 0.7 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    % B. Plot Strip Chart (Jittered Scatter)
    jitter_width = 0.2; % Total width of the strip
    
    % Generate distinct colors for each fitted parameter using the 'lines' colormap
    param_colors = lines(num_params_fitted);

    for p = 1:num_params_fitted
        % Extract values for this parameter across all runs
        y_vals = normalized_matrix(:, p);
        
        % Generate random x-offsets (jitter) centered around the integer index p
        % (rand - 0.5) gives [-0.5, 0.5]. Multiply by width to scale.
        x_vals = p + (rand(num_runs, 1) - 0.5) * jitter_width;
        
        % Use the specific color for this parameter (p)
        scatter(x_vals, y_vals, 50, param_colors(p,:), 'filled', 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
    end
    
    % C. Plot Box and Whiskers
    % We use a dark grey [0.2 0.2 0.2] for the box outlines to frame the colored dots nicely
    h = boxplot(normalized_matrix, 'Positions', 1:num_params_fitted, 'Colors', [0.2 0.2 0.2], 'Symbol', 'r+', 'Widths', 0.5);
    %h = boxplot(normalized_matrix, 'Positions', 1:num_params_fitted, 'Colors', [0.2 0.2 0.2], 'Symbol', '', 'Widths', 0.5, 'Whisker', Inf);% Don't remove outliers
    set(h, 'LineWidth', 2);
    % D. Plot Baseline
    %yline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Standard Fit Baseline');
    
    % E. Plot Mean of Ensemble
    %mean_normalized = mean(normalized_matrix, 1);
    %plot(1:num_params_fitted, mean_normalized, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'Mean of Ensemble');
    
    % Formatting
    ylabel('Normalized Parameter Value');
    %title(sprintf('Parameter Convergence Robustness (%d Runs, Fitted Params Only)', num_runs));
    
    xlim([0.5, num_params_fitted + 0.5]);
    xticks(1:num_params_fitted);
    xticklabels(fitted_names);
    xtickangle(45);
    
    grid on;
    %legend('show', 'Location', 'best');
    
    % Adjust Y-limits
    ymin = min(normalized_matrix(:));
    ymax = max(normalized_matrix(:));
    range = ymax - ymin;
    if range == 0, range = 0.1; end
    %ylim([ymin - range*0.1, ymax + range*0.1]);
    ylim([0.7, 1.3]);
    
    hold off;
    
    % --- Report ---
    disp('Printing Fit Info using Robustness Standard Deviations (Sensitivity to Initial Guess):');
    printGreyFitInfo(sys_baseline, 0, full_robustness_sd);

end