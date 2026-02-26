function printCorrelationAnalysis(fit_report, free_param_names, correlation_threshold, verbose)
% printCorrelationAnalysis Calculates and displays parameter correlation info.
% This helper function takes an estimation report and displays a correlation
% analysis for the free parameters, including the full matrix and a list of
% highly correlated pairs based on a threshold.
%
%   Inputs:
%       fit_report              - The model's estimation report object.
%       free_param_names        - Cell array of strings with free parameter names.
%       correlation_threshold   - Numeric threshold for reporting high correlation.
%       verbose                 - Verbosity level (>=3 displays full matrix).

    fprintf('\n--- Correlation Matrix of Free Parameters ---\n');
    free_par_cov = fit_report.Parameters.FreeParCovariance;

    % Check if correlation is possible
    if isempty(free_par_cov) || size(free_par_cov, 1) <= 1
        fprintf('Correlation matrix requires at least two free parameters.\n');
        return; % Exit the function
    end

    % Calculate the correlation matrix from the covariance matrix
    [~, correlation_matrix] = cov2corr(free_par_cov);

    % --- Display full correlation table (verbose >= 3) ---
    if verbose >= 3
        % Round the matrix for a cleaner display
        correlation_matrix_rounded = round(correlation_matrix, 2, 'significant');
        % Create and display a labeled table
        correlation_table = array2table(correlation_matrix_rounded, 'RowNames', free_param_names, 'VariableNames', free_param_names);
        disp(correlation_table);
    end

    % --- Check for and report high correlations ---
    highly_correlated_pairs = {}; % Use a cell array to store findings

    num_params = length(free_param_names);
    for i = 1:num_params
        for j = (i + 1):num_params % Upper triangle to avoid duplicates/diagonal
            corr_value = correlation_matrix(i, j);
            if abs(corr_value) > correlation_threshold
                highly_correlated_pairs{end+1} = {free_param_names{i}, free_param_names{j}, corr_value};
            end
        end
    end

    % Print the results
    if ~isempty(highly_correlated_pairs)
        fprintf('The following parameter pairs are highly correlated (Threshold > %.2f):\n', correlation_threshold);
        for k = 1:length(highly_correlated_pairs)
            pair = highly_correlated_pairs{k};
            fprintf('  - "%s" and "%s" (Correlation: %.2f)\n', pair{1}, pair{2}, pair{3});
        end
        fprintf('\n');
    end
end