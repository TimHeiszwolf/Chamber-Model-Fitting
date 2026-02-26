function printGreyFitInfo(sys_grey, verbose, external_sd)
% printGreyFitInfo Displays parameters and fit info from a grey-box model.
%
%   Inputs:
%       sys_grey    - The grey-box model object.
%       verbose     - (Optional) Verbosity level (0, 1, 2, 3).
%       external_sd - (Optional) A vector of standard deviations to use 
%                     INSTEAD of the ones stored in sys_grey. Use this when
%                     printing accumulated batch uncertainties.

    % --- Check for optional inputs ---
    if nargin < 2
        verbose = 0;
    end
    
    % Check if external_sd is provided as the 3rd argument
    use_external_sd = (nargin >= 3 && ~isempty(external_sd));

    correlation_threshold = 0.75; 

    % --- Get parameter vectors and their standard deviations ---
    [pvec, pvec_sd_internal] = getpvec(sys_grey);
    
    % Decide which SD to use
    if use_external_sd
        pvec_sd = external_sd;
    else
        pvec_sd = pvec_sd_internal;
    end

    % --- Get parameter names from the model structure ---
    param_names = {sys_grey.Structure.Parameters.Name};
    
    % --- Print the formatted table header ---
    fprintf('\n--- Estimated Parameters, Uncertainties, and Bounds ---\n');
    fprintf('%-25s | %-12s | %-12s | %-15s | %-15s\n', 'Parameter', 'Value', 'Std.', 'Rel. Unc. (%)', 'Bound Pos. (%)');
    fprintf('%s\n', repmat('-', 1, 91)); 

    % --- Loop through each parameter and print its info ---
    for i = 1:length(pvec)
        % --- Handle case where standard deviations are not computed ---
        % Note: If using external_sd, we check if the specific element is > 0
        if ~isempty(pvec_sd) && pvec_sd(i) > 0
            relative_uncertainty = abs(pvec_sd(i) / pvec(i)) * 100;
            relative_uncertainty_str = sprintf('%.2f%%', relative_uncertainty);
            sd_val_str = sprintf('%.4e', pvec_sd(i));
        else
            relative_uncertainty_str = 'N/A';
            sd_val_str = 'N/A';
        end

        % --- Calculate position within bounds ---
        min_bound = sys_grey.Structure.Parameters(i).Minimum;
        max_bound = sys_grey.Structure.Parameters(i).Maximum;
        range = max_bound - min_bound;

        if range > 0 && isfinite(range)
            bound_position_percent = ((pvec(i) - min_bound) / range) * 100;
            bound_str = sprintf('%.2f%%', bound_position_percent);
        else
            bound_str = 'N/A';
        end

        % --- Print row ---
        fprintf('%-25s | %-12.4e | %-12s | %-15s | %-15s\n', ...
            param_names{i}, pvec(i), sd_val_str, relative_uncertainty_str, bound_str);
    end
    fprintf('\n');

    % --- Verbose Reporting Section ---
    if verbose >= 1
        fit_report = sys_grey.Report;
 
        % --- Correlation Matrix Analysis (verbose >= 2) ---
        if verbose >= 2 && ~use_external_sd && ~isempty(pvec_sd_internal)
            % We only do correlation analysis if we are looking at the INTERNAL
            % fit report, because correlation matrices in sys.Report are 
            % tied to the specific batch, not the cumulative SDs.
            
            all_params = sys_grey.Structure.Parameters;
            is_free = [all_params.Free];
            free_param_names = {all_params(is_free).Name};
            
            printCorrelationAnalysis(fit_report, free_param_names, correlation_threshold, verbose);
            
        elseif verbose >= 2
            fprintf('--- Correlation Analysis Skipped ---\n');
            fprintf('Cannot analyze correlation on accumulated batch SDs (requires full covariance matrix).\n');
        end
 
        % --- Display Fit and Termination Info (verbose >= 1) ---
        fprintf('\n--- Fit and Termination Report ---\n');
        if isfield(fit_report, 'fit')
            disp(fit_report.fit);
        end
        if isfield(fit_report, 'Termination')
            disp(fit_report.Termination);
        end
    end
end