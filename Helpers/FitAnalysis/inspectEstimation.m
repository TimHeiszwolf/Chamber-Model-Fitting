function inspectEstimation(sys_fitted, data, settings, fit_type_name)
    % inspectEstimation Visualizes and compares system identification results.
    %
    %   inspectEstimation(sys_fitted, data, settings, fit_type_name)
    %
    %   Inputs:
    %       sys_fitted    - A single model object, or a cell array of models.
    %       data          - An iddata object for comparison.
    %       settings      - A struct with a function handle to get matrices.
    %       fit_type_name - A string or cell array of strings matching the 
    %                       number of systems in sys_fitted.

    % --- Input Formatting ---
    % Ensure sys_fitted is a cell array for consistent looping
    if ~iscell(sys_fitted)
        systems_list = {sys_fitted};
    else
        systems_list = sys_fitted;
    end

    % Ensure fit_type_name is a cell array
    if ~iscell(fit_type_name)
        titles_list = {fit_type_name};
    else
        titles_list = fit_type_name;
    end

    % Verify dimensions match
    if length(systems_list) ~= length(titles_list)
        error('The number of provided fit names must match the number of systems.');
    end

    data.Name = 'Data';% This ensures the legend says "Data" instead of long strings like "Validation Data: ...".
    % --- Iterate through systems and generate plots ---
    for i = 1:length(systems_list)
        current_sys = systems_list{i};
        current_title = titles_list{i};
        
        % 1. Time-Domain Comparison
        figure('Name', sprintf('Time-Domain %s', current_title));
        opt_compare = compareOptions('InitialCondition', 'e');% Configure options (e.g., estimated initial conditions)
        compare(data, current_sys, Inf, opt_compare);
        styleTimeDomainComparison(data, settings.validation_split, current_title);

        % 2. Bode Plots (Using existing external function)
        inspectFrequencyResponse(data, settings.LPM_data, current_sys, current_title);

        % 3. Residual Analysis Plot
        figure('Name', sprintf('Residual %s', current_title));
        resid(data, current_sys);
        %title(sprintf('Residual Analysis: %s', current_title), 'Interpreter', 'none');%Title could be commented out for if you want it easier use it in the report.
        grid on;
    end
end