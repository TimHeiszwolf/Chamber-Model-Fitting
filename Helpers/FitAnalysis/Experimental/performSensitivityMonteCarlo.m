function performSensitivityMonteCarlo(sys_grey, num_simulations)
% performSensitivityMonteCarlo performs a Monte Carlo sensitivity analysis.
% This function assesses the impact of parameter uncertainty on the model's
% step response. It runs multiple simulations with parameters sampled from
% normal distributions defined by their estimated values and standard
% deviations.
%
%   Inputs:
%       sys_grey        - The identified grey-box model (idgrey object).
%       num_simulations - The number of Monte Carlo simulations to run.

    % --- Get original parameters and their uncertainties ---
    [pvec, pvec_sd] = getpvec(sys_grey);
    param_struct = sys_grey.Structure.Parameters;
    num_params = length(pvec);
    
    % --- Get step response of the original system for comparison ---
    % Let the step function determine the time vector automatically
    [y_orig, t_sim] = step(sys_grey);
    num_outputs = size(y_orig, 2);
    
    % --- Initialize storage for all simulation results ---
    % This 3D matrix will store all step responses.
    % Dimensions: [time_points x outputs x simulations]
    all_step_responses = zeros(length(t_sim), num_outputs, num_simulations);

    fprintf('\n--- Starting Monte Carlo Sensitivity Analysis (%d simulations) ---\n', num_simulations);
    
    % --- Main Monte Carlo simulation loop ---
    for i = 1:num_simulations
        % Generate a new parameter set by sampling from a normal distribution
        % for each parameter.
        pvec_perturbed = normrnd(pvec, pvec_sd);

        % Enforce the parameter bounds defined in the model structure.
        % This prevents the sampler from choosing physically impossible values.
        for j = 1:num_params
            min_bound = param_struct(j).Minimum;
            max_bound = param_struct(j).Maximum;
            if pvec_perturbed(j) < min_bound
                pvec_perturbed(j) = min_bound;
            elseif pvec_perturbed(j) > max_bound
                pvec_perturbed(j) = max_bound;
            end
        end

        % Create a temporary model and set its parameters to the new
        % perturbed values.
        sys_perturbed = sys_grey;
        setpvec(sys_perturbed, pvec_perturbed);

        % Calculate the step response for the perturbed system using the
        % pre-determined time vector.
        y_perturbed = step(sys_perturbed, t_sim);
        
        % Store the resulting step response.
        all_step_responses(:, :, i) = y_perturbed;

        % Display progress to the user.
        if mod(i, 100) == 0
            fprintf('Completed %d simulations...\n', i);
        end
    end
    fprintf('All simulations completed.\n');

    % --- Analyze and Plot the Results ---
    % Calculate the mean and confidence intervals across all simulations.
    % Using percentiles is a robust way to define the confidence interval.
    mean_response = mean(all_step_responses, 3);
    lower_bound = prctile(all_step_responses, 2.5, 3); % 2.5th percentile
    upper_bound = prctile(all_step_responses, 97.5, 3); % 97.5th percentile

    % --- Create the plot ---
    figure('Name', 'Monte Carlo Sensitivity Analysis: Step Response');
    output_names = sys_grey.OutputName;
    
    for k = 1:num_outputs
        subplot(num_outputs, 1, k);
        hold on;

        % Plot the 95% confidence interval as a shaded area.
        fill([t_sim; flipud(t_sim)], [lower_bound(:, k); flipud(upper_bound(:, k))], ...
             [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

        % Plot the mean response from all the simulations.
        plot(t_sim, mean_response(:, k), 'b-', 'LineWidth', 1.5);

        % Plot the original, unperturbed system's response for reference.
        plot(t_sim, y_orig(:, k), 'r--', 'LineWidth', 1.5);
        
        hold off;
        grid on;
        title(['Step Response: ' output_names{k}]);
        xlabel('Time (s)');
        ylabel('Amplitude');
        if k == 1
             legend('95% Confidence Interval', 'Mean Perturbed Response', 'Original Model Response', 'Location', 'best');
        end
    end
    
    sgtitle(sprintf('Monte Carlo Analysis with %d Simulations', num_simulations));
end