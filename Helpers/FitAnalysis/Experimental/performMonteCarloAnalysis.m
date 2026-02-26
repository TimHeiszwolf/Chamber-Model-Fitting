function performMonteCarloAnalysis(sys_grey, data, settings, num_simulations)
% performMonteCarloAnalysis Performs Monte Carlo sensitivity analysis on a grey-box model.
%   This function simulates the model response multiple times using parameter sets
%   sampled from normal distributions defined by the estimated parameters and
%   their standard deviations. It then plots the mean response, a 95%
%   confidence interval, and the original measured data for comparison.
%
%   Inputs:
%       sys_grey        - The fitted idgrey model object.
%       data            - The iddata object containing the input/output data.
%       settings        - The main settings struct, containing the chamber_model object.
%       num_simulations - (Optional) The number of Monte Carlo iterations. Defaults to 200.

    % ---- Argument Validation ----
    if nargin < 4
        num_simulations = 200; % Default number of simulations
    end

    fprintf('--- Starting Monte Carlo Sensitivity Analysis with %d simulations ---\n', num_simulations);

    % ---- 1. Extract Model Information ----
    [pvec, pvec_sd] = getpvec(sys_grey);
    params_struct = sys_grey.Structure.Parameters;

    num_params = length(pvec);
    num_outputs = size(data.y, 2);
    time_vector = data.SamplingInstants;
    num_time_steps = length(time_vector);

    % ---- 2. Run Simulations ----
    % Pre-allocate a matrix to store all output trajectories
    all_outputs = zeros(num_time_steps, num_outputs, num_simulations);
    
    % Setup progress bar
    h_waitbar = waitbar(0, 'Running Monte Carlo Simulations...');

    for i = 1:num_simulations
        % Update progress bar
        waitbar(i / num_simulations, h_waitbar);

        % ---- 2a. Sample New Parameter Set ----
        % This loop ensures that sampled parameters respect the min/max bounds
        sampled_pvec = zeros(size(pvec));
        for j = 1:num_params
            is_valid_sample = false;
            while ~is_valid_sample
                % Sample from a normal distribution
                sample = pvec(j) + pvec_sd(j) * randn();
                
                % Check if the sample is within the defined bounds
                min_bound = params_struct(j).Minimum;
                max_bound = params_struct(j).Maximum;
                if sample >= min_bound && sample <= max_bound
                    sampled_pvec(j) = sample;
                    is_valid_sample = true;
                end
            end
        end

        % ---- 2b. Create and Simulate System ----
        try
            % Convert the vector of sampled parameters to a cell array
            params_cell = num2cell(sampled_pvec);
            
            % Get matrices for the sampled parameters using the model's method
            [A_s, B_s, C_s, D_s] = settings.chamber_model.getMatrices(params_cell{:});
            
            % Create a temporary continuous-time state-space system
            sys_temp = idss(A_s, B_s, C_s, D_s, 'Ts', 0);
            
            % Simulate the response using 'zero' initial condition for speed and stability
            opt = simOptions('InitialCondition', 'zero'); 
            y_sim = sim(sys_temp, data, opt);

            % Store the results
            all_outputs(:, :, i) = y_sim.y;

        catch ME
            warning('Simulation failed for one parameter set: %s. Skipping this iteration.', ME.message);
            % Store NaNs for failed simulations to ignore them later
            all_outputs(:, :, i) = NaN;
        end
    end
    
    close(h_waitbar); % Close the progress bar

    % ---- 3. Analyze and Plot Results ----
    figure;
    sgtitle(sprintf('Monte Carlo Analysis (%d Simulations)', num_simulations), 'FontSize', 14, 'FontWeight', 'bold');
    
    output_names = data.OutputName;
    time_for_fill = [time_vector; flipud(time_vector)];

    for k = 1:num_outputs
        subplot(num_outputs, 1, k);
        hold on;

        % Get the outputs for the current channel, ignoring any failed (NaN) runs
        output_k = squeeze(all_outputs(:, k, :));
        output_k = output_k(:, all(~isnan(output_k), 1));

        if isempty(output_k)
            title(sprintf('%s - All simulations failed', output_names{k}));
            grid on;
            continue;
        end

        % Calculate mean and 95% confidence interval
        mean_response = mean(output_k, 2);
        lower_bound = prctile(output_k, 2.5, 2);
        upper_bound = prctile(output_k, 97.5, 2);

        % Plot 95% confidence interval as a shaded area
        fill(time_for_fill, [upper_bound; flipud(lower_bound)], [0.8 0.8 1], 'EdgeColor', 'none', 'DisplayName', '95% Confidence Interval');

        % Plot measured data and mean response
        plot(time_vector, data.y(:, k), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Measured Data');
        plot(time_vector, mean_response, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean Simulated Response');
        
        hold off;
        grid on;
        title(output_names{k});
        xlabel('Time (s)');
        ylabel('Amplitude'); % Using generic 'Amplitude' as it might be normalized
        legend('Location', 'best');
        xlim([time_vector(1), time_vector(end)]);
    end
    
    fprintf('--- Monte Carlo Analysis Complete ---\n');
end