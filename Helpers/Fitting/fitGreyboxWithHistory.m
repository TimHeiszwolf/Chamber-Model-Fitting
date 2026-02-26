function [sys_final, history] = fitGreyboxWithHistory(training_data, sys_init, settings, step_size, verbose, validation_data)
    %fitGreyboxWithHistory Iteratively fits a grey-box model and logs history.
    %
    %   This function provides a wrapper around MATLAB's 'greyest' function to
    %   perform iterative model fitting in discrete steps. It is useful for
    %   long-running optimizations, allowing for:
    %     - Logging of parameter values and MSE (Training & Validation) at each step.
    %     - Real-time console output with progress and Estimated Time Remaining (ETR).
    %     - Generation of plots summarizing the fitting process.
    %
    %   Syntax:
    %     [sys_final, history] = fitGreyboxWithHistory(training_data, sys_init, settings)
    %     [sys_final, history] = fitGreyboxWithHistory(..., step_size, verbose, validation_data)
    %
    %   Outputs:
    %     sys_final: The final 'idgrey' model after all iterations.
    %     history:   A struct containing logs of:
    %                .fit:        Fit report from each step
    %                .train_mse:  MSE on training data over iterations
    %                .val_mse:    (Optional) MSE on validation data
    %                .parameters: Parameter values over iterations
    %                .fitting_time: Struct with .iterations and .time logs
    
    % --- Input Arguments ---
    arguments
        training_data   iddata      % The training data (iddata object)
        sys_init        idgrey      % The initial grey-box model (idgrey object)
        settings        struct      % Struct with settings, must contain:
                                    %   .max_iterations (total iterations to run)
                                    %   .search_method  (e.g., 'auto', 'lm', 'gn')
        % Optional inputs
        step_size       double = 10 % Iterations per loop step (default: 10)
        verbose         double = 1  % Verbosity level:
                                    %   0: Off
                                    %   1: Print status to console
                                    %   2: +Plot MSE/Loss
                                    %   3: +Plot Parameter evolution
                                    %   4: +Plot Fitting time
        validation_data iddata = [] % Optional validation data (iddata object)
    end
    
    % --- Initialization ---
    p_fit_order = 2; % Order of the polynomial fit for time estimation (ETR)
    
    sys_current = sys_init;

    % Configure options for the 'greyest' estimation function
    opt = greyestOptions(); 
    opt.display = 'off';
    opt.InitialState = 'estimate'; % 'estimate', 'backcast', 'backcast', 'model'
    opt.SearchMethod = settings.search_method; 
    if isfield(settings, 'output_weight')
        opt.OutputWeight = settings.output_weight;
    else
        opt.OutputWeight = [];% How much weight each output has. Can be "noise", can also be diag([1, 1, 20, 1]); (third channel has 20 times more weight). Or an empty array [] if you want it all equally.
    end
    opt.SearchOptions.MaxIter = step_size;
    %opt.SearchOptions.Tolerance = 1e-5;
    
    num_steps = ceil(settings.max_iterations / step_size); % Calculate the total number of loops (steps) needed
    current_iteration = 0;
    
    % --- History/Log Initialization ---
    history_fit = [];
    e_train = pe(training_data, sys_current); % Calculate initial training MSE
    
    % Handle Multi-Experiment (Cell Array) OutputData
    err_train_data = e_train.OutputData;
    if iscell(err_train_data)
        % Force vertical concatenation (stacking time series on top of each other)
        % This works even if experiments have different lengths.
        err_train_data = vertcat(err_train_data{:});
    end
    history_train_mse = [mean(err_train_data(:).^2)];
    
    if ~isempty(validation_data) % has_val_data
        e_val = pe(validation_data, sys_current); 
        
        err_val_data = e_val.OutputData;
        if iscell(err_val_data)
            err_val_data = vertcat(err_val_data{:});
        end
        history_val_mse = [mean(err_val_data(:).^2)];
    end
    
    history_parameters = [getpvec(sys_init)]; % Add the initial parameters
    history_time.iterations = [0];
    history_time.time = [0];
    
    tic; % Start the timer before the loop
    
    % --- Main Fitting Loop ---
    for i = 1:num_steps
        sys_current = greyest(training_data, sys_current, opt); % Run the 'greyest' estimation for 'step_size' iterations
        current_iteration = current_iteration + step_size;
        
        % --- Log History ---
        fit_report = sys_current.Report;
        history_fit = [history_fit fit_report.fit];
        
        estimated_params = getpvec(sys_current);
        history_parameters = [history_parameters estimated_params];
        
        e_train = pe(training_data, sys_current);
        
        % Handle Multi-Experiment OutputData
        err_train_data = e_train.OutputData;
        if iscell(err_train_data)
            err_train_data = vertcat(err_train_data{:});
        end
        current_train_mse = mean(err_train_data(:).^2);
        history_train_mse = [history_train_mse current_train_mse];
        
        % --- Timing and ETA Calculations ---
        time_elapsed_sec = toc;
        history_time.iterations = [history_time.iterations current_iteration];
        history_time.time = [history_time.time time_elapsed_sec];
        
        avg_iters_per_min = current_iteration / time_elapsed_sec * 60;
        
        if i >= 5 % After 5 steps, we have enough data to fit a curve and do average.
            [p_fit, S_fit] = polyfit(history_time.iterations, history_time.time, p_fit_order);
            [estimated_end_time, ~] = polyval(p_fit, settings.max_iterations, S_fit);
            
            etr_min = (estimated_end_time - time_elapsed_sec) / 60;
            iters_per_min = step_size * 5 / (history_time.time(end) - history_time.time(end - 5)) * 60;
        else % Else use simple linear estimation
            iters_per_min = avg_iters_per_min;
            etr_min = (settings.max_iterations - current_iteration) / iters_per_min;
        end
        
        % --- Console Output ---
        iter_char_width = length(num2str(settings.max_iterations));
        
        if ~isempty(validation_data) % has_val_data
            % Calculate and log validation MSE
            e_val = pe(validation_data, sys_current);
            
            err_val_data = e_val.OutputData;
            if iscell(err_val_data)
                err_val_data = vertcat(err_val_data{:});
            end
            current_val_mse = mean(err_val_data(:).^2);
            history_val_mse = [history_val_mse current_val_mse];
            
            if verbose >= 1
                fprintf('Iter %*d/%*d | Iter/min: %.1f | ETA (min): %.1f | Train MSE: %.4e | Val MSE: %.4e | Loss: %.4e\n', iter_char_width, current_iteration, iter_char_width, settings.max_iterations, iters_per_min, etr_min, current_train_mse, current_val_mse, fit_report.fit.LossFcn);
            end
        elseif verbose >= 1
            fprintf('Iter %*d/%*d | Iter/min: %.1f | ETA (min): %.1f | Train MSE: %.4e | Loss: %.4e\n', iter_char_width, current_iteration, iter_char_width, settings.max_iterations, iters_per_min, etr_min, current_train_mse, fit_report.fit.LossFcn);
        end
        
        % --- Early Stopping Check ---
        if strcmp(fit_report.Termination.WhyStop, 'Near (local) minimum, (norm(g) < tol).') 
            break
        end
    end
    
    % --- Finalization ---
    sys_final = sys_current;
    
    % --- Plotting (based on verbosity level) ---
    if verbose >= 2
        % Plot 1.1: MSE
        Loss_values = [history_fit.LossFcn];
        num_actual_steps = size(history_fit, p_fit_order);
        iterations = (1:num_actual_steps) * step_size;
        
        figure
        subplot(2, 1, 1)
        semilogy([0 iterations], history_train_mse) 
        title('Calculated MSE versus Iterations')
        xlabel('Iterations');
        ylabel('MSE');
        grid on;
        
        if ~isempty(validation_data) % has_val_data
            hold on
            semilogy([0 iterations], history_val_mse, '--')
            hold off
            legend('Training MSE', 'Validation MSE')
        else
            legend('Training MSE')
        end
        
        % Plot 1.2: Loss
        subplot(2, 1, 2)
        semilogy(iterations, Loss_values)
        title('Loss function versus Iterations')
        xlabel('Iterations');
        ylabel('Loss Value');
        grid on;
    end
    
    if verbose >= 3
        % Plot 2: Normalized Parameter Evolution
        params_init = getpvec(sys_init); 
        params_init(params_init == 0) = eps; 
        normalized_params = history_parameters ./ params_init; 
        param_names = {sys_final.Structure.Parameters.Name}; 
        
        figure
        semilogy([0 iterations], normalized_params)
        title('Normalized Parameter Evolution (p / p_{init})')
        xlabel('Iterations');
        ylabel('Normalized Value (p/p_0)');
        grid on;
        legend(param_names, 'Interpreter', 'none'); 
    end
    
    if verbose >= 4
        % Plot 3: Fitting Time vs. Iterations
        p_fit = polyfit(history_time.iterations, history_time.time, 2); 
        y_fit = polyval(p_fit, history_time.iterations);
        
        figure
        hold on;
        plot(history_time.iterations, history_time.time)
        plot(history_time.iterations, y_fit, 'r--')
        hold off;
        title('Fitting time versus iterations')
        xlabel('Iterations');
        ylabel('Time (s)');
        grid on;
        legend('Data', 'Polynomial fit')
    end
    
    % --- Assemble Final Output Struct ---
    history.fit = history_fit;
    history.train_mse = history_train_mse;
    if ~isempty(validation_data)
        history.val_mse = history_val_mse;
    end
    history.parameters = history_parameters;
    history.fitting_time = history_time;
    
end