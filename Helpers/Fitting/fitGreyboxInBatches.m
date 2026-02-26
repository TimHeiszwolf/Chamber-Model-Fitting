function [final_model, full_history] = fitGreyboxInBatches(sys_init, batch_strategy, settings, validation_data)
% fitGreyboxInBatches Performs iterative, batched fitting of an idgrey model.
%
%   Inputs:
%       sys_init        - The initial idgrey model.
%       batch_strategy  - A cell array of structs defining the fitting strategy.
%                         Each struct must contain:
%                           .data: The iddata object for this batch.
%                           .params_to_fit: Cell array of param names to fit.
%                                             Use {'ALL'} to fit all free params.
%                           .iterations: Number of iterations for this batch.
%       settings        - The main settings struct (must contain
%                         .search_method). Can also contain .batch_cycles.
%       fit_history_settings - Struct with settings for fitGreyboxWithHistory
%                            (step_size, verbose, validation_data).
%
%   Outputs:
%       final_model   - The idgrey model after all batches and cycles.
%       full_history  - A cell array containing the history from each
%                       individual batch.

    arguments
        sys_init idgrey
        batch_strategy cell
        settings struct
        validation_data iddata = [] 
    end
    
    fprintf('\n--- Starting Batch Grey-Box Fitting ---\n');
    
    current_model = sys_init;
    full_history = {};
    history_idx = 1;
    
    % Store the original "globally free" status.
    globally_free_status = [sys_init.Structure.Parameters.Free];
    all_param_names = {sys_init.Structure.Parameters.Name};
    num_params = length(all_param_names);
    
    up_to_date_pvec_sd = zeros(num_params, 1);

    num_cycles = 1;
    if isfield(settings, 'batch_cycles')
        num_cycles = settings.batch_cycles;
    end
    
    num_batches = length(batch_strategy);
    
    for cycle = 1:num_cycles
        fprintf('\n \n--- Starting Cycle %d / %d ---\n', cycle, num_cycles);
        
        for batch_num = 1:num_batches
            batch_info = batch_strategy{batch_num};
            
            if isfield(batch_info, 'max_cycle') && batch_info.max_cycle <= cycle
                continue
            end

            fprintf('  --- Running Batch %d / %d (Cycle %d / %d) ---\n', batch_num, num_batches, cycle, num_cycles);
            
            % --- 1. Prepare the model for this batch ---
            sys_batch_input = current_model;
            params_to_fit_names = batch_info.params_to_fit;
            
            is_fit_all = any(strcmpi(params_to_fit_names, 'ALL'));
            
            fprintf('    Parameters to fit: ');
            for j = 1:num_params
                is_in_batch = is_fit_all || any(strcmpi(all_param_names{j}, params_to_fit_names));
                
                if is_in_batch
                    sys_batch_input.Structure.Parameters(j).Free = true;
                    fprintf('"%s" ', all_param_names{j});
                else
                    sys_batch_input.Structure.Parameters(j).Free = false;
                end
            end
            fprintf('\n');
            
            % --- 2. Prepare settings for this batch ---
            settings_for_batch = settings;
            settings_for_batch.max_iterations = batch_info.iterations;
            data_for_batch = batch_info.data;

            % --- 3. Run the batch ---
            [fitted_batch_model, history] = fitGreyboxWithHistory(...
                data_for_batch, ...
                sys_batch_input, ...
                settings_for_batch, ...
                settings.step_size, ...
                0, ...
                validation_data);
            
            % --- 4. Update Model and SD Tracking ---
            current_model = fitted_batch_model;
            
            % Get current batch SDs
            [~, batch_sd] = getpvec(current_model);
            
            if ~isempty(batch_sd)
                % Identify indices where we have a valid estimate (non-zero SD)
                % Parameters fixed in this batch will have SD = 0
                valid_sd_indices = batch_sd > 0;
                
                % Update only the valid ones
                up_to_date_pvec_sd(valid_sd_indices) = batch_sd(valid_sd_indices);
            end
            
            % Store in history
            history.up_to_date_pvec_sd = up_to_date_pvec_sd;
            full_history{history_idx} = struct('cycle', cycle, ...
                                               'batch_num', batch_num, ...
                                               'history', history);%, ...
                                               %'up_to_date_pvec_sd', up_to_date_pvec_sd);
            history_idx = history_idx + 1;

            if settings.verbose >= 0
                % Pass the cumulative SDs to the print function
                printGreyFitInfo(current_model, 0, up_to_date_pvec_sd);
            end
        end
        
        fprintf('--- Cycle %d / %d Complete ---\n\n', cycle, num_cycles);
    end
    
    final_model = current_model;

    plotBatchHistory(full_history, sys_init, settings.verbose)

    fprintf('--- Batch Grey-Box Fitting Finished ---\n');
    
end
