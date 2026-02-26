function [training_data, validation_data] = splitData(data, validation_split_ratio)
%SPLIT_IDDATA Splits an iddata object into training and validation sets.
%
%   Syntax:
%       [training_data, validation_data] = split_iddata(data, validation_split_ratio)
%
%   Description:
%       This function takes an iddata object and splits it into two parts:
%       a training set and a validation set, based on the specified ratio.
%       The split is performed along the time axis (samples).
%
%   Inputs:
%       data - The input iddata object to be split.
%
%       validation_split_ratio - A scalar value between 0 and 1 that 
%       defines the fraction of the data to be used for validation.
%           - If 0, all data is used for training.
%           - If 1, all data is used for validation.
%           - If 0.3, 30% of the data is used for validation and 70% for training.
%
%   Outputs:
%       training_data   - An iddata object containing the training data.
%                         Will be empty if validation_split_ratio is 1.
%
%       validation_data - An iddata object containing the validation data.
%                         Will be empty if validation_split_ratio is 0.
    arguments
        data {mustBeA(data, 'iddata')}
        validation_split_ratio = 0.5
    end

    % Determine the number of experiments
    % size(data) returns [Samples, Outputs, Inputs, Experiments]
    num_experiments = size(data, 4);
    
    training_sets = {};
    validation_sets = {};
    
    % Loop through each experiment by index
    for k = 1:num_experiments
        % Extract single experiment using index
        exp_data = getexp(data, k);
        num_samples = size(exp_data, 1);
        
        % Calculate split point for this specific experiment
        % floor ensures we get a whole number index
        split_point = floor(num_samples * (1 - validation_split_ratio));
        
        % Logic to handle slicing
        if split_point >= num_samples
            % Case: 100% Training (ratio 0)
            train_part = exp_data;
            val_part = []; 
        elseif split_point <= 0
            % Case: 100% Validation (ratio 1)
            train_part = [];
            val_part = exp_data;
        else
            % Case: Normal Split
            % Note: For single exp, exp_data(1:split_point) works standardly
            train_part = exp_data(1:split_point);
            val_part = exp_data(split_point+1:end);
        end
        
        % Store parts if they aren't empty
        if ~isempty(train_part)
            training_sets{end+1} = train_part; 
        end
        
        if ~isempty(val_part)
            validation_sets{end+1} = val_part; 
        end
    end
    
    % Merge results back into combined iddata objects
    if isempty(training_sets)
        training_data = [];
    else
        training_data = merge(training_sets{:});
    end
    
    if isempty(validation_sets)
        validation_data = [];
    else
        validation_data = merge(validation_sets{:});
    end

end