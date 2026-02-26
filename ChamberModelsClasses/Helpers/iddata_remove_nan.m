function data_out = iddata_remove_nan(data_in)
    % This function checks an iddata object for NaN values,
    % interpolates them if found, and reports the count.
    
    % 1. Count NaN values in Input and Output
    % We use ismissing() as it's the modern standard and 'all' sums
    % over all dimensions, perfect for MIMO data too.
    nan_count_input = sum(ismissing(data_in.InputData), 'all');
    nan_count_output = sum(ismissing(data_in.OutputData), 'all');
    
    total_nans = nan_count_input + nan_count_output;
    
    % 2. Check if any NaNs were found
    if total_nans > 0
        % 3. If NaNs exist, correct them using 'interpolate'
        %fprintf('Data contains NaN values. Correcting...\n');
        data_out = misdata(data_in);
        
        % 4. Issue a warning with the details
        warning('IDDATA:NaNsCorrected', ...
                'Corrected %d total NaN value(s): %d in InputData and %d in OutputData.', ...
                total_nans, nan_count_input, nan_count_output);
    else
        % 5. If no NaNs, just return the original object
        %fprintf('Data is clean. No NaN values found.\n');
        data_out = data_in;
    end
end