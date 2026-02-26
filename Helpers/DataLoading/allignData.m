function all_outputs_fix = allignData(all_inputs, normalize_data, detrend_data, filter_data, smooth_data, new_sample_time, minimum_start_time, maximum_end_time)
% Synchronizes multiple time-series datasets to a common time vector.
    %
    % Syntax:
    %   all_outputs = allignData(all_inputs)
    %   all_outputs = allignData(all_inputs, normalize_data, detrend_data, ...)
    %
    % Description:
    %   This function takes a cell array of structures (each with .time and .data)
    %   and aligns them to a shared time vector based on their overlap. It handles
    %   resampling, outlier filtering, smoothing, and normalization.
    %
    % Inputs:
    %   all_inputs         - Cell array of structs {struct1, struct2, ...}. 
    %                        Each struct must have .time and .data fields.
    %                        Supports nested cells for division: {numerator_struct, denominator_struct}.
    %   normalize_data     - 0: None, 1: Max to 1, 2: Min-Max (0 to 1). 
    %                        Also accepts a [offset, scale] matrix per input.
    %   detrend_data       - Struct with .setting (bool) and .f0 (cutoff freq).
    %   filter_data        - Boolean; if true, removes outliers using quartiles.
    %   smooth_data        - Window size in seconds for movmean; false to skip.
    %   new_sample_time    - "Highest" (lowest freq), "Lowest" (highest freq), or numeric value.
    %   minimum_start_time - Hard limit for the start of the aligned vector.
    %   maximum_end_time   - Hard limit for the end of the aligned vector.
    %
    % Outputs:
    %   all_outputs_fix    - Cell array of synchronized structs with identical .time vectors.

    arguments% Use this to specify inputs and set defaults.
        all_inputs% A set of structures containing a .time and .data with time and values. Can contain nested cells like {numerator, denominator} for division.
        normalize_data = 1;% false/0 is no normalisation, true/1 is normalisation of maximum to one, 2 is normalisation of minimum to zero and maximum to one.
        detrend_data = struct('setting', false, 'f0', 1)
        filter_data = true;
        smooth_data = false;% If false, don't smooth, else give value in (seconds) of the window size.
        new_sample_time = "Highest"% If numeric is the new sample time, depending on the string (either "Highest" or "Lowest") the highest/lowest sample time in the data will be picked.
        minimum_start_time = -Inf
        maximum_end_time = +Inf
    end

    %% Pre-process inputs for division operations
    % This section checks for nested cell arrays which signify a division operation.
    % It aligns the numerator and denominator, divides them, and replaces the
    % nested cell with the resulting struct for further processing.
    processed_inputs = {}; % Use a temporary cell array to store the processed inputs
    for i = 1:length(all_inputs)
        if iscell(all_inputs{i})
            % This block handles a division request, e.g., {numerator, denominator}.
            
            % Validate that the nested cell is correctly formatted for division.
            if numel(all_inputs{i}) ~= 2
                error('Nested cell for division at index %d must contain exactly two elements: {numerator, denominator}.', i);
            end
            numerator_struct = all_inputs{i}{1};
            denominator_struct = all_inputs{i}{2};
            if ~isstruct(numerator_struct) || ~isfield(numerator_struct, 'time') || ~isfield(numerator_struct, 'data') || ...
               ~isstruct(denominator_struct) || ~isfield(denominator_struct, 'time') || ~isfield(denominator_struct, 'data')
                error(['Division cell at index ', num2str(i), ' must contain two valid structs, each with .time and .data fields.']);
            end

            % Recursively call allignData to align just the numerator and denominator.
            % This is done without normalization/filtering and at the highest resolution ('Lowest' sample time)
            % to preserve data integrity before division.
            %aligned_pair = allignData({numerator_struct, denominator_struct}, 0, detrend_data, false, false, "Lowest", -Inf, +Inf);
            aligned_pair = allignData({numerator_struct, denominator_struct}, 0, struct('setting', false, 'f0', 1), false, false, "Lowest", -Inf, +Inf);
            
            % Perform the element-wise division.
            divided_struct = struct();
            divided_struct.time = aligned_pair{1}.time;
            
            %denominator_data = aligned_pair{2}.data;
            if any(aligned_pair{2}.data == 0)
                warning('Division by zero encountered at input index %d. Resulting values will be Inf.', i);
            end
            divided_struct.data = aligned_pair{1}.data ./ aligned_pair{2}.data;
            
            % Add the newly created divided struct to the list for main processing.
            processed_inputs{end+1} = divided_struct;
            
        else
            % If the input is a regular struct, add it directly to the list.
            processed_inputs{end+1} = all_inputs{i};
        end
    end
    
    % Replace the original inputs with the pre-processed ones.
    all_inputs = processed_inputs;
    num_structures = length(all_inputs);

    %% Validate structures within the cell array
    for k = 1:num_structures
        if ~isstruct(all_inputs{k}) || ~isfield(all_inputs{k}, 'time') || ~isfield(all_inputs{k}, 'data')
            error(['Input structure at index ', num2str(k), ' is invalid. It must have .time and .data fields.']);
        end
        if isempty(all_inputs{k}.time) || isempty(all_inputs{k}.data)
            error(['Input structure at index ', num2str(k), ' has empty .time or .data.']);
        elseif size(all_inputs{k}.time,1) ~= size(all_inputs{k}.data,1)
             error(['For input structure at index ', num2str(k), ', the number of rows in .data must match the length of .time.']);
        end
    end


    %% Optional detrending
    if detrend_data.setting
        for k = 1:num_structures
            current_struct = all_inputs{k};
            dt_array = diff(current_struct.time);
            delta_time = median(dt_array);
    
            yc = frf_cutsignal(current_struct.data, current_struct.time, minimum_start_time-0.01, detrend_data.f0, 1/delta_time, 1, maximum_end_time+0.01);%-0.01 added to prevent nans.
            %yc  = frf_cutsignal(data.y(:,i), data.SamplingInstants, settings.start_time, f_exc(1), fs, P, settings.end_time);
            yc.time = yc.time';
            yc.data = yc.data_detrend';
            all_inputs{k} = yc;
        end
    end

    %% Configuration
    delta_time_of_arrays = zeros(1, num_structures);% Always calculate original sample times to support anti-aliasing/smoothing
    for k = 1:num_structures
        dt_array = diff(all_inputs{k}.time);
        delta_time_of_arrays(k) = median(dt_array);
    end

    if ~isnumeric(new_sample_time)%isstring(new_sample_time)||ischar(new_sample_time)% Checks if the new_sample_time needs to be determined or if it already is a value.
        if strcmpi(new_sample_time, "Highest")
            new_sample_time = max(delta_time_of_arrays);
        elseif strcmpi(new_sample_time, "Lowest")
            new_sample_time = min(delta_time_of_arrays);
        else
            error("new_sample_time is not well defined.")
        end
    end

    % Iterate through the rest of the structures to find the common overlap
    overall_overlap_start_time = min(all_inputs{1}.time(:));
    overall_overlap_end_time = max(all_inputs{1}.time(:));
    for k = 1:num_structures
        current_struct = all_inputs{k};
        current_min_time = min(current_struct.time(:));
        current_max_time = max(current_struct.time(:));
        
        overall_overlap_start_time = max(overall_overlap_start_time, current_min_time);
        overall_overlap_end_time = min(overall_overlap_end_time, current_max_time);
        %disp(['After structure ', num2str(k), ' (Range: [', num2str(current_min_time), ', ', num2str(current_max_time), ']s) - Current Overlap: [', num2str(overall_overlap_start_time), ', ', num2str(overall_overlap_end_time), ']s']);
    end

    overall_overlap_start_time = max(overall_overlap_start_time, minimum_start_time);% Shorten the common time if overlap is larger than desired range.
    overall_overlap_end_time = min(overall_overlap_end_time, maximum_end_time);

    common_time_vector = (overall_overlap_start_time:new_sample_time:overall_overlap_end_time)';% Make the common time vector.

    % Adjust common_time_vector if the last point slightly exceeds overlap_end_time
    if ~isempty(common_time_vector) && common_time_vector(end) > overall_overlap_end_time && length(common_time_vector) > 1
        % If only one point and it's > end_time, or many points and last is > end_time
        if (common_time_vector(end) - overall_overlap_end_time) > new_sample_time * 0.5 % If significantly over
             common_time_vector(end) = [];
        else % if only slightly over, can adjust the last point, or remove. Removing is safer.
             common_time_vector(end) = []; % Simplest is to remove
        end
    end

    if isempty(common_time_vector) || ( (overall_overlap_end_time - overall_overlap_start_time) >= new_sample_time && length(common_time_vector) < 2 )
        error('Overlap duration is smaller than the new sample time.');
    end

    %% filter data for outliers
    if filter_data
        for k = 1:num_structures
                current_struct_original = all_inputs{k};
                
                time_mask = (current_struct_original.time >= (overall_overlap_start_time - 2*new_sample_time)) & (current_struct_original.time <= (overall_overlap_end_time + 2*new_sample_time));
                
                original_data = current_struct_original.data(time_mask);
                original_time = current_struct_original.time(time_mask);
    
                % Replace outliers with NaN instead of removing them to maintain time structure
                [clean_data, outlier_idx] = filloutliers(original_data, 'linear', 'quartiles');
    
                % Ensure time is unique for interpolation
                [unique_time, ia] = unique(original_time);
                current_struct_new.data = clean_data(ia);
                current_struct_new.time = unique_time;

                all_inputs{k} = current_struct_new;
    
        end
    end

    %% Smooth data% Mostly useless now that below anti-aliasing is applied. Can still be used but I don't see many cases in which you want to use it.
    if ~(smooth_data==false)
        for k = 1:num_structures
            %TODO, this doesn't perfectly work with the filtered out data since it is based on window size not time range.
            current_struct_original = all_inputs{k};

            delta_time = delta_time_of_arrays(k);
            %window_size = 2*ceil(smooth_data/delta_time)+1;
            window_size = ceil(smooth_data/delta_time);

            if window_size < 2
                warning("Warning small window size for smoothing.")
            end

            current_struct_original.data = movmean(current_struct_original.data, window_size);

            all_inputs{k} = current_struct_original;
        end
    end

    %% Do interpolation
    all_outputs_fix = cell(1, num_structures);
    
    for k = 1:num_structures
        current_struct_original = all_inputs{k};
        all_outputs_fix{k} = struct('time', [], 'data', []); % Initialize output structure
        
        % Ensure time vectors are column vectors
        original_time_col = current_struct_original.time(:);
        original_data = current_struct_original.data;

        %{
        % Anti-aliasing: If down-sampling, average data within the new sample window
        dt_original = median(diff(original_time_col));
        if new_sample_time > dt_original * 1.05 % Allow 5% margin
            %disp("Doing anti aliasing");disp(k);disp(new_sample_time);disp(dt_original);
            window_size = round(new_sample_time / dt_original);
            if window_size > 1
                original_data = movmean(original_data, window_size, 1);
            end
        end
        %}

        % Uses polyphase anti-aliasing filter
        dt_original = delta_time_of_arrays(k);
        if new_sample_time > dt_original * 1.05
            % Convert rates to integers for resample(data, p, q)
            %disp("Doing anti aliasing");disp(k);disp(new_sample_time);disp(dt_original);
            [p, q] = rat(dt_original / new_sample_time);
            original_data = resample(original_data, p, q);
            % Resample creates its own time vector; use interp1 below to ensure alignment with common_time_vector.
            original_time_col = (0:length(original_data)-1)' * new_sample_time + original_time_col(1);
        end
        
        % Perform interpolation
        try
            % interp1 interpolates each column of `original_data` independently if it's a matrix.
            interpolated_data = interp1(original_time_col, original_data, common_time_vector, 'linear');% Maybe try another interpolation method? https://nl.mathworks.com/help/matlab/ref/double.interp1.html#btwp6lt-1-method
            
            all_outputs_fix{k}.time = common_time_vector;
            all_outputs_fix{k}.data = interpolated_data;
        catch ME
            error(['Could not interpolate data for structure at index ', num2str(k), '. Error: ', ME.message]);
        end
    end

    %% Do optional normalization
    for k = 1:num_structures
        current_struct = all_outputs_fix{k};
        data_to_normalize = current_struct.data;
        
        offset = 0; % The value to be subtracted (min value or 0).
        scaling = 1;  % The value to divide by (the range or the max value).

        % --- Part 1: Determine the offset and scale based on the method ---
        if ismatrix(normalize_data) && ~isscalar(normalize_data)
            offset = normalize_data(k, 1);
            scaling = normalize_data(k, 2);

        elseif normalize_data == 2
            offset = min(data_to_normalize(:));
            max_value = max(data_to_normalize(:));
            offset = offset;
            scaling = max_value - offset;

        elseif normalize_data == 1
            max_value = max(data_to_normalize(:));
            scaling = max_value;
        end
        
        % --- Part 2: Apply the calculated normalization ---
        if scaling == 0
            warning('Normalization scale is zero for structure %d. Data may become NaN or Inf.', k);
        end
        
        current_struct.data = (data_to_normalize - offset) / scaling;
        current_struct.normalisation_offset_used = offset;
        current_struct.normalisation_scaling_used = scaling;
        
        % Update the cell array with the modified structure.
        all_outputs_fix{k} = current_struct;
    end

end