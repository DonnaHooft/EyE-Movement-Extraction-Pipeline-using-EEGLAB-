%% Step 1 Data extraction from tables, adding flags for subject and condition
! filename_excel1 = "subject_condition_data.xlsx. "; %adjust for the correct name
! filename_excel2 = "combined_eye_stim_freq.xlsx"; %adjust for the correct name
data_table1 = readtable(filename_excel1); %=trigger data
data_table2 = readtable(filename_excel2);% = eye data

%% Cell error
% Something odd was going one with the coolumn 100-120 seconds, for some
% reason the numbers were sored as a cell, causing issue slater on in the
% pipeline, this section restores that error
% Assuming the column name is 'x100_120Seconds'
column_name = 'x100_120Seconds';
% Step 1: Extract the column data
column_data = data_table2.(column_name);
% Step 2: Check if the column data is a cell array
if iscell(column_data)
    % Convert cell array of strings to numeric doubles
    numeric_data = cellfun(@str2double, column_data);

    % Step 3: Update the table with the numeric data
    data_table2.(column_name) = numeric_data;
else
    error('The column is not a cell array of strings.');
end

column_name = 'x440_460Seconds';
% Step 1: Extract the column data
column_data = data_table2.(column_name);

% Step 2: Check if the column data is a cell array
if iscell(column_data)
    % Convert cell array of strings to numeric doubles
    numeric_data = cellfun(@str2double, column_data);

    % Step 3: Update the table with the numeric data
    data_table2.(column_name) = numeric_data;
else
    error('The column is not a cell array of strings.');
end

%% % Step 2: Parse filenames to extract subject, condition, and frequency, complemented data arrays with this information
% information
filenames_table1 = data_table1(:, 1);
filenames_table2 = data_table2(:, 1);
% Initialize arrays to store extracted information
subjects_table1 = cell(size(filenames_table1));
conditions_table1 = cell(size(filenames_table1));

% Loop through filenames in table 1
for i = 1:size(filenames_table1, 1)
    % Extract subject, condition, and frequency information from filename
    % -> ADJUST FOR YOUR SUBJECTS!
    subject_match= regexp(filenames_table1{i,1}, '(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|016|017|018)', 'match', 'ignorecase');
    if ~isempty(subject_match)
        subjects_table1{i} = lower(char(subject_match{1})); % Store the first match
    else
        subjects_table1{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
    condition_match = regexp(filenames_table1{i,1}, '(PreWake|High|Low|Mid)', 'match', 'ignorecase');
    if ~isempty(condition_match)
        conditions_table1{i} = lower(char(condition_match{1})); % Store the first match
    else
        conditions_table1{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
end

data_table1=[data_table1, subjects_table1, conditions_table1];
subjects_table2 = cell(size(filenames_table2));
conditions_table2 = cell(size(filenames_table2));
% Loop through filenames in table 1
for i = 1:size(filenames_table2, 1)
    % Extract subject, condition, and frequency information from filename
    subject_match= regexp(filenames_table2{i,1}, '(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|016|017|018)', 'match', 'ignorecase');
    if ~isempty(subject_match)
        subjects_table2{i} = lower(char(subject_match{1})); % Store the first match
    else
        subjects_table2{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
    condition_match = regexp(filenames_table2{i,1}, '(PreWake|High|Low|Mid)', 'match', 'ignorecase');
    if ~isempty(condition_match)
        conditions_table2{i} = lower(char(condition_match{1})); % Store the first match
    else
        conditions_table2{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
end
data_table2=[data_table2, subjects_table2, conditions_table2];
%% Step 3:  filtered_data is formed, is data structure whith 12 fields, for each conidtions, for each paramater 
% Extract rows for sum_nr_eye_movements
nr_eye_movements_data = data_table2(strcmp(data_table2.ArrayType, 'sum_nr_eye_movements'), :);
% Extract rows for sum_avg_abs_amplitude_change
avg_abs_amplitude_change_data = data_table2(strcmp(data_table2.ArrayType, 'sum_avg_abs_amplitude_change'), :);
% Extract rows for sum_combined measure
combined_measure_data = data_table2(strcmp(data_table2.ArrayType, 'sum_combined_measure_interval'), :);
% Define the conditions
conditions = {'prewake', 'low', 'mid', 'high'};

% Initialize a structure to hold the filtered data
filtered_data = struct();

% Loop through each condition and filter the data
for i = 1:length(conditions)
    condition = conditions{i};
    
    % Filter for nr_eye_movements_data
    filtered_data.(['nr_eye_movements_data_' condition]) = ...
        nr_eye_movements_data(strcmp(nr_eye_movements_data.Var48, condition), :);
     % Filter for avg_abs_amplitude_change_data
    filtered_data.(['avg_abs_amplitude_change_data_' condition]) = ...
        avg_abs_amplitude_change_data(strcmp(avg_abs_amplitude_change_data.Var48, condition), :);
    
    % Filter for combined_measure_data
    filtered_data.(['combined_measure_data_' condition]) = ...
        combined_measure_data(strcmp(combined_measure_data.Var48, condition), :);
end

%% repeated measures anova
% Combine data_table1 and data_table2 to get the necessary columns in one table
merged_data = vertcat(combined_measure_data);
merged_data2 = vertcat(nr_eye_movements_data);
merged_data3 = vertcat(avg_abs_amplitude_change_data);

% Initialize arrays to store data for ANOVA
subjects = {};
conditions = {};
frequencies = [];
eye_movements = [];
avg_abs_amplitude_changes = [];
combined_measures = [];

% Loop through each row of merged_data
for i = 1:size(merged_data, 1)
    subject = merged_data.Var47{i};
    condition = merged_data.Var48{i};
    flicker_data = data_table1(:, 2:end);
    eye_movement_data = merged_data2{i, 3:46};
    amplitude_data=  merged_data3{i, 3:46}
    combined_data = merged_data{i, 3:46}
    % Eye movement data from data_table2

    for freq_idx = 1:(width(flicker_data)-3)
        frequency = flicker_data{i,freq_idx};
        eye_movement = eye_movement_data(freq_idx);
        ampls= amplitude_data(freq_idx);
        cm = combined_data(freq_idx);
        
        subjects = [subjects; subject];
        conditions = [conditions; condition];
        frequencies = [frequencies; frequency];
        eye_movements = [eye_movements; eye_movement];
        avg_abs_amplitude_changes = [avg_abs_amplitude_changes;ampls];
        combined_measures = [ combined_measures;cm];
    end
end

% Check the structure of the variables before creating the table
disp(struct2table(struct('subjects', {subjects}, 'conditions', {conditions}, 'frequencies', frequencies, 'eye_movements', eye_movements)));

% Create a table for ANOVA
anova_data = table(subjects, conditions, frequencies, eye_movements, avg_abs_amplitude_changes, combined_measures, ...
                   'VariableNames', {'Subject', 'Condition', 'Frequency', 'EyeMovements', 'AvgAbsAmplitudeChanges', 'CombinedMeasures'});




%% All conditions in one figure
% Define the conditions and corresponding colors
conditions = {'prewake', 'low', 'mid', 'high'};
colors = {'b', 'g', 'y', 'r'};
markers = {'o', 's', 'd', '^'}; % Different markers for different conditions

% Define array types and their corresponding variable names
array_types = {'nr_eye_movements_data', 'avg_abs_amplitude_change_data', 'combined_measure_data'};
array_labels = {'Number of Eye Movements', 'Average Amplitude Change', 'Combined Measure'};

% Loop through each array type and condition to generate the plots
for i = 1:length(array_types)
    array_type = array_types{i};
    array_label = array_labels{i};
    
    % Create a new figure for the current array type
    figure;
    hold on;

    % Store plot handles for the legend
    scatter_handles = [];
    trendline_handles = [];

    for j = 1:length(conditions)
        condition = conditions{j};
        data_key = [array_type '_' condition];
        
        % Get the filtered data for the current array type and condition
        current_data = filtered_data.(data_key);
        
        % Initialize arrays to store frequencies and eye movement data
        unique_frequencies = [];
        eye_movement_data = [];
        
        % Loop through each row in the current data
        for k = 1:size(current_data, 1)
            subject = current_data.Var47{k};
            condition = current_data.Var48{k};
            
            % Find the corresponding row in data_table1
            idx = strcmp(data_table1.Var46, subject) & strcmp(data_table1.Var47, condition);
            
            if any(idx)
                flicker_data = data_table1{idx, 2:45};
                for freq_idx = 1:size(flicker_data, 2)
                    frequency = flicker_data(freq_idx);
                    
                    % Append the frequency to the list of unique frequencies
                    unique_frequencies = [unique_frequencies; frequency];
                    
                    % Sum the corresponding eye movements
                    eye_movement_count = current_data{k, 3:46}(freq_idx);
                    eye_movement_data = [eye_movement_data; eye_movement_count];
                end
            end
        end
        
        % Count the occurrences of each unique frequency
        [unique_frequencies, ~, freq_indices] = unique(unique_frequencies);
        frequency_counts = accumarray(freq_indices, 1);
        
        % Sum eye movement data for each frequency
        eye_movement_data_sum = accumarray(freq_indices, eye_movement_data, [], @sum);
        
        % Divide by the frequency counts to get the average
        eye_movement_data_avg = eye_movement_data_sum ./ frequency_counts;
        
        % Remove baseline (first value) from eye movement data
        eye_movement_data_avg = eye_movement_data_avg - eye_movement_data_avg(1);
        
        % Remove the baseline (first value) from frequencies
        unique_frequencies = unique_frequencies(2:end);
        eye_movement_data_avg = eye_movement_data_avg(2:end);
            
        % Create a subplot for the current condition
        scatter_handle = scatter(unique_frequencies, eye_movement_data_avg, colors{j}, 'filled');
        scatter_handles = [scatter_handles; scatter_handle];

        % Fit a polynomial trendline
        p = polyfit(unique_frequencies, eye_movement_data_avg, 2);
        fitted_values = polyval(p, unique_frequencies);
        
        % Plot the trendline
        trendline_handle = plot(unique_frequencies, fitted_values, colors{j}, 'LineWidth', 1);
        trendline_handles = [trendline_handles; trendline_handle];
        
        xlabel('Visual Stimulation Frequency');
        ylabel(array_label);
        grid on;
    end
    
    hold off;
    % Manually set the legend entries
    legend([scatter_handles; trendline_handles], ...
        {'Data prewake', 'Data low', 'Data mid', 'Data high', ...
         'Trendline prewake', 'Trendline low', 'Trendline mid', 'Trendline high'});
    
    % Add a main title for the entire figure
    sgtitle('');
end
%% SUBPLOTS clear all + rerun (else data is removed)
nr_eye_movements_data = data_table2(strcmp(data_table2.ArrayType, 'sum_nr_eye_movements'), :);
% Extract rows for sum_avg_abs_amplitude_change
avg_abs_amplitude_change_data = data_table2(strcmp(data_table2.ArrayType, 'sum_avg_abs_amplitude_change'), :);
% Extract rows for sum_combined measure
combined_measure_data = data_table2(strcmp(data_table2.ArrayType, 'sum_combined_measure_interval'), :);
% Define the conditions
conditions = {'prewake', 'low', 'mid', 'high'};
% Initialize a structure to hold the filtered data
filtered_data = struct();

% Loop through each condition and filter the data
for i = 1:length(conditions)
    condition = conditions{i};
    
    % Filter for nr_eye_movements_data
    filtered_data.(['nr_eye_movements_data_' condition]) = ...
        nr_eye_movements_data(strcmp(nr_eye_movements_data.Var48, condition), :);
    % Filter for avg_abs_amplitude_change_data
    filtered_data.(['avg_abs_amplitude_change_data_' condition]) = ...
        avg_abs_amplitude_change_data(strcmp(avg_abs_amplitude_change_data.Var48, condition), :);
    % Filter for combined_measure_data
    filtered_data.(['combined_measure_data_' condition]) = ...
        combined_measure_data(strcmp(combined_measure_data.Var48, condition), :);
end
% Define array types and their corresponding variable names
array_types = {'nr_eye_movements_data', 'avg_abs_amplitude_change_data', 'combined_measure_data'};
array_labels = {'Number of Eye Movements', 'Average Amplitude Change', 'Combined Measure'};

% Loop through each array type and condition to generate the plots
% Loop through each array type and condition to generate the plots
for i = 1:length(array_types)
    array_type = array_types{i};
    array_label = array_labels{i};
    
    % Create a new figure for the current array type
    figure;
    
    for j = 1:length(conditions)
        condition = conditions{j};
        data_key = [array_type '_' condition];
        
        % Get the filtered data for the current array type and condition
        current_data = filtered_data.(data_key);
        
        % Initialize arrays to store frequencies and eye movement data
        unique_frequencies = [];
        eye_movement_data = [];
        
        % Loop through each row in the current data
        for k = 1:size(current_data, 1)
            subject = current_data.Var47{k};
            condition = current_data.Var48{k};
            
            % Find the corresponding row in data_table1
            idx = strcmp(data_table1.Var46, subject) & strcmp(data_table1.Var47, condition);
            
            if any(idx)
                flicker_data = data_table1{idx, 2:45};
                for freq_idx = 1:size(flicker_data, 2)
                    frequency = flicker_data(freq_idx);
           % Append the frequency to the list of unique frequencies
                    unique_frequencies = [unique_frequencies; frequency];
             % Sum the corresponding eye movements
                    eye_movement_count = current_data{k, 3:46}(freq_idx);
                    eye_movement_data = [eye_movement_data; eye_movement_count];
                end
            end
        end
        
        % Count the occurrences of each unique frequency
        [unique_frequencies, ~, freq_indices] = unique(unique_frequencies);
        frequency_counts = accumarray(freq_indices, 1);
        
        % Sum eye movement data for each frequency
        eye_movement_data_sum = accumarray(freq_indices, eye_movement_data, [], @sum);
        % Divide by the frequency counts to get the average
        eye_movement_data_avg = eye_movement_data_sum ./ frequency_counts;
        % Remove baseline (first value) from eye movement data
        eye_movement_data_avg = eye_movement_data_avg - eye_movement_data_avg(1);
        % Remove the baseline (first value) from frequencies
        unique_frequencies = unique_frequencies(2:end);
        eye_movement_data_avg = eye_movement_data_avg(2:end);
        
% % %         % Remove the second value (outlier) from each array
% %         if strcmp(condition, 'mid')
% %             unique_frequencies(1) = [];
% %             eye_movement_data_avg(1) = [];
% %         end
% %         
% %         if strcmp(condition, 'prewake')
% %             unique_frequencies(end) = [];
% %             eye_movement_data_avg(end) = [];
% %         end

        % Create a subplot for the current condition
        subplot(2, 2, j);
        scatter(unique_frequencies, eye_movement_data_avg, 'filled');
        hold on;
        % Fit a polynomial trendline
        p = polyfit(unique_frequencies, eye_movement_data_avg, 2);
        fitted_values = polyval(p, unique_frequencies);
        % Find the maximum point of the parabola
        max_value = max(fitted_values);
        max_index = find(fitted_values == max_value);
        max_frequency = unique_frequencies(max_index);
        
        % Plot the trendline
        plot(unique_frequencies, fitted_values, 'r-', 'LineWidth', 2);
        hold on;
        
        % Mark the maximum point with a star
        plot(max_frequency, max_value, 'r*', 'MarkerSize', 10);
        
        xlabel('Unique Frequencies');
        ylabel(array_label);
        title([array_label ' vs Flicker Frequencies (' condition ')']);
        legend(array_label, 'Trendline');
        grid on;
        hold off;
    end
    
    % Add a main title for the entire figure
    sgtitle([array_label ' Across Different Conditions']);
end




