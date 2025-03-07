%% Step 1 Data extraction from tables, adding flags for subject and condition
! filename_excel1 = "subject_condition_data.xlsx. "; %adjust for the correct name
! filename_excel2 = "combined_eye_stim_freq.xlsx"; %adjust for the correct name
data_table1 = readtable(filename_excel1); %=trigger data
data_table2 = readtable(filename_excel2);% = eye data

%% cell error
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
%% Step 2: Parse filenames to extract subject, condition, and frequency, complemented data arrays with this information
% Initialize arrays to store extracted information
filenames_table1 = data_table1{:, 1};
filenames_table2 = data_table2{:, 1};
% Initialize arrays to store extracted information
subjects_table1 = cell(size(filenames_table1));
conditions_table1 = cell(size(filenames_table1));

% Loop through filenames in table 1
for i = 1:length(filenames_table1)
    % Extract subject, condition, and frequency information from filename
    subject_match = regexp(filenames_table1{i}, '(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|016|017|018)', 'match', 'ignorecase');
    if ~isempty(subject_match)
        subjects_table1{i} = lower(char(subject_match{1})); % Store the first match
    else
        subjects_table1{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
    condition_match = regexp(filenames_table1{i}, '(PreWake|High|Low|Mid)', 'match', 'ignorecase');
    if ~isempty(condition_match)
        conditions_table1{i} = lower(char(condition_match{1})); % Store the first match
    else
        conditions_table1{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
end

data_table1.subjects_table1 = subjects_table1;
data_table1.conditions_table1 = conditions_table1;

% Initialize arrays to store extracted information
subjects_table2 = cell(size(filenames_table2));
conditions_table2 = cell(size(filenames_table2));

% Loop through filenames in table 2
for i = 1:length(filenames_table2)
    % Extract subject, condition, and frequency information from filename
    subject_match = regexp(filenames_table2{i}, '(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|016|017|018)', 'match', 'ignorecase');
    if ~isempty(subject_match)
        subjects_table2{i} = lower(char(subject_match{1})); % Store the first match
    else
        subjects_table2{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
    condition_match = regexp(filenames_table2{i}, '(PreWake|High|Low|Mid)', 'match', 'ignorecase');
    if ~isempty(condition_match)
        conditions_table2{i} = lower(char(condition_match{1})); % Store the first match
    else
        conditions_table2{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
end

data_table2.subjects_table2 = subjects_table2;
data_table2.conditions_table2 = conditions_table2;
%%

% Step 3: filtered_data is formed, is data structure with 12 fields, for each condition, for each parameter 
% Extract rows for sum_nr_eye_movements
nr_eye_movements_data = data_table2(strcmp(data_table2.ArrayType, 'sum_nr_eye_movements'), :);
% Extract rows for sum_avg_abs_amplitude_change
avg_abs_amplitude_change_data = data_table2(strcmp(data_table2.ArrayType, 'sum_avg_abs_amplitude_change'), :);
% Extract rows for sum_combined_measure
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
        nr_eye_movements_data(strcmp(nr_eye_movements_data.conditions_table2, condition), :);
    % Filter for avg_abs_amplitude_change_data
    filtered_data.(['avg_abs_amplitude_change_data_' condition]) = ...
        avg_abs_amplitude_change_data(strcmp(avg_abs_amplitude_change_data.conditions_table2, condition), :);

    % Filter for combined_measure_data
    filtered_data.(['combined_measure_data_' condition]) = ...
        combined_measure_data(strcmp(combined_measure_data.conditions_table2, condition), :);
end

% Step 4: Repeated measures ANOVA
% Combine data_table1 and data_table2 to get the necessary columns in one table
merged_data = combined_measure_data;



%%
% Initialize arrays to store data for ANOVA
subjects = {};
conditions = {};
frequencies = [];
eye_movements = [];

% Loop through each row of merged_data
for i = 1:size(merged_data, 1)
    subject = merged_data.subjects_table2{i};
    condition = merged_data.conditions_table2{i};
    flicker_data = data_table1{strcmp(data_table1.subjects_table1, subject) & strcmp(data_table1.conditions_table1, condition), 2:45}; % Frequency data from data_table1
    eye_movement_data = merged_data{i, 3:46}; % Eye movement data from data_table2

    for freq_idx = 1:length(flicker_data)
        frequency = flicker_data(freq_idx);
        eye_movement = eye_movement_data(freq_idx);

        subjects = [subjects; subject];
        conditions = [conditions; condition];
        frequencies = [frequencies; frequency];
        eye_movements = [eye_movements; eye_movement];
    end
end
%%
% Create a table for ANOVA
anova_data = table(subjects, conditions, frequencies, eye_movements, ...
                   'VariableNames', {'Subject', 'Condition', 'Frequency', 'CombinedMeasure'});


% Group by Subject, Condition, and Frequency and calculate the mean of CombinedMeasure
grouped_data = groupsummary(anova_data, {'Subject', 'Condition', 'Frequency'}, 'mean', 'CombinedMeasure');
% Rename 'GroupCount' column to 'NumberOfObservations' and then remove it
grouped_data = removevars(grouped_data, 'GroupCount');

% Specify the filename for the Excel file
filename = 'grouped_data.xlsx';

% Save the table as an Excel file
writetable(grouped_data, filename);


%% Subsequent frequency stimulation categorization done with simple formulae in excell
grouped_data_updated = readtable('grouped_data_updated.xlsx')

%% Across different frequencies
% Convert 'FrequencyCategory' and 'Subject' to categorical variables
grouped_data_updated.Frequency_Category = categorical(grouped_data_updated.Frequency_Category);
grouped_data_updated.Subject = categorical(grouped_data_updated.Subject);

% Get unique subjects and frequency categories
freqCategories = categories(grouped_data_updated.Frequency_Category);
subjects = categories(grouped_data_updated.Subject);

% Convert frequency categories and subjects to cell arrays of character vectors
freqCategoriesCell = cellstr(freqCategories);
subjectsCell = cellstr(subjects);

% Initialize the pivot table
numSubjects = numel(subjectsCell);
numFreqCategories = numel(freqCategoriesCell);

% Create a cell array to hold the data
dataArray = cell(numSubjects, numFreqCategories + 1);
dataArray(:, 1) = subjectsCell; % Set the first column to subjects

% Fill the cell array with NaN values
dataArray(:, 2:end) = num2cell(NaN(numSubjects, numFreqCategories));

% Fill the dataArray with the mean_CombinedMeasure values
for i = 1:numSubjects
    for j = 1:numFreqCategories
        freqCat = freqCategoriesCell{j};
        % Retrieve values for each subject-frequency category pair
        value = grouped_data_updated.mean_CombinedMeasure(grouped_data_updated.Subject == subjects(i) & ...
                                                          grouped_data_updated.Frequency_Category == freqCat);
        if ~isempty(value)
            dataArray{i, j + 1} = value;
        end
    end
end

% Convert dataArray to a table
pivotData2 = cell2table(dataArray);
pivotData2.Properties.VariableNames{'dataArray1'} = 'Subject';
for k = 2:numFreqCategories + 1
    pivotData2.Properties.VariableNames{sprintf('dataArray%d', k)} = freqCategoriesCell{k-1};
end

% Flatten any nested arrays in the table (if necessary)
% Here we assume that each cell contains either a numeric value or NaN
for i = 2:width(pivotData2)  % Start from column 2 to exclude 'Subject'
    pivotData2{:, i} = cellfun(@(x) mean(x, 'omitnan'), pivotData2{:, i}, 'UniformOutput', false);
end

% New desired column order
newOrder = {'baseline', 'low', 'medium', 'high', 'other'};

% Reorder the columns in pivotData2
pivotData2 = pivotData2(:, ['Subject', newOrder]);

% Convert the table to a matrix, excluding the 'Subject' column
dataMatrix = cell2mat(pivotData2{:, 2:end});

% Replace NaN values with zero
dataMatrix(isnan(dataMatrix)) = 0;

% Perform Friedman test
[p_friedman, tbl_friedman, stats_friedman] = friedman(dataMatrix, 1, 'off');

% Display Friedman Test results
disp('Friedman Test results:');
disp(tbl_friedman);
disp(['p-value = ', num2str(p_friedman)]);

% Perform post-hoc test (Dunn's test)
c_dunn = multcompare(stats_friedman, 'Display', 'off');

% Display Dunn's Test results
disp('Dunn''s Test results:');
disp(c_dunn);

% Plotting

% Extract the unique subjects
subjects = pivotData2.Subject;

% Create a matrix with each row corresponding to a subject and each column to a condition
conditionLabels = pivotData2.Properties.VariableNames(2:end);
numConditions = numel(conditionLabels);

% Create a new figure
figure;

% Plot boxplot for each condition
boxplot(dataMatrix, 'Labels', conditionLabels);
title('Friedman Test Results');
xlabel('Frequency Sitmulation Category');
ylabel('Combinatory Measure');
hold on;

% Plot individual data points and lines connecting conditions for each subject
for i = 1:height(pivotData2)
    % Get the data for the current subject
    subjectData = dataMatrix(i, :);
    
    % Define x positions for the conditions
    x_positions = 1:numConditions;
    
    % Plot the individual data points
    scatter(x_positions, subjectData, 'filled');
    
    % Plot lines connecting the conditions for the same subject
    if numConditions > 1
        plot(x_positions, subjectData, '-o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5, 'LineWidth', 0.5);
    end
end

hold off;

%%


% Convert to categorical
% grouped_data_updated.Subject = categorical(grouped_data_updated.Subject);
% grouped_data_updated.Condition = categorical(grouped_data_updated.Condition);
% grouped_data_updated.Frequency = categorical(grouped_data_updated.Frequency_Category);
% 
% % Define unique levels of each within-subject factor
% unique_conditions = unique(grouped_data_updated.Condition);
% unique_frequencies = unique(grouped_data_updated.Frequency);
% 
% % Define the within-subjects design
% within_design = table(grouped_data_updated.Condition, 'VariableNames', {'Condition'});
% 
% % Fit the repeated measures model and perform the ANOVA for NrEyeMovements
% rm = fitrm(anova_data, 'mean_CombinedMeasure~1+Frequency+Condition+Condition*Frequency', 'WithinDesign', within_design);
% ranova_results = ranova(rm, 'WithinModel', 'Condition*Frequency');
% disp(ranova_results);


%MORGEN - van alle catgorien --> per flciekr conidtie , per prop, oer
%beide????


%% This code was used to analyse the differences between the effect of the prewake IAF or the condition IAF on the eye movement. 
data_IAF_prewake =readtable('IAF prewake CM.xlsx');
data_IAF_condition = readtable('IAF Condition CM.xlsx');

data=data_IAF_condition;

% Extract conditions and subjects
conditions = data{:, 2};  % Assuming conditions are in the 2nd column
subjects = data{:, 1};    % Assuming subjects are in the 1st column
values = data{:, 3};      % Assuming the values are in the 3rd column

% Extract data for 'prewake' condition
prewake_data =  data(strcmp(data{:, 2}, 'prewake'), :);
% Extract data for 'low' condition
low_data = data(strcmp(data{:, 2}, 'low'), :);
% Extract data for 'high' condition
high_data = data(strcmp(data{:, 2}, 'high'), :);
% Extract data for 'mid' condition
mid_data = data(strcmp(data{:, 2}, 'mid'), :);

% Get unique subjects
subjects = unique(data.Subject);

% Initialize arrays to store data for each condition
prewake_values = nan(length(subjects), 1);
low_values = nan(length(subjects), 1);
mid_values = nan(length(subjects), 1);
high_values = nan(length(subjects), 1);

% Fill arrays with average values for each subject
for i = 1:length(subjects)
    subject = subjects{i};
    
    % Compute mean values for each condition
    prewake_values(i) = mean(prewake_data.mean_CombinedMeasure(strcmp(prewake_data.Subject, subject)));
    low_values(i) = mean(low_data.mean_CombinedMeasure(strcmp(low_data.Subject, subject)));
    mid_values(i) = mean(mid_data.mean_CombinedMeasure(strcmp(mid_data.Subject, subject)));
    high_values(i) = mean(high_data.mean_CombinedMeasure(strcmp(high_data.Subject, subject)));
end

% Create a matrix with each row corresponding to a subject and each column to a condition
data_matrix = [prewake_values, low_values, mid_values, high_values];

% Replace NaN values with zero
data_matrix(isnan(data_matrix)) = 0;
% Perform Friedman test
[p_friedman, tbl_friedman, stats_friedman] = friedman(data_matrix, 1, 'off');

% Display Friedman Test results
disp('Friedman Test results:');
disp(tbl_friedman);
disp(['p-value = ', num2str(p_friedman)]);
% Perform Friedman test
[p_friedman, tbl_friedman, stats_friedman] = friedman(data_matrix, 1, 'off');

% Perform post-hoc test (Dunn's test)
c_dunn = multcompare(stats_friedman, 'Display', 'off');

% Display results
disp('Dunn''s Test results:');
disp(c_dunn);

%% per subject plot
% Extract unique subjects
subjects = unique(data.Subject);

% Plot boxplot for mean_CombinedMeasure1
figure;
boxplot([prewake_data.mean_CombinedMeasure; low_data.mean_CombinedMeasure; mid_data.mean_CombinedMeasure; high_data.mean_CombinedMeasure], ...
    [ones(size(prewake_data, 1), 1); 2*ones(size(low_data, 1), 1); 3*ones(size(mid_data, 1), 1); 4*ones(size(high_data, 1), 1)]);
title('Flicker Stimulation'); %visual flciker stimulatin or resting state dep on current data
% Combinatory Measure Resting State (per Subject)'
xlabel('Anesthesia Level');
ylabel('Combinatory Measure (NxI)');
set(gca, 'XTickLabel', {'prewake', 'low', 'mid', 'high'});
hold on;

% Plot individual data points and lines connecting the conditions
for i = 1:length(subjects)
    subject = subjects{i};
    
    % Get the data for the current subject
    prewake_value = mean(prewake_data.mean_CombinedMeasure(strcmp(prewake_data.Subject, subject)));
    low_value = mean(low_data.mean_CombinedMeasure(strcmp(low_data.Subject, subject)));
    mid_value = mean(mid_data.mean_CombinedMeasure(strcmp(mid_data.Subject, subject)));
    high_value = mean(high_data.mean_CombinedMeasure(strcmp(high_data.Subject, subject)));
    
    % Combine values into an array and remove empty values
    values = [prewake_value; low_value; mid_value; high_value];
    x_positions = 1:4; % Define x positions for the points
    
    % Plot the individual data points
    scatter(x_positions, values, 'filled');
    
    % Plot lines connecting the conditions for the same subject
    if length(values) > 1
        plot(x_positions, values, '-o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5, 'LineWidth', 0.3);
    end
end















