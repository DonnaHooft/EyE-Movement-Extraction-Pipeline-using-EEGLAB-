%% for resting state analysis with time weighing factor 
filename_excel = "resting_state_data.xlsx";
data = readtable(filename_excel);
data = data(1:end-2, :);

% Extract conditions from filenames and duration
filenames = data(:, 1);
conditions = cell(size(filenames));
durations = zeros(size(filenames));
subjects_table = cell(size(filenames));
for i = 1:size(filenames, 1)
    % Extract condition from filename using regular expression
    match = regexp(filenames{i,1}, '(prewake|low|mid|high)', 'match', 'ignorecase');
    if ~isempty(match)
        conditions{i} = lower(char(match{1})); % Store the first match
    else
        conditions{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
    % Extract duration from filename using regular expression
    duration_match = regexp(filenames{i,1}, '(\d+)[mM]in', 'tokens', 'once');
    if ~isempty(duration_match)
        durations(i) = str2double(duration_match{1});
    else
        durations(i) = 1; % Default to 1 if no duration found
    end
    
    % Extract subject, condition, and frequency information from filename
    subject_match= regexp(filenames{i,1}, '(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015|0160|017|018)', 'match', 'ignorecase');
    if ~isempty(subject_match)
        subjects_table{i} = lower(char(subject_match{1})); % Store the first match
    else
        subjects_table{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
end

% Concatenate conditions_table and durations_table with your existing data table
conditions_table = table(conditions, 'VariableNames', {'Conditions'});
durations_table = table(durations, 'VariableNames', {'Duration'});
subjects_table2 = table(subjects_table, 'VariableNames', {'Subject'});
data = [data, conditions_table, durations_table, subjects_table2];

expanded_data = [];
for i = 1:size(data, 1)
    duration = data.Duration(i);
    for j = 1:duration
        expanded_data = [expanded_data; data(i, :)];
    end
end

% Now expanded_data contains the weighted data
prewake_data = expanded_data(strcmp(expanded_data.Conditions, 'prewake'), :);
low_data = expanded_data(strcmp(expanded_data.Conditions, 'low'), :);
high_data = expanded_data(strcmp(expanded_data.Conditions, 'high'), :);
mid_data = expanded_data(strcmp(expanded_data.Conditions, 'mid'), :);

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
    prewake_values(i) = mean(prewake_data.NrOfEyeMovements(strcmp(prewake_data.Subject, subject)));
    low_values(i) = mean(low_data.NrOfEyeMovements(strcmp(low_data.Subject, subject)));
    mid_values(i) = mean(mid_data.NrOfEyeMovements(strcmp(mid_data.Subject, subject)));
    high_values(i) = mean(high_data.NrOfEyeMovements(strcmp(high_data.Subject, subject)));
end
%= NrOfEyeMovements1
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

% Extract unique subjects
subjects = unique(data.Subject);

% Plot boxplot for NrOfEyeMovements1
figure;
boxplot([prewake_data.NrOfEyeMovements; low_data.NrOfEyeMovements; mid_data.NrOfEyeMovements; high_data.NrOfEyeMovements], ...
    [ones(size(prewake_data, 1), 1); 2*ones(size(low_data, 1), 1); 3*ones(size(mid_data, 1), 1); 4*ones(size(high_data, 1), 1)]);
title('Resting State - Combinatory Measure (per Subject)');
xlabel('Anesthesia Level');
ylabel('Combinatory Measure (NxI)');
set(gca, 'XTickLabel', {'prewake', 'low', 'mid', 'high'});
hold on;

% Plot individual data points and lines connecting the conditions
for i = 1:length(subjects)
    subject = subjects{i};
    
    % Get the data for the current subject
    prewake_values = prewake_data.NrOfEyeMovements(strcmp(prewake_data.Subject, subject));
    low_values = low_data.NrOfEyeMovements(strcmp(low_data.Subject, subject));
    mid_values = mid_data.NrOfEyeMovements(strcmp(mid_data.Subject, subject));
    high_values = high_data.NrOfEyeMovements(strcmp(high_data.Subject, subject));
    
    % Combine values into an array
    values = [mean(prewake_values), mean(low_values), mean(mid_values), mean(high_values)];
    x_positions = 1:4; % Define x positions for the points
    
    % Plot the individual data points
    scatter(x_positions, values, 'filled', 'DisplayName', subject);
    
    % Plot lines connecting the conditions for the same subject
    if length(values) > 1
        plot(x_positions, values, '-o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5, 'LineWidth', 0.3);
    end
end

hold off;

%% for 15 min flicker general analysis per conditions 
filename_excel ="flicker_data.xlsx";
data = readtable(filename_excel);
data = data(1:end-2, :);
% Extract conditions from filenames
filenames = data(:, 1);
conditions = cell(size(filenames));
subjects_table = cell(size(filenames));
for i = 1:size(filenames)
    % Extract condition from filename using regular expression
    match = regexp(filenames{i,1}, '(prewake|low|mid|high)', 'match', 'ignorecase');
    if ~isempty(match)
        conditions{i} = lower(char(match{1})); % Store the first match
    else
        conditions{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
    
    % Extract subject, condition, and frequency information from filename
    subject_match= regexp(filenames{i,1}, '(001|002|003|004|005|006|007|008|009|010|011|012|013|014|015)', 'match', 'ignorecase');
    if ~isempty(subject_match)
        subjects_table{i} = lower(char(subject_match{1})); % Store the first match
    else
        subjects_table{i} = 'Unknown'; % If no match found, mark as 'Unknown'
    end
end
% Concatenate conditions_table with your existing data table
conditions_table = table(conditions, 'VariableNames', {'Conditions'});
subjects_table2 = table(subjects_table, 'VariableNames', {'Subject'});
data = [data, conditions_table, subjects_table2];

% Extract data for 'prewake' condition
prewake_data =  data(strcmp(data.Conditions, 'prewake'), :);
% Extract data for 'low' condition
low_data = data(strcmp(data.Conditions, 'low'), :);
% Extract data for 'high' condition
high_data = data(strcmp(data.Conditions, 'high'), :);
% Extract data for 'mid' condition
mid_data = data(strcmp(data.Conditions, 'mid'), :);


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
    prewake_values(i) = mean(prewake_data.NrOfEyeMovements(strcmp(prewake_data.Subject, subject)));
    low_values(i) = mean(low_data.NrOfEyeMovements(strcmp(low_data.Subject, subject)));
    mid_values(i) = mean(mid_data.NrOfEyeMovements(strcmp(mid_data.Subject, subject)));
    high_values(i) = mean(high_data.NrOfEyeMovements(strcmp(high_data.Subject, subject)));
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

% Extract unique subjects
subjects = unique(data.Subject);

% Plot boxplot for NrOfEyeMovements1
figure;
boxplot([prewake_data.NrOfEyeMovements; low_data.NrOfEyeMovements; mid_data.NrOfEyeMovements; high_data.NrOfEyeMovements], ...
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
    prewake_value = mean(prewake_data.NrOfEyeMovements(strcmp(prewake_data.Subject, subject)));
    low_value = mean(low_data.NrOfEyeMovements(strcmp(low_data.Subject, subject)));
    mid_value = mean(mid_data.NrOfEyeMovements(strcmp(mid_data.Subject, subject)));
    high_value = mean(high_data.NrOfEyeMovements(strcmp(high_data.Subject, subject)));
    
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

