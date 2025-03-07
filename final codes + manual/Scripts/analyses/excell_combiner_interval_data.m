% Define the root folder
root_folder = '/home/aventurero/DonnaH'; ! correct for your own main path

% Specify the conditions
conditions = {'PreWake', 'PropLow', 'PropHigh', 'PropMid'}; !correct if necessary

% Initialize cell array to store combined data
!combined_data = readtable('/home/aventurero/DonnaH/Sub005/PreWake/eye_mov_per_interval.xlsx'); %correct for your path and file name

% Loop through each subject number
for subj = 5:15
    % Generate the subject folder name
    subject_folder = fullfile(root_folder, sprintf('Sub%03d', subj));

    % Loop through each condition
    for cond = 1:length(conditions)
        condition_folder = fullfile(subject_folder, conditions{cond});
        excel_file = fullfile(condition_folder, 'eye_mov_per_interval.xlsx');

        % Check if the condition folder exists
        if isfolder(condition_folder)
            % Check if the Excel file exists
            if isfile(excel_file)
                % Read the data from the Excel file
                data = readtable(excel_file);

                % Ensure both combined_data and data have the same number of variables
                if ~isempty(combined_data)
                    combined_var_names = combined_data.Properties.VariableNames;
                    data_var_names = data.Properties.VariableNames;

                    % Add missing columns to combined_data
                    for i = 1:length(data_var_names)
                        if ~ismember(data_var_names{i}, combined_var_names)
                            combined_data.(data_var_names{i}) = NaN(height(combined_data), 1);
                        end
                    end

                    % Add missing columns to data
                    for i = 1:length(combined_var_names)
                        if ~ismember(combined_var_names{i}, data_var_names)
                            data.(combined_var_names{i}) = NaN(height(data), 1);
                        end
                    end

                    % Reorder columns in data to match combined_data
                    data = data(:, combined_data.Properties.VariableNames);
                end

                % Append the data to the combined data
                combined_data = [combined_data; data];
            else
                % If the Excel file doesn't exist, display a warning
                warning('File %s does not exist.', excel_file);
            end
        else
            % If the condition folder doesn't exist, display a warning
            warning('Folder %s does not exist.', condition_folder);
        end
    end
end

%%
% Combine all tables into one
if ~isempty(combined_data)
    % Write the combined data to an Excel file
    output_filename = 'combined_eye_stim_freq.xlsx';
    writetable(combined_data, output_filename);

    disp(['Combined Excel file created: ' output_filename]);
else
    disp('No data found to combine.');
end
