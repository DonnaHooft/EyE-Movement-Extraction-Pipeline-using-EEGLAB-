root_folder = '/home/aventurero/Documents/klifoanae/propofolstudy/Experiments' ;

% Specify the conditions
conditions = {'PreWake_InsideScannerfMRI', 'PropofolHigh', 'PropofolLow', 'PropofolMid'};

% Create a cell array to store the data
data = cell(0);

% Loop through each subject folder
for sub = 1:16
    % Generate the subject folder name
    ! subject_folder = sprintf('sub-%03d', sub); %adjust for number of subjects
   
    % Loop through each condition
    for cond = 1:length(conditions)
        % Generate the condition folder path
        condition_folder = fullfile(root_folder, subject_folder, 'CalculateIAF', conditions{cond});
       
        % Read the ConditionOrder.txt file
        txt_file = fullfile(condition_folder, 'ConditionOrder.txt');
        if exist(txt_file, 'file')
            % Read the data from the text file
            condition_data = importdata(txt_file, ',');
           
            % Combine subject and condition name
            combined_name = sprintf('%s_%s', subject_folder, conditions{cond});
           
            % Add the combined name and data to the cell array
            data = [data; {combined_name}, condition_data];
        else
            % If the file doesn't exist, display a warning
            warning('File %s does not exist.', txt_file);
        end
    end
end

% Create a cell array for the header
header = cell(1, size(data, 2));
header{1} = 'Subject_Condition';
for i = 2:length(header)
    header{i} = sprintf('Column%d', i-1);
end

% Write the data to an Excel file
filename_excel = 'subject_condition_data.xlsx';
writecell([header; data], filename_excel);
disp(['Excel file created: ' filename_excel]);