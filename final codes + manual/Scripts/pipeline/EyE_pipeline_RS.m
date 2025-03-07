clc;
clear all;
%% Allocate main path
%This will be the path wherein the script, subject folders and the 'Matlab_functions' or 'Matlab_plugins' (wherein 'CWRegrTool' and 'eeglab2023.0') folders are stored. 
% As well as your 'average_eye_kernel.mat' and your% 'subject_information.mat'
! mainpath = '/home/zwergwidder/DonnaH' ; mainpath2 = '/home/zwergwidder/DonnaH';

%% List all subject folders
AllSubFolders = dir(fullfile(mainpath, 'Sub*'));
%adjust if the prefix of your subjects folder is not 'Sub'

%% Select specific subjects and sessions or analyze all (if empty)
! selectedSubjects = {'Sub001', 'Sub002','Sub013','Sub016'}; % Add or remove subjects as needed (f.e. 'Sub001' folder name(s))
! selectedSessions = {'PreWake', 'PropofolHigh','PropofolLow','PropofolMid'}; % Add or remove sessions as needed (f.e. 'PreWake' folder name(s)); Note that naming could also be 'PropLow', 'PropHigh', etc. 

if ~isempty(selectedSubjects)
    AllSubFolders = AllSubFolders(ismember({AllSubFolders.name}, selectedSubjects));
end

%% Read the Excel file with channels to be removed 
% adjust if you name or path is different
! load('subject_information.mat');
subject_info_data = data;
%% Create header for Excel file
header = {'Filename', 'Nr of Eye Movements', 'Nr Eye movements per Minute', 'Average Intensity', ...
    'Combinatory Measure 1', 'Combinatory Measure 2', 'Maximum Eye Score'};

% Initialize cell array to store data
data1 = cell(length(AllSubFolders)*8, length(header));
%Load average eye kernel in your workspace to use alter on in convolving
! load('average_eye_kernel.mat');

%ensure to adjust this to the corresponding paths, sometimes some debugging
%is required to ultimately get this to work, in my case adding the path to 
%the plugins here in the pipeline worked out. 
! addpath(genpath('/home/zwergwidder/DonnaH/Matlab_functions/'));
! addpath('/home/zwergwidder/DonnaH/software/CWRegrTool/CWRegrTool/')
! addpath('/home/zwergwidder/DonnaH/software/eeglab2023.0/');

%% Loop through each subject folder
for sub = 1:length(AllSubFolders)
    currentSubjectFolder = fullfile(mainpath, AllSubFolders(sub).name);
    
    % List all sessions within the current subject folder
    sessions = dir(fullfile(currentSubjectFolder, '*'));
    sessions = sessions([sessions.isdir] & ~startsWith({sessions.name}, '.') & ~strcmp({sessions.name}, '.DS_Store'));
    
    if ~isempty(selectedSessions)
        sessions = sessions(ismember({sessions.name}, selectedSessions));
    end
    
    % Loop through each session
    for session = 1:length(sessions)
        currentSessionFolder = fullfile(currentSubjectFolder, sessions(session).name);
        
        % Check if the session folder is not empty
        if numel(dir(fullfile(currentSessionFolder, '*.mat'))) > 0
            % Form the current session folder path
            SessionFolder = currentSessionFolder;
            
            %% start eeglab
            eeglab;
            
            %% Import Data
            cd(currentSessionFolder);
            files = dir(fullfile(currentSessionFolder, '*.mat'));
            % Loop through each file in session folder
            for i = 1:length(files)
                fileName = files(i).name;
                
                %                 EEG = pop_loadbv(SessionFolder, fileName);
                fullFilePath = fullfile(SessionFolder, fileName);
                load(fullFilePath);
                
                % Extract information
                recording_time = EEG.xmax - EEG.xmin;
                num_volume_triggers = sum(strcmp({EEG.event.type}, 'R148'));
                sampling_rate_before = EEG.srate;
                ecg_channel_number = find(strcmp({EEG.chanlocs.labels}, 'ECG'));
                
                %% This script assumes the preprocessign steps as described in the manual, in case you get raw .eeg data, pleas uncomment the following lines
                
                %                 % Remove gradient artifact
                %                 ECG = FindChannelNumberFromLabel(EEG.chanlocs, 'ECG');
                %                 CW1 = FindChannelNumberFromLabel(EEG.chanlocs, 'CW1');
                %                 CW2 = FindChannelNumberFromLabel(EEG.chanlocs, 'CW2');
                %                 CW3 = FindChannelNumberFromLabel(EEG.chanlocs, 'CW3');
                %                 CW4 = FindChannelNumberFromLabel(EEG.chanlocs, 'CW4');
                %
                %                 EEG = pop_fmrib_fastr(EEG,70,5,10,'R148',0,0,0,[],[],0.03,[ECG CW1 CW2 CW3 CW4],'auto');
                %                 % Remove cardioballistic artifact
                %                 assert(ECG == 19);
                %                 EEG = pop_cwregression(EEG,EEG.srate,4,0.021,1,'hann',[CW1 CW2 CW3 CW4],[1:ECG-1 ECG+1:64],'taperedhann',0);
                %
                
                % Remove any data at beginning of recording that occurred before first fMRI volume trigger
                trigger_onsets = FindTriggerOnsets(EEG, 'R148');
                TR = mode(diff(trigger_onsets))/EEG.srate; %fMRI TR (repetition time) was 1 second for single echo EPI; this changed to 1.5 when we changed to multi echo EPI sequence
                EEG = pop_select(EEG, 'point', [trigger_onsets(3) trigger_onsets(end-2)]);
                
                % Remove ECG and CWL channels - RN before filtering
                EEG = pop_select(EEG,'rmchannel',{'ECG','CW1','CW2','CW3','CW4'});

                % Remove the specified substring from fileName
                % in the files received after preprocessing, the filenames
                % had this addition, please comment out the following few
                % rows if this is not the case, or adjsut to the particular
                % addition in your case. 
 !               fileName_cleaned = erase(fileName, '_GS_classic+SVD_CWregression.mat');
                % Replace the file extension with .vhdr
 !               fileName_cleaned = replace(fileName_cleaned, '.mat', '.vhdr');
                % Find the row corresponding to the cleaned fileName
 !               file_row = strcmp(subject_info_data(:, 1), fileName_cleaned);
                
                % Assuming the file names are stored in the first column
                % this sections removes bad channels, if they are store in
                % the 'subject_information.mat'. This is to optimize the
                % ICA, if this is causing issues, comment it out. 
                channels_to_remove = subject_info_data{file_row, 7};   % Assuming the channels to remove are in the 8th column
                % Remove single quotes and leading/trailing spaces
                % Check if channels_to_remove is not empty
                if ~isempty(channels_to_remove)
                    % Split the string into individual parts
                    parts = strsplit(channels_to_remove, ', ');
                    % Join the parts back together with commas
                    channels_to_remove_formatted = strjoin(parts, ', ');
                    % Perform the operation only if channels_to_remove is not empty
                    EEG = pop_select(EEG, 'rmchannel', {channels_to_remove_formatted});
                end

                % Change sampling rate
                EEG = pop_resample(EEG,1000);
                EEGc = EEG;
                % Apply high pass filter to EEG data to remove noise artifacts below 0.5 Hz and above 100 Hz
                EEGc = pop_eegfiltnew(EEGc, 0.5, []);

                %                 % If desired could uncomment to manually inspect the interim file. 
                %                 interimFilename = extractBefore(fileName, '.vhdr');
                %                 save(fullfile(SessionFolder, ['interim_', interimFilename, '.mat']), 'EEGc');
                %                 pop_saveset(EEGc, ['interim_', interimFilename, '.mat'], SessionFolder);
                
                 %% Running ICA, extracting eye movement and metrics
                EEGc = pop_runica(EEGc,'extended',1);
                EEGc = iclabel(EEGc);
                ic_scores = EEGc.etc.ic_classification.ICLabel;
                % Extract components which include eye movement and plotting  | Removing channel noise and bad components
                % Multiply the ICA weights with the ICA sphere to get the unmixing matrix
                unmixing_matrix = EEGc.icaweights * EEGc.icasphere;
                % source activation (Components) = unmixing * channel data
                EEGc.icaact = unmixing_matrix * EEGc.data;
                
                % Extract components which include eye movement and plotting  | Removing channel noise and bad components
                eye_class_scores = EEGc.etc.ic_classification.ICLabel.classifications(:, 3);  % Extract scores for the eye class (class 3)
                max_eye_score = max(eye_class_scores);
                % Find the maximum score for the eye class
                threshold = 0.8*max_eye_score; %0.8 is variable, but sensible from looking at multiple ICAs
                eye_components = (eye_class_scores >= threshold); %from this the final eye
                eye_component_activations = zeros(1, size(EEGc.data, 2));
                
                total_eye_movements = 0;
                avg_amplitude_change_mean = 0;
                avg_eye_movements_per_min = 0;
                total_CM_1 = 0;
                total_CM_2 = 0;
                
 !               if max_eye_score > 0.20  % adjust accordingly to the quality of your data, 
                % if it is better, please increase, if most of your ICA is compormised, 
                % you could opt to lower this, however lower than 0.10 is not recommended. 
                    eye_component_activations = EEGc.icaact(eye_components,:);
                    % Extract the traces of the components labeled as eye movements
                    % eye_component_activations = EEGc.icaact(eye_components,:);
                    [num_eye_components, num_samples] = size(eye_component_activations);
                    timepoints = EEGc.times;
                    % Initialize a matrix to store the rescaled component traces
                    rescaled_eye_component=zeros(num_eye_components,length(timepoints));
                    
                    timepoints = EEGc.times;  % Assuming EEGc.times holds the timepoints
                    sampling_frequency = 1000; % Assuming a sampling frequency of 1000 Hz, adjust as needed
                    % Design & apply bandpass filteraverage eye movement component
                    low_cutoff = 0.5; % Hz
                    high_cutoff = 50; % Hz
                    filter_order = 4; % Filter order (adjust as needed)
                    [b, a] = butter(filter_order, [low_cutoff, high_cutoff] / (sampling_frequency / 2), 'bandpass');
                    
% Iterate over each eye component
% applies average eye movement kernel to data  subsequently rescales to -1 to 1 values
                    
                    for j = 1:num_eye_components
                        % Select individual eye component
                        eye_component = eye_component_activations(j, :);
                        
                        % Apply filtering and smoothing to individual eye component
                        eye_component = double(eye_component);
                        % Define sampling frequency
                        nan_indices = isnan(eye_component);
                        eye_component(nan_indices) = 0;
                        
                        % Apply filtering
                        filtered_eye_component = filtfilt(b, a, eye_component);
                        kernel_filtered_trace= conv(filtered_eye_component, average_eye_kernel, 'same');
                        
                        % Rescale the eye component between -1 and 1
                        max_val = max(abs( kernel_filtered_trace));
                        rescaled_eye_component(j,:) =  kernel_filtered_trace / max_val;
                    end

                    time_interval=10;
                    x_1000 = time_interval * 1000;
                    
                    % Initialize arrays to store final metrics per component
                    num_components = size(rescaled_eye_component, 1);
                    total_eye_movements_per_component = zeros(num_components, 1);
                    avg_amplitude_change_mean_per_component = zeros(num_components, 1);
                    combined_measure_per_component = zeros(num_components, floor(length(rescaled_eye_component(1,:)) / x_1000));
                    
                    for c = 1:num_components
                        % Extract the rescaled eye component data for the current component
                        component = rescaled_eye_component(c, :);
                        min_value = min(component);
                        max_value = max(component);
                        % Rescale smoothed_eye_mean between -1 and 1, this to in order to work wth
                        % different values of eye traces
                        rescaled_eye_mean = -1 + 2 * (component - min_value) / (max_value - min_value);
                        %Ensure 0-centered data
                        rescaled_eye_mean = rescaled_eye_mean - mean(rescaled_eye_mean);
                        % Define the threshold values
                        peak_threshold = 0.5 * max(abs(rescaled_eye_mean));
                        % Set values below the threshold to zero
                        rescaled_eye_mean(abs(rescaled_eye_mean) <= peak_threshold) = 0;
                          
                        % Initialize arrays to store intermediary metrics per interval
                        nr_eye_movements = zeros(1, floor(length(rescaled_eye_mean) / x_1000));
                        avg_abs_amplitude_change = zeros(1, length(nr_eye_movements));
                        combined_measure_interval = zeros(1, length(nr_eye_movements));
                        
                        % Calculates number of eye movements, average amplitude change, and combined measure per time interval
                        for k = 1:length(nr_eye_movements)
                            % Define the start and end indices of the current interval
                            start_idx = (k - 1) * x_1000 + 1;
                            end_idx = min(k * x_1000, length(rescaled_eye_mean));
                            interval_data = rescaled_eye_mean(start_idx:end_idx);
                            % Find local maxima in the interval data
                            [maxima, max_locs] = findpeaks(interval_data, 'MinPeakHeight', peak_threshold);
                            % Find local minima in the interval data
                            [minima, min_locs] = findpeaks(-interval_data, 'MinPeakHeight', peak_threshold);
                            minima = -minima; % Invert back to positive values
                            
                            if ~isempty(max_locs) | ~isempty(min_locs)
                                % Find zero-crossings in the interval data
                                % Combine peak locations of maxima, minima, and zero crossings
                                zero_crossings = find(diff(sign(interval_data)) == -2);
                                event_locs = sort([max_locs, min_locs, zero_crossings]);
                                % Calculate the time difference between adjacent peaks
                                peak_time_diff = diff(event_locs);
                                
                                % Count the number of peaks with time difference greater than 500 ms
                                if length(event_locs) == 1
                                    nr_eye_movements(k) = 1;
                                else
                                    nr_eye_movements(k) = sum(peak_time_diff > 500) + 1;
                                    % Add 1 for the first peak in the interval
                                end
                                % Calculate the absolute amplitude changes for each peak
                                abs_amplitude_changes = sum(abs(maxima)) + sum(abs(minima));
                                % Calculate the average absolute amplitude change per peak
                                avg_abs_amplitude_change(k) = mean(abs_amplitude_changes);
                                % Calculate the combined measure for the interval
                                combined_measure_interval(k) = nr_eye_movements(k) * avg_abs_amplitude_change(k);
                            end
                            
                        end
                        
                        % Store the final metrics per component
                        total_eye_movements_per_component(c) = sum(nr_eye_movements);
                        avg_amplitude_change_mean_per_component(c) = mean(avg_abs_amplitude_change(avg_abs_amplitude_change > 0));
                        combined_measure_per_component(c, :) = combined_measure_interval;
                    end
                    
                    time_measure =60;
                    % Calculate overall metrics if needed
                    total_eye_movements = sum(total_eye_movements_per_component);
                    avg_amplitude_change_mean = mean(avg_amplitude_change_mean_per_component);
                    avg_eye_movements_per_min = (total_eye_movements / (length(EEGc.data) / 1000)) * time_measure;
                    total_CM_1 = avg_amplitude_change_mean * avg_eye_movements_per_min;
                    total_CM_2 = sum(combined_measure_per_component(:)); 
                    
                    % Depending on how you do the post analysis you could
                    % use CM1 or CM2, CM2 is the accumulated value for the
                    % whole dataset, CM1 is the average combined measure
                    % per minute. 
                end
                
                %% Write data into excell file after each iteration
%               nr_eye_movements, avg_amplitude_change_mean, avg_eye_movements_per_min,
%               total_CM1, total_CM2 & max_eye_score, 
                data1{sub, 1} = fileName_cleaned; % Filename
                data1{sub, 2} = total_eye_movements; % Nr of movements
                data1{sub, 3} = avg_eye_movements_per_min; % Nr of mov per minute
                data1{sub, 4} = avg_amplitude_change_mean; % Intensity/amplitude
                data1{sub, 5} = total_CM_1; % Combinatory measure 1
                data1{sub, 6} = total_CM_2; % Combinatory measure 1
                data1{sub, 7} = max_eye_score; %Highest score for eye from ICLabel

 !               filename_excel2 = '/home/zwergwidder/DonnaH/resting_state_data.xlsx'; %adjust to your desired filename and path
                T = cell2table(data1, 'VariableNames', header);
                if exist(filename_excel2, 'file')
                    % If the file already exists, load it
                    T_existing = readtable(filename_excel2);
                    % Ensure both tables have the same variable names
                    T_existing.Properties.VariableNames = header;
                    T_existing.('Nr of Eye Movements') = num2cell(T_existing.('Nr of Eye Movements'));
                    T_existing.('Nr Eye movements per Minute') = num2cell(T_existing.('Nr Eye movements per Minute'));
                    T_existing.('Average Intensity') = num2cell(T_existing.('Average Intensity'));
                    T_existing.('Combinatory Measure 1') = num2cell(T_existing.('Combinatory Measure 1'));
                    T_existing.('Combinatory Measure 2') = num2cell(T_existing.('Combinatory Measure 2'));
                    T_existing.('Maximum Eye Score') = num2cell(T_existing.('Maximum Eye Score'));
                    
                    % Append new data to the existing table
                    T_updated = [T_existing; T];
                    
                else
                    % If the file doesn't exist, create a new table
                    T_updated = T;
                end
                % Write the updated table to the Excel file
                writetable(T_updated, filename_excel2);
            end
        end
    end
end

