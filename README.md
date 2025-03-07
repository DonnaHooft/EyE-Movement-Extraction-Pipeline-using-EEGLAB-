# EyE Movement Extraction Pipeline

## Overview
The EyE Movement Extraction Pipeline is a MATLAB-based tool for extracting eye movement data from EEG recordings in fMRI-compatible studies. It is designed for use with BrainProducts EEG caps and utilizes EEGLAB plugins.

This pipeline was developed for the study: **"Simultaneous EEG-fMRI study in healthy humans during induction of propofol anesthesia to investigate the dynamics of thalamocortical functional connectivity in the alpha frequency."**

## Prerequisites
- MATLAB (2019b or later)
- EEGLAB plugin (included in package)
- Minimum 32GB RAM recommended
- Sufficient storage for .eeg files
- Access to Webdisk (Klifoanae) recommended

## Package Contents
- **Manual (PDF)**
- **MATLAB Scripts**:
  - EyE_pipeline_RS.mat (Resting state data analysis)
  - EyE_pipeline_FLCK.mat (Flicker data analysis)
  - condition_order_excell_file_creator.mat
  - excell_combiner_interval_data.mat
  - general_data_analysis.mat
  - flicker_analysis.mat
  - within_subj_flicker_analysis.mat
- **Excel Examples**:
  - combined_eye_stim_freq_example.xlsx
  - flicker_data_example.xlsx
  - resting_state_data_example.xlsx
  - subject_condition_data_example.xlsx
- **MATLAB Plugins** (CWRegrTool and EEGLAB v2023.0)

## Pipeline Execution
1. Ensure correct folder structure:
   - Main folder contains EEG files and subject information.
   - Required files: `subject_information.mat` and `average_eye_kernel.mat`.
2. Adjust script parameters (marked with `!` in the code).
3. Run:
   - `EyE_pipeline_RS.mat` for resting-state data.
   - `EyE_pipeline_FLCK.mat` for flicker data.

## Output Files
### **EyE_pipeline_RS.mat**
- `resting_state_data.xlsx`
  - Total number of eye movements
  - Average movements per minute
  - Amplitude changes & combinatory measures
  - Maximum eye score from ICLabel

### **EyE_pipeline_FLCK.mat**
- `flicker_data.xlsx` (similar to above)
- `eye_mov_per_interval.xlsx` (eye movements per 20s interval)
- `eeg_times.mat` (stimulation frequency change times)
- `mri_volumes.mat` (corresponding MRI volume numbers)

## Additional Scripts
- `condition_order_excell_file_creator.mat`: Creates stimulation frequency protocol files.
- `excell_combiner_interval_data.mat`: Merges interval-based data across subjects.
- `general_data_analysis.mat`: Performs Friedman’s test and post-hoc Dunn test for eye movement differences.
- `flicker_analysis.mat`: Generates preliminary scatter plots of eye movement data.
- `within_subj_flicker_analysis.mat`: Analyzes stimulation frequency effects on eye movement.

## Acknowledgments
Developed under the supervision of **Prof. Afra Wohlschläger, Anna Gehrig, and Juliana Zimmermann** at **Technical University of Munich, Department of Neuroradiology**.

For further inquiries, consult:
- EEGLAB Documentation: [https://sccn.ucsd.edu/eeglab/](https://sccn.ucsd.edu/eeglab/)
- MATLAB Helpdesk: [https://www.mathworks.com/help/](https://www.mathworks.com/help/)

