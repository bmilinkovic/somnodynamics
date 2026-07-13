% analyze_mi_mannwhitney.m
%
% This script performs the following tasks:
%
% 1. Loads the Mutual Information (MI) results for each participant and each sleep stage (N1, N2, N3, R)
%    and the wake condition (W_).
%
% 2. Prepares the data for statistical comparison (Mann-Whitney U test) between the wake condition (W_)
%    and other sleep stages for each brain region (variable) across participants.
%
% 3. For each region, condition comparison, and frequency band, the script computes the Mann-Whitney U test
%    and calculates the effect size.
%
% Inputs:
%   - Directory containing the MI and CMI data files.
%   - The following frequency bands are analyzed:
%       * Delta: 1-4 Hz
%       * Theta: 4-8 Hz
%       * Alpha: 8-12 Hz
%       * Sigma: 12-16 Hz
%       * Beta: 15-30 Hz
%
% Outputs:
%   - Time-domain MI effect sizes saved in 'mi_time_effect_sizes.mat'.
%   - Frequency-domain CMI effect sizes saved in 'mi_freq_effect_sizes.mat'.
%
% Author: Borjan Milinkovic
% Date: 17.09.2024

%% Step 1: Define parameters
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/';

% Conditions: sleep stages and wake
conditions = {'N1', 'N2', 'N3', 'W_', 'R_'};

% Frequency bands
bands = {'delta', 'theta', 'alpha', 'sigma', 'beta'};

% Step 2: Load MI and CMI data
load(fullfile(results_dir, 'mi_data.mat'), 'mi_data');  % Time-domain MI data
load(fullfile(results_dir, 'freq_cmi_data.mat'), 'freq_cmi_data');  % Frequency-domain CMI data

% Step 3: Initialize storage for Mann-Whitney U test results
mi_time_effect_sizes = struct();  % Store the time-domain MI results
mi_freq_effect_sizes = struct();  % Store the frequency-domain CMI results

%% Step 4: Perform Mann-Whitney U test on time-domain MI values
participants = fieldnames(mi_data);  % List of participants

for c = 1:length(conditions)
    condition = conditions{c};
    if strcmp(condition, 'W_'), continue; end  % Skip wake as a condition

    % Number of regions based on one participant (assuming all have the same number of regions)
    num_regions = length(mi_data.(participants{1}).W_);  % Assuming MI values are available for each region

    % Initialize arrays for storing MI values across participants for each region
    group1 = cell(1, num_regions);  % For W_ (wake) data
    group2 = cell(1, num_regions);  % For sleep condition (e.g., N1, N2, etc.)

    % Loop through each region and collect data across participants
    for r = 1:num_regions
        group1{r} = cell(1, length(participants));  % Wake MI values
        group2{r} = cell(1, length(participants));  % Sleep stage MI values

        for p = 1:length(participants)
            participant_id = participants{p};

            % Extract MI values for the wake and sleep stage conditions
            wake_mi = mi_data.(participant_id).W_(r);  % MI for wake
            cond_mi = mi_data.(participant_id).(condition)(r);  % MI for sleep stage

            % Store the MI values for this region and participant
            group1{r}{p} = wake_mi;  % Store wake MI
            group2{r}{p} = cond_mi;  % Store sleep stage MI

            % Debugging: Print the data being appended to each group
            fprintf('Participant: %s, Region: %d, Wake MI: %.4f, Condition MI: %.4f\n', ...
                    participant_id, r, wake_mi, cond_mi);
        end

        % Debugging: Print the sizes of group1 and group2 for each region
        fprintf('Region %d, Wake Group Size: %d, %s Group Size: %d\n', ...
            r, length(group1{r}), condition, length(group2{r}));
    end

    % Step 5: Perform Mann-Whitney U test for each region
    for r = 1:num_regions
        if length(group1{r}) == length(group2{r})
            % Perform Mann-Whitney U test on the MI values for this region
            [z, U] = mann_whitney_group(group1{r}, group2{r});

            % Calculate the p-value and effect size
            p_value = erfc(abs(z) / sqrt(2));  % Two-tailed p-value

            % Step 6: Store the results for the region, converting `r` to a string
            region_field = sprintf('Region_%d', r);  % Convert `r` to a valid field name
            mi_time_effect_sizes.(condition).(region_field).z = z;
            mi_time_effect_sizes.(condition).(region_field).U = U;
            mi_time_effect_sizes.(condition).(region_field).p_value = p_value;

            fprintf('Region %d: W_ vs. %s (Time-domain MI): Z=%.4f, U=%.4f, p=%.4f\n', ...
                    r, condition, z, U, p_value);
        else
            fprintf('Mismatch in group sizes for region %d, skipping Mann-Whitney U test.\n', r);
        end
    end
end

%% Step 6: Perform Mann-Whitney U test on frequency-domain CMI values for each band
for c = 1:length(conditions)
    condition = conditions{c};
    if strcmp(condition, 'W_'), continue; end  % Skip wake as a condition

    for b = 1:length(bands)
        band_name = bands{b};

        % Number of regions based on one participant (assuming all have the same number of regions)
        num_regions = length(freq_cmi_data.(participants{1}).W_.(band_name));  % Assuming CMI values are available for each region

        % Initialize arrays for storing CMI values across participants for each region
        group1 = cell(1, num_regions);  % For W_ (wake) data
        group2 = cell(1, num_regions);  % For sleep condition (e.g., N1, N2, etc.)

        % Loop through each region and collect data across participants
        for r = 1:num_regions
            group1{r} = cell(1, length(participants));  % Wake CMI values
            group2{r} = cell(1, length(participants));  % Sleep stage CMI values

            for p = 1:length(participants)
                participant_id = participants{p};

                % Extract CMI values for the wake and sleep stage conditions
                wake_cmi = freq_cmi_data.(participant_id).W_.(band_name)(r);  % CMI for wake
                cond_cmi = freq_cmi_data.(participant_id).(condition).(band_name)(r);  % CMI for sleep stage

                % Store the CMI values for this region and participant
                group1{r}{p} = wake_cmi;  % Store wake CMI
                group2{r}{p} = cond_cmi;  % Store sleep stage CMI
            end
        end

        % Step 7: Perform Mann-Whitney U test for each region
        for r = 1:num_regions
            if length(group1{r}) == length(group2{r})
                % Perform Mann-Whitney U test on the CMI values for this region and frequency band
                [z, U] = mann_whitney_group(group1{r}, group2{r});

                % Calculate the p-value and effect size
                p_value = erfc(abs(z) / sqrt(2));  % Two-tailed p-value

                % Step 8: Store the results for the region and frequency band, converting `r` to a string
                region_field = sprintf('Region_%d', r);  % Convert `r` to a valid field name
                mi_freq_effect_sizes.(condition).(band_name).(region_field).z = z;
                mi_freq_effect_sizes.(condition).(band_name).(region_field).U = U;
                mi_freq_effect_sizes.(condition).(band_name).(region_field).p_value = p_value;

                fprintf('Region %d: W_ vs. %s (%s CMI): Z=%.4f, U=%.4f, p=%.4f\n', ...
                        r, condition, band_name, z, U, p_value);
            else
                fprintf('Mismatch in group sizes for region %d, skipping Mann-Whitney U test.\n', r);
            end
        end
    end
end

%% Step 9: Save the Mann-Whitney U test results for time-domain MI and frequency-domain CMI
save(fullfile(results_dir, 'mi_time_effect_sizes.mat'), 'mi_time_effect_sizes');  % Save time-domain MI results
save(fullfile(results_dir, 'mi_freq_effect_sizes.mat'), 'mi_freq_effect_sizes');  % Save frequency-domain CMI results
