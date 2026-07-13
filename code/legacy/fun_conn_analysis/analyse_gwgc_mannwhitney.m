% analyze_gwgc_mannwhitney.m
%
% This script performs the following tasks:
%
% 1. Loads the Granger causality (GC) results for each participant and each sleep stage (N1, N2, N3, R) 
%    and the wake condition (W_) from the specified directory.
%
% 2. The transfer function (H) is computed to calculate spectral Granger causality across the entire 
%    frequency range. These GC values are then integrated over specific frequency bands.
%
% 3. Prepares the data for statistical comparison (Mann-Whitney U test) between the wake condition (W_) 
%    and other sleep stages for each brain region (variable) across participants.
%
% 4. For each region, condition comparison, and frequency band, the script computes the Mann-Whitney U test
%    and calculates the effect size.
%
% Inputs:
%   - Directory containing the .mat files with GC parameters for each participant and condition.
%   - The following frequency bands are analyzed:
%       * Delta: 1-4 Hz
%       * Theta: 4-8 Hz
%       * Alpha: 8-13 Hz
%       * Beta: 15-30 Hz
%       * Gamma: 30-45 Hz
%
% Outputs:
%   - Mann-Whitney U test results (Z-scores, U-statistics, p-values) and effect sizes saved in 'mw_effect_sizes.mat'.
%
% Example Usage:
%   This script compares the wake condition (W_) against each sleep stage for each frequency band and region.
%
% Author: Borjan Milinkovic
% Date: 17.09.2024


%% Step 1: Define parameters
data_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/';

% Conditions: sleep stages and wake
conditions = {'N1', 'N2', 'N3', 'W_', 'R_'};

% Frequency bands and their ranges (in Hz)
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
freq_bands = [1 4; 4 8; 8 13; 15 30; 30 45];  % [low, high] for each band

% Sampling rate (assumed to be 256 Hz)
fs = 256;

% Step 2: Initialize storage for Granger causality data and Mann-Whitney U results
gc_data = struct();  % Store Granger causality for each participant and condition
mw_effect_sizes = struct();  % Store Mann-Whitney U test results and effect sizes

%% Step 3: Load files for all participants and all conditions
files = dir(fullfile(data_dir, 'pwcgc_params_EPCTL*_*.mat'));  % List of .mat files for participants

for f = 1:length(files)
    % Extract file name, participant ID, and condition
    file_name = files(f).name;
    fprintf('Loading file: %s\n', file_name);

    % Extract participant ID and condition (e.g., 'EPCTL01', 'N1')
    tokens = regexp(file_name, 'pwcgc_params_(EPCTL\d{2})_(\w+)\.mat', 'tokens');
    participant_id = tokens{1}{1};
    condition = tokens{1}{2};  % e.g., 'N1', 'N2', 'W_', 'R_'

    % Load the parameters (V, CAK)
    data = load(fullfile(data_dir, file_name), 'V', 'CAK', 'Vss');

    % Frequency resolution for computing Granger causality
    fres = 256;

    % Compute the transfer function (H) for the full frequency range
    H = var2trfun(data.CAK, fres);  % Compute transfer function

    % Get the corresponding frequency vector
    freqs = sfreqs(fres, fs);

    % Step 4: Compute Granger causality over all frequencies
    gc_all_freqs = var_to_siogc(data.CAK, data.V, 'out', fres, H);

    % Extract the Granger causality for each frequency band
    for b = 1:length(bands)
        band_name = bands{b};
        freq_range = freq_bands(b, :);

        % Get indices corresponding to the frequency range
        freq_idx = find(freqs >= freq_range(1) & freqs <= freq_range(2));

        % Extract GC values for the specified frequency range
        gc_band = gc_all_freqs(:, freq_idx);

        % Step 5: Integrate the GC values over the frequency range (mean over frequencies)
        gc_band_integrated = mean(gc_band, 2);  % Take the mean across frequencies (for each region)

        % Step 6: Store the integrated GC values in the gc_data structure as a cell array
        if ~isfield(gc_data, participant_id)
            gc_data.(participant_id) = struct();
        end
        gc_data.(participant_id).(condition).(band_name) = num2cell(gc_band_integrated);  % Store as a cell array
    end
end

% Save GC data for future use
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/';
save(fullfile(results_dir, 'gcout_subject_data_struct.mat'), 'gc_data');

%% Step 7: Perform Mann-Whitney U test for each region (variable)
participants = fieldnames(gc_data);  % List of participants

for c = 1:length(conditions)
    condition = conditions{c};
    if strcmp(condition, 'W_'), continue; end  % Skip wake as a condition

    for b = 1:length(bands)
        band_name = bands{b};

        % Get the number of brain regions (assuming all participants have the same number)
        num_regions = length(gc_data.(participants{1}).W_.(band_name));  % Number of regions is the number of cells

        % Initialize arrays for storing GC values across subjects for each region
        group1 = cell(1, num_regions);  % For W_ (wake) data
        group2 = cell(1, num_regions);  % For sleep condition (e.g., N1, N2, etc.)

        % Loop through each region and collect data across participants
        for r = 1:num_regions
            group1{r} = cell(1, length(participants));  % Ensure it's a cell array for wake GC values
            group2{r} = cell(1, length(participants));  % Ensure it's a cell array for sleep stage GC values

            for p = 1:length(participants)
                participant_id = participants{p};

                % Extract GC values for the specified region and band (already integrated over frequencies)
                wake_gc = gc_data.(participant_id).W_.(band_name){r};  % GC for wake (cell value)
                cond_gc = gc_data.(participant_id).(condition).(band_name){r};  % GC for sleep stage (cell value)

                % Store the GC values for this region and participant into the corresponding cell array
                group1{r}{p} = wake_gc;  % Store participant's wake GC value in group1
                group2{r}{p} = cond_gc;  % Store participant's sleep stage GC value in group2

                % Debugging: Print the data being appended to each group
                fprintf('Participant: %s, Region: %d, Wake GC: %.4f, Condition GC: %.4f\n', ...
                        participant_id, r, wake_gc, cond_gc);
            end

            % Debugging: Print the sizes of group1 and group2 for each region
            fprintf('Region %d, Wake Group Size: %d, %s Group Size: %d\n', ...
                r, length(group1{r}), condition, length(group2{r}));
        end

        % Step 8: Perform Mann-Whitney U test for each region
        for r = 1:num_regions
            % Check if both groups have the same number of samples before running the test
            if length(group1{r}) == length(group2{r})
                % Run the Mann-Whitney U test for the GC values of each region
                [z, U] = mann_whitney_group(group1{r}, group2{r});

                % Calculate the p-value and effect size
                p_value = erfc(abs(z) / sqrt(2));  % Two-tailed p-value

                % Step 9: Store the results for the region
                mw_effect_sizes.(condition).(band_name)(r).z = z;
                mw_effect_sizes.(condition).(band_name)(r).U = U;
                mw_effect_sizes.(condition).(band_name)(r).p_value = p_value;

                fprintf('Region %d: W_ vs. %s (%s band): Z=%.4f, U=%.4f, p=%.4f\n', ...
                        r, condition, band_name, z, U, p_value);
            else
                % If there is a mismatch, report it
                fprintf('Mismatch in group sizes for region %d, skipping Mann-Whitney U test.\n', r);
            end
        end
    end
end

% Step 10: Save the Mann-Whitney U test results and effect sizes
save(fullfile(results_dir, 'gc_out_mw_effect_sizes.mat'), 'mw_effect_sizes');
