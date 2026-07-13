% analyze_gwgc_mannwhitney_in_out_timedomain.m
%
% This script performs the following tasks:
%
% 1. Loads the Granger causality (GC) results for each participant and each sleep stage (N1, N2, N3, R)
%    and the wake condition (W_) from the specified directory.
%
% 2. The transfer function (H) is computed to calculate spectral Granger causality across the entire
%    frequency range for both inbound and outbound GC. The GC values are integrated over specific frequency bands.
%
% 3. **Time-domain Granger causality** is calculated by averaging the frequency-domain GC across all frequencies.
%
% 4. Prepares the data for statistical comparison (Mann-Whitney U test) between the wake condition (W_)
%    and other sleep stages for each brain region (variable) across participants for both time-domain and frequency-domain GC.
%
% 5. For each region, condition comparison, and frequency band, the script computes the Mann-Whitney U test
%    and calculates the effect size for both time-domain and frequency-domain GC.
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
%   - Frequency-domain GC Mann-Whitney U test results saved in separate files for inbound and outbound GC.
%   - Time-domain GC Mann-Whitney U test results saved in separate files for inbound and outbound GC.
%
% Author: Borjan Milinkovic
% Date: 17.09.2024

%% Step 1: Define parameters
data_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/';

% Conditions: sleep stages and wake
conditions = {'N1', 'N2', 'N3', 'W_', 'R_'};

% Frequency bands and their ranges (in Hz)
bands = {'delta', 'theta', 'alpha', 'sigma', 'beta'};
freq_bands = [1 4; 4 8; 8 12; 12 16; 16 30];  % [low, high] for each band

% Sampling rate (assumed to be 256 Hz)
fs = 256;

% Step 2: Initialize storage for Granger causality data and Mann-Whitney U results
gc_in_data = struct();  % Store inbound Granger causality for each participant and condition (frequency-domain)
gc_out_data = struct();  % Store outbound Granger causality for each participant and condition (frequency-domain)
gc_in_time_data = struct();  % Store inbound time-domain Granger causality (averaged frequency-domain)
gc_out_time_data = struct();  % Store outbound time-domain Granger causality (averaged frequency-domain)
mw_in_effect_sizes = struct();  % Store Mann-Whitney U test results and effect sizes for inbound GC (frequency-domain)
mw_out_effect_sizes = struct();  % Store Mann-Whitney U test results and effect sizes for outbound GC (frequency-domain)
mw_in_time_effect_sizes = struct();  % Store Mann-Whitney U test results for inbound time-domain GC
mw_out_time_effect_sizes = struct();  % Store Mann-Whitney U test results for outbound time-domain GC

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

    % Step 4: Compute Granger causality for both inbound and outbound GC
    gc_in_all_freqs = var_to_siogc(data.CAK, data.V, 'in', fres, H);
    gc_out_all_freqs = var_to_siogc(data.CAK, data.V, 'out', fres, H);

    % Step 5: Compute **time-domain Granger causality** by averaging the frequency-domain GC across all frequencies
    gc_in_time = mean(gc_in_all_freqs, 2);  % Time-domain inbound GC (mean across all frequencies)
    gc_out_time = mean(gc_out_all_freqs, 2);  % Time-domain outbound GC (mean across all frequencies)

    % Store the time-domain GC in the gc_in_time_data and gc_out_time_data structures
    gc_in_time_data.(participant_id).(condition) = gc_in_time;  % Time-domain inbound GC
    gc_out_time_data.(participant_id).(condition) = gc_out_time;  % Time-domain outbound GC

    % Extract the Granger causality for each frequency band for both in and out (frequency-domain)
    for b = 1:length(bands)
        band_name = bands{b};
        freq_range = freq_bands(b, :);

        % Get indices corresponding to the frequency range
        freq_idx = find(freqs >= freq_range(1) & freqs <= freq_range(2));

        % Extract GC values for the specified frequency range
        gc_in_band = gc_in_all_freqs(:, freq_idx);
        gc_out_band = gc_out_all_freqs(:, freq_idx);

        % Integrate the GC values over the frequency range (mean over frequencies)
        gc_in_band_integrated = mean(gc_in_band, 2);  % For inbound GC
        gc_out_band_integrated = mean(gc_out_band, 2);  % For outbound GC

        % Store the integrated GC values in the gc_data structures as cell arrays
        if ~isfield(gc_in_data, participant_id)
            gc_in_data.(participant_id) = struct();
        end
        if ~isfield(gc_out_data, participant_id)
            gc_out_data.(participant_id) = struct();
        end

        gc_in_data.(participant_id).(condition).(band_name) = num2cell(gc_in_band_integrated);  % Inbound
        gc_out_data.(participant_id).(condition).(band_name) = num2cell(gc_out_band_integrated);  % Outbound
    end
end

% Save GC data for future use (time-domain and frequency-domain)
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/';
save(fullfile(results_dir, 'gcin_freq_subject_data_struct.mat'), 'gc_in_data');
save(fullfile(results_dir, 'gcout_freq_subject_data_struct.mat'), 'gc_out_data');
save(fullfile(results_dir, 'gcin_time_subject_data_struct.mat'), 'gc_in_time_data');
save(fullfile(results_dir, 'gcout_time_subject_data_struct.mat'), 'gc_out_time_data');

%% Step 6: Perform Mann-Whitney U test for both time-domain and frequency-domain GC
participants = fieldnames(gc_in_data);  % List of participants

for c = 1:length(conditions)
    condition = conditions{c};
    if strcmp(condition, 'W_'), continue; end  % Skip wake as a condition

    for b = 1:length(bands)
        band_name = bands{b};

        % Get the number of brain regions (assuming all participants have the same number)
        num_regions = length(gc_in_data.(participants{1}).W_.(band_name));  % Number of regions

        % Initialize arrays for storing GC values across subjects for each region
        group1_in = cell(1, num_regions);  % For W_ (wake) data for inbound GC
        group2_in = cell(1, num_regions);  % For sleep condition for inbound GC
        group1_out = cell(1, num_regions);  % For W_ (wake) data for outbound GC
        group2_out = cell(1, num_regions);  % For sleep condition for outbound GC

        % Initialize arrays for time-domain data
        group1_in_time = cell(1, num_regions);  % Inbound time-domain GC (wake)
        group2_in_time = cell(1, num_regions);  % Inbound time-domain GC (sleep stage)
        group1_out_time = cell(1, num_regions);  % Outbound time-domain GC (wake)
        group2_out_time = cell(1, num_regions);  % Outbound time-domain GC (sleep stage)

        % Loop through each region and collect data across participants (frequency-domain)
        for r = 1:num_regions
            group1_in{r} = cell(1, length(participants));  % Wake GC values (inbound)
            group2_in{r} = cell(1, length(participants));  % Sleep GC values (inbound)
            group1_out{r} = cell(1, length(participants));  % Wake GC values (outbound)
            group2_out{r} = cell(1, length(participants));  % Sleep GC values (outbound)

            for p = 1:length(participants)
                participant_id = participants{p};

                % Extract GC values for the specified region and band (frequency-domain)
                wake_gc_in = gc_in_data.(participant_id).W_.(band_name){r};  % Inbound GC for wake
                cond_gc_in = gc_in_data.(participant_id).(condition).(band_name){r};  % Inbound GC for sleep stage
                wake_gc_out = gc_out_data.(participant_id).W_.(band_name){r};  % Outbound GC for wake
                cond_gc_out = gc_out_data.(participant_id).(condition).(band_name){r};  % Outbound GC for sleep stage

                % Store frequency-domain GC
                group1_in{r}{p} = wake_gc_in;  % Store inbound wake GC
                group2_in{r}{p} = cond_gc_in;  % Store inbound sleep stage GC
                group1_out{r}{p} = wake_gc_out;  % Store outbound wake GC
                group2_out{r}{p} = cond_gc_out;  % Store outbound sleep stage GC

                % Store time-domain GC in cell arrays for Mann-Whitney test compatibility
                group1_in_time{r}{p} = gc_in_time_data.(participant_id).W_(r);
                group2_in_time{r}{p} = gc_in_time_data.(participant_id).(condition)(r);
                group1_out_time{r}{p} = gc_out_time_data.(participant_id).W_(r);
                group2_out_time{r}{p} = gc_out_time_data.(participant_id).(condition)(r);
            end

            % Perform Mann-Whitney U test for both time-domain and frequency-domain GC
            if length(group1_in{r}) == length(group2_in{r})
                [z_in, U_in] = mann_whitney_group(group1_in{r}, group2_in{r});
                p_value_in = erfc(abs(z_in) / sqrt(2));  % Two-tailed p-value
                mw_in_effect_sizes.(condition).(band_name)(r).z = z_in;
                mw_in_effect_sizes.(condition).(band_name)(r).U = U_in;
                mw_in_effect_sizes.(condition).(band_name)(r).p_value = p_value_in;
            end

            if length(group1_out{r}) == length(group2_out{r})
                [z_out, U_out] = mann_whitney_group(group1_out{r}, group2_out{r});
                p_value_out = erfc(abs(z_out) / sqrt(2));  % Two-tailed p-value
                mw_out_effect_sizes.(condition).(band_name)(r).z = z_out;
                mw_out_effect_sizes.(condition).(band_name)(r).U = U_out;
                mw_out_effect_sizes.(condition).(band_name)(r).p_value = p_value_out;
            end

            % Time-domain Mann-Whitney U test
            if length(group1_in_time{r}) == length(group2_in_time{r})
                [z_in_time, U_in_time] = mann_whitney_group(group1_in_time{r}, group2_in_time{r});
                mw_in_time_effect_sizes.(condition)(r).z = z_in_time;
                mw_in_time_effect_sizes.(condition)(r).U = U_in_time;
                mw_in_time_effect_sizes.(condition)(r).p_value = erfc(abs(z_in_time) / sqrt(2));
            end

            if length(group1_out_time{r}) == length(group2_out_time{r})
                [z_out_time, U_out_time] = mann_whitney_group(group1_out_time{r}, group2_out_time{r});
                mw_out_time_effect_sizes.(condition)(r).z = z_out_time;
                mw_out_time_effect_sizes.(condition)(r).U = U_out_time;
                mw_out_time_effect_sizes.(condition)(r).p_value = erfc(abs(z_out_time) / sqrt(2));
            end
        end
    end
end

%% Step 7: Save the Mann-Whitney U test results for time-domain and frequency-domain GC
save(fullfile(results_dir, 'gc_in_freq_mw_effect_sizes.mat'), 'mw_in_effect_sizes');  % Save inbound frequency-domain GC results
save(fullfile(results_dir, 'gc_out_freq_mw_effect_sizes.mat'), 'mw_out_effect_sizes');  % Save outbound frequency-domain GC results
save(fullfile(results_dir, 'gc_in_time_mw_effect_sizes.mat'), 'mw_in_time_effect_sizes');  % Save inbound time-domain GC results
save(fullfile(results_dir, 'gc_out_time_mw_effect_sizes.mat'), 'mw_out_time_effect_sizes');  % Save outbound time-domain GC results
