% analyze_mi_spectralmi.m
%
% This script performs the following tasks:
%
% 1. Loads the Granger causality (GC) and VAR results for each participant and condition (N1, N2, N3, R, W_).
%
% 2. Computes the Mutual Information (MI) between each region (variable) and the rest of the system
%    using the covariance matrix from the VAR model.
%
% 3. Computes the Spectral Conditional Mutual Information (CMI) across the frequency domain
%    for each region using the cross-power spectral density (CPSD) matrix.
%
% 4. Stores and saves the MI and CMI results for later use.
%
% Inputs:
%   - Directory containing the .mat files with VAR parameters (e.g., covariance matrices).
%   - The following frequency bands are analyzed:
%       * Delta: 1-4 Hz
%       * Theta: 4-8 Hz
%       * Alpha: 8-13 Hz
%       * Beta: 15-30 Hz
%       * Gamma: 30-45 Hz
%
% Outputs:
%   - MI results (time domain) and CMI results (frequency domain) for each region and condition.
%   - Saved MI and CMI values for later analysis.
%
% Example Usage:
%   This script computes MI and CMI across conditions and frequency bands for each region.
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

% Step 2: Initialize storage for MI, CMI, and results
mi_data = struct();  % Store Mutual Information (MI) for each participant and condition
cmi_data = struct();  % Store Spectral Conditional Mutual Information (CMI)
freq_cmi_data = struct();  % Store CMI for each frequency band

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

    % Step 4: Compute Mutual Information (MI) in the time domain using covariance matrix V
    mi_values = cov_to_iomi(data.V);  % Compute MI between each variable and the rest of the system

    % Step 5: Store the MI values in the mi_data structure
    if ~isfield(mi_data, participant_id)
        mi_data.(participant_id) = struct();
    end
    mi_data.(participant_id).(condition) = mi_values;

    % Frequency resolution for computing CMI
    fres = 256;

    % Step 6: Compute the transfer function H and cross-power spectral density (CPSD)
    H = var2trfun(data.CAK, fres);  % Compute the transfer function for the VAR model
    S = var_to_cpsd(H, data.V, fres);  % Compute the CPSD matrix

    % Step 7: Compute Spectral Conditional Mutual Information (CMI) using CPSD
    cmi_values = cpsd_to_scwiomi(S);  % CMI over the entire frequency domain

    % Store the spectral CMI data
    if ~isfield(cmi_data, participant_id)
        cmi_data.(participant_id) = struct();
    end
    cmi_data.(participant_id).(condition) = cmi_values;

    % Step 8: Extract and integrate CMI over specific frequency bands
    freqs = sfreqs(fres, fs);  % Get the frequency vector corresponding to the CPSD

    for b = 1:length(bands)
        band_name = bands{b};
        freq_range = freq_bands(b, :);

        % Get indices corresponding to the frequency range
        freq_idx = find(freqs >= freq_range(1) & freqs <= freq_range(2));

        % Extract and integrate CMI over the specified frequency range
        cmi_band_integrated = mean(cmi_values(:, freq_idx), 2);  % Take the mean over frequencies for each region

        % Step 9: Store the integrated CMI values in freq_cmi_data
        if ~isfield(freq_cmi_data, participant_id)
            freq_cmi_data.(participant_id) = struct();
        end
        freq_cmi_data.(participant_id).(condition).(band_name) = cmi_band_integrated;
    end
end

% Save MI and CMI data for future use
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/';
save(fullfile(results_dir, 'mi_data.mat'), 'mi_data');  % Save Mutual Information data
save(fullfile(results_dir, 'cmi_data.mat'), 'cmi_data');  % Save Spectral Conditional MI data
save(fullfile(results_dir, 'freq_cmi_data.mat'), 'freq_cmi_data');  % Save integrated CMI by frequency bands

