% analyze_coarse_grained_mi_gc.m
%
% This script performs the following tasks:
%
% 1. Loads the Desikan-Killiany atlas region labels and coarse-grains them into larger anatomical regions.
% 2. Loads the source reconstructed time-series data (X) from the preprocessed data directory.
% 3. Concatenates epochs, performs coarse-graining, and re-epochs the data.
% 4. Performs Dickey-Fuller test, demeaning, VAR model building, GC significance testing.
% 5. Stores the results in the 'results/coarse-grained' directory.
%
% Author: Borjan Milinkovic
% Date: 18.09.2024


% Step 1: Load the labels and define coarse-grained regions
load('labels.mat', 'labels');  % Load the Desikan-Killiany labels

% Convert the labels from character matrix to a cell array of strings
labels = cellstr(labels);  % Converts the character matrix into a cell array of strings

% Coarse-graining based on anatomical regions from Desikan-Killiany atlas
groups = { ...
    find(contains(labels, {'superiorfrontal', 'rostralmiddlefrontal', 'caudalmiddlefrontal', ...
                           'parsopercularis', 'parsorbitalis', 'parstriangularis', ...
                           'lateralorbitofrontal', 'medialorbitofrontal', 'frontalpole', ...
                           'precentral'})),  % Frontal Lobe
    find(contains(labels, {'superiorparietal', 'inferiorparietal', 'supramarginal', ...
                           'precuneus', 'postcentral', 'paracentral'})),  % Parietal Lobe
    find(contains(labels, {'superiortemporal', 'middletemporal', 'inferiortemporal', ...
                           'fusiform', 'temporalpole', 'transversetemporal', 'bankssts'})),  % Temporal Lobe
    find(contains(labels, {'lateraloccipital', 'cuneus', 'lingual', 'pericalcarine'})),  % Occipital Lobe
    find(contains(labels, {'insula', 'entorhinal', 'parahippocampal', 'posteriorcingulate', ...
                           'isthmuscingulate', 'caudalanteriorcingulate', 'rostralanteriorcingulate'}))  % Limbic System
};

% Display the new grouping structure
disp('Updated Coarse-Grained Regions:');
for i = 1:length(groups)
    fprintf('Group %d: %s\n', i, strjoin(labels(groups{i}), ', '));
end

% Ensure that the groups are not empty and are within the valid range
for i = 1:length(groups)
    assert(~isempty(groups{i}), sprintf('Group %d is empty', i));
    assert(all(groups{i} >= 1 & groups{i} <= 68), sprintf('Group %d contains out-of-range indices', i));
end

% Step 2: Define directories
preprocessed_data_dir = '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data/';
data_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/';
coarse_grained_dir = '/Users/borjan/code/healthy-sleep-project/code/results/coarse-grained/';

% Create results directory if it doesn't exist
if ~exist(coarse_grained_dir, 'dir')
    mkdir(coarse_grained_dir);
end

% List of participants and sleep stages
files = dir(fullfile(data_dir, 'pwcgc_params_EPCTL*_*.mat'));  % List of all .mat files for participants

% Step 3: Loop through each participant and condition to load the necessary data
for f = 1:length(files)
    % Get the file name and extract participant ID and condition (e.g., 'EPCTL01', 'N1')
    file_name = files(f).name;
    fprintf('Loading file: %s\n', file_name);
    
    % Extract participant ID and sleep stage from the filename
    tokens = regexp(file_name, 'pwcgc_params_(EPCTL\d{2})_(\w+)\.mat', 'tokens');
    participant_id = tokens{1}{1};  % e.g., 'EPCTL01'
    sleep_stage = tokens{1}{2};     % e.g., 'N1', 'N2', etc.

    % Clean up sleep stage variable by removing any extra underscores (if present)
    sleep_stage = strrep(sleep_stage, '_', '');  % This removes any extra underscores
    
    
    % Step 4: Load the source reconstructed time-series data (X) from the corresponding file
    source_data_file = fullfile(preprocessed_data_dir, participant_id, 'mat', ...
        sprintf('%s_%s_source_time_series.mat', participant_id, sleep_stage));
    
    % Load the source reconstructed time series (X)
    load(source_data_file, 'source_ts');  % Load the time series variable 'source_ts'
    X = source_ts;  % Dimensions: 68 x 513 x 74 (regions x time points x epochs)
    
    % Step 5: Concatenate epochs into a single long time series for coarse-graining
    [n_regions, n_timepoints, n_epochs] = size(X);
    X_concat = reshape(X, [n_regions, n_timepoints * n_epochs]);  % Concatenate epochs along time dimension
    
    % Compute the average time-series for each coarse-grained group
    X_coarse_concat = nan(length(groups), size(X_concat, 2));  % Initialize an array for coarse-grained time-series
    for i = 1:length(groups)
        group_indices = groups{i};  % Get indices of regions in the group
        X_coarse_concat(i, :) = mean(X_concat(group_indices, :), 1);  % Average time series over the group regions
    end
    
    % Step 6: Re-epoch the coarse-grained data back into 513 time-point chunks
    X_coarse = reshape(X_coarse_concat, [length(groups), n_timepoints, n_epochs]);

    % Step 7: Perform Dickey-Fuller test to check stationarity for each region and epoch
    fprintf('Performing Dickey-Fuller test for participant %s, sleep stage %s...\n', participant_id, sleep_stage);
    for region = 1:size(X_coarse, 1)
        for epoch = 1:size(X_coarse, 3)
            ts = squeeze(X_coarse(region, :, epoch));  % Extract time series for the current region and epoch
            [h, pValue, stat, cValue] = adftest(ts);  % Perform the Augmented Dickey-Fuller test
            if h == 1
                fprintf('Region %d, Epoch %d: Stationary (p = %.4f)\n', region, epoch, pValue);
            else
                fprintf('Region %d, Epoch %d: Non-stationary (p = %.4f)\n', region, epoch, pValue);
            end
        end
    end
    
    % Step 8: Demean the data if necessary (function 'demean' assumed)
    if exist('demean', 'file')
        X_coarse = demean(X_coarse, true);
    end
    
    % Step 9: Model order selection and VAR model building
    momax = 20;  % Maximum model order
    regmode = 'LWR';  % LWR
    [moaic, mobic] = tsdata_to_varmo(X_coarse, momax, regmode, [], [], []);
    morder = moselect(sprintf('VAR model order selection, max = %d', momax), 'AIC', 'AIC', moaic, 'BIC', mobic);
    [A, V] = tsdata_to_var(X_coarse, morder, regmode);
    fprintf('VAR model order selected: %d\n', morder);
    
    % Step 10: Pairwise Granger Causality and significance testing
    F = var_to_pwcgc(A, V);
    test_stat = 'F';  % F-statistic
    tstat = var_to_pwcgc_tstat(X_coarse, V, morder, regmode, test_stat);

    % p-value and significance testing
    alpha = 0.05;
    pval = mvgc_pval(tstat, test_stat, 1, 1, size(A,1)-2, morder, size(X_coarse,2), size(X_coarse,3));
    sig = significance(pval, alpha, 'FDR');  % FDR correction

    % Step 11: Save results
    eweights = F / nanmax(F(:));  % Normalize Granger Causality
    eweights_file = fullfile(coarse_grained_dir, sprintf('coarse_eweights_matrix_%s_%s.mat', participant_id, sleep_stage));
    save(eweights_file, 'eweights');
    
    % Convert to SS parameters for Directed Influence (DI) analysis
    [Ass, Vss] = transform_var(A, V);
    [A, C, K] = var_to_ss(Ass);
    CAK = Ass;
    fres = 256;
    H = var2trfun(Ass, fres);

    % Save the model parameters
    modfile = fullfile(coarse_grained_dir, sprintf('coarse_pwcgc_params_%s_%s.mat', participant_id, sleep_stage));
    save(modfile, 'X_coarse', 'V', 'CAK', 'H', 'Vss');
    fprintf('Saved Granger Causality model for %s, sleep stage %s.\n', participant_id, sleep_stage);
end
