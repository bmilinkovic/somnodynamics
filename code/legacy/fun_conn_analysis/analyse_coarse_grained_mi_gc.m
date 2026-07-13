% analyze_cross_subject_differences.m
%
% This script performs the following tasks:
%
% 1. Loads the precomputed coarse-grained time-series data.
% 2. Extracts Granger causality (GC) and conditional mutual information (CMI) values for each participant.
% 3. Performs Mann-Whitney group testing between wake and each sleep stage (N1, N2, N3, REM).
% 4. Applies FDR correction for multiple comparisons.
% 5. Plots effect sizes (z-scores) and marks significant differences with white dots.
%
% Author: Borjan Milinkovic
% Date: 18.09.2024

% Step 1: Define directories for the coarse-grained data
coarse_grained_dir = '/Users/borjan/code/healthy-sleep-project/code/results/coarse-grained/';
results_dir = '/Users/borjan/code/healthy-sleep-project/code/results/cross_subject_results/';
data_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/';

% Create results directory if it doesn't exist
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% List of coarse-grained data files
files = dir(fullfile(coarse_grained_dir, 'coarse_pwcgc_params_EPCTL*_*.mat'));

% Initialize storage for GC and MI results, grouped by condition
gc_wake = [];
gc_sleep_N1 = [];
gc_sleep_N2 = [];
gc_sleep_N3 = [];
gc_sleep_R = [];

cmi_wake = [];
cmi_sleep_N1 = [];
cmi_sleep_N2 = [];
cmi_sleep_N3 = [];
cmi_sleep_R = [];

% Define sleep stages of interest
sleep_stages = {'N1', 'N2', 'N3', 'R'};

% Step 2: Loop through each participant and condition to load the necessary data
for f = 1:length(files)
    % Get the file name and extract participant ID and condition (e.g., 'EPCTL01', 'N1')
    file_name = files(f).name;
    fprintf('Loading file: %s\n', file_name);
    
    % Extract participant ID and sleep stage from the filename
    tokens = regexp(file_name, 'coarse_pwcgc_params_(EPCTL\d{2})_(\w+)\.mat', 'tokens');
    participant_id = tokens{1}{1};  % e.g., 'EPCTL01'
    condition = tokens{1}{2};       % e.g., 'N1', 'W', etc.
    
    % Step 3: Load the precomputed coarse-grained data (GC and CMI)
    coarse_data_file = fullfile(coarse_grained_dir, file_name);
    load(coarse_data_file, 'V', 'CAK');  % Load the covariance matrix and VAR coefficients
    
    % Step 4: Compute Granger causality and conditional mutual information (CMI)
    F = var_to_pwcgc(CAK, V);  % Compute Granger causality
    cmi = cov_to_pwcmi(V);     % Compute conditional mutual information
    
    % Step 5: Group the results by condition
    if strcmp(condition, 'W')  % Wake condition
        gc_wake = cat(3, gc_wake, F);  % Append GC results for wake
        cmi_wake = cat(3, cmi_wake, cmi); % Append CMI results for wake
    elseif strcmp(condition, 'N1')  % Sleep stage N1
        gc_sleep_N1 = cat(3, gc_sleep_N1, F);  % Append GC results for N1
        cmi_sleep_N1 = cat(3, cmi_sleep_N1, cmi); % Append CMI results for N1
    elseif strcmp(condition, 'N2')  % Sleep stage N2
        gc_sleep_N2 = cat(3, gc_sleep_N2, F);  % Append GC results for N2
        cmi_sleep_N2 = cat(3, cmi_sleep_N2, cmi); % Append CMI results for N2
    elseif strcmp(condition, 'N3')  % Sleep stage N3
        gc_sleep_N3 = cat(3, gc_sleep_N3, F);  % Append GC results for N3
        cmi_sleep_N3 = cat(3, cmi_sleep_N3, cmi); % Append CMI results for N3
    elseif strcmp(condition, 'R')  % Sleep stage R (REM)
        gc_sleep_R = cat(3, gc_sleep_R, F);  % Append GC results for R
        cmi_sleep_R = cat(3, cmi_sleep_R, cmi); % Append CMI results for R
    end
end

% Step 6: Perform Mann-Whitney group testing (Wake vs each sleep stage) and FDR correction
% Step 6: Perform Mann-Whitney group testing (Wake vs each sleep stage) and FDR correction

% Initialize storage for z-scores and p-values for each comparison
z_gc = [];
p_gc = [];
z_cmi = [];
p_cmi = [];


% Perform Mann-Whitney Group test for Granger causality
[z_gc_N1, p_gc_N1] = run_mann_whitney_group(gc_wake, gc_sleep_N1);
[z_gc_N2, p_gc_N2] = run_mann_whitney_group(gc_wake, gc_sleep_N2);
[z_gc_N3, p_gc_N3] = run_mann_whitney_group(gc_wake, gc_sleep_N3);
[z_gc_R, p_gc_R] = run_mann_whitney_group(gc_wake, gc_sleep_R);

% Perform Mann-Whitney Group test for CMI
[z_cmi_N1, p_cmi_N1] = run_mann_whitney_group(cmi_wake, cmi_sleep_N1);
[z_cmi_N2, p_cmi_N2] = run_mann_whitney_group(cmi_wake, cmi_sleep_N2);
[z_cmi_N3, p_cmi_N3] = run_mann_whitney_group(cmi_wake, cmi_sleep_N3);
[z_cmi_R, p_cmi_R] = run_mann_whitney_group(cmi_wake, cmi_sleep_R);

% Perform FDR correction
alpha = 0.05;  % Significance level

[sig_gc_N1, pcrit_gc_N1] = significance(p_gc_N1, alpha, 'FDRD');
[sig_gc_N2, pcrit_gc_N2] = significance(p_gc_N2, alpha, 'FDRD');
[sig_gc_N3, pcrit_gc_N3] = significance(p_gc_N3, alpha, 'FDRD');
[sig_gc_R, pcrit_gc_R] = significance(p_gc_R, alpha, 'FDRD');

[sig_cmi_N1, pcrit_cmi_N1] = significance(p_cmi_N1, alpha, 'FDRD');
[sig_cmi_N2, pcrit_cmi_N2] = significance(p_cmi_N2, alpha, 'FDRD');
[sig_cmi_N3, pcrit_cmi_N3] = significance(p_cmi_N3, alpha, 'FDRD');
[sig_cmi_R, pcrit_cmi_R] = significance(p_cmi_R, alpha, 'FDRD');


% Step 7: Plot effect sizes (z-scores) and mark significant differences with white dots


% Plot Granger causality effect sizes and significance
plot_effect_size_diff(z_gc_N1, sig_gc_N1, 'Granger Causality (Wake vs N1)', 'gc_z_N1.png');
plot_effect_size_diff(z_gc_N2, sig_gc_N2, 'Granger Causality (Wake vs N2)', 'gc_z_N2.png');
plot_effect_size_diff(z_gc_N3, sig_gc_N3, 'Granger Causality (Wake vs N3)', 'gc_z_N3.png');
plot_effect_size_diff(z_gc_R, sig_gc_R, 'Granger Causality (Wake vs REM)', 'gc_z_R.png');

% Plot CMI effect sizes and significance
plot_effect_size_diff(z_cmi_N1, sig_cmi_N1, 'Conditional Mutual Information (Wake vs N1)', 'cmi_z_N1.png');
plot_effect_size_diff(z_cmi_N2, sig_cmi_N2, 'Conditional Mutual Information (Wake vs N2)', 'cmi_z_N2.png');
plot_effect_size_diff(z_cmi_N3, sig_cmi_N3, 'Conditional Mutual Information (Wake vs N3)', 'cmi_z_N3.png');
plot_effect_size_diff(z_cmi_R, sig_cmi_R, 'Conditional Mutual Information (Wake vs REM)', 'cmi_z_R.png');

% Step 8: Save Data
save(fullfile(results_dir, 'gc_cmi_results.mat'), 'z_gc_N1', 'z_gc_N2', 'z_gc_N3', 'z_gc_R', ...
    'z_cmi_N1', 'z_cmi_N2', 'z_cmi_N3', 'z_cmi_R', ...
    'sig_gc_N1', 'sig_gc_N2', 'sig_gc_N3', 'sig_gc_R', ...
    'sig_cmi_N1', 'sig_cmi_N2', 'sig_cmi_N3', 'sig_cmi_R');

fprintf('Cross-subject statistical analysis completed and results saved.\n');

%% HELPER FUNCTIONS

% Define a helper function to run Mann-Whitney Group test
function [z_vals, p_vals] = run_mann_whitney_group(gc_wake, gc_sleep)
    z_vals = [];
    p_vals = [];
    n_regions = size(gc_wake, 1);
    n_subjects = size(gc_wake, 3);  % Number of subjects
    
    % Loop through all pairs of regions and perform the Mann-Whitney U test
    for region = 1:n_regions
        wake_data = squeeze(gc_wake(region, :, :));  % Extract wake data for this region pair
        sleep_data = squeeze(gc_sleep(region, :, :)); % Extract sleep data for this region pair
        
        % Convert the wake and sleep data into cell arrays for Mann-Whitney Group
        wake_cell = cell(1, n_subjects);
        sleep_cell = cell(1, n_subjects);
        
        for subj = 1:n_subjects
            wake_cell{subj} = wake_data(:, subj);  % Each subject's data in a separate cell
            sleep_cell{subj} = sleep_data(:, subj); % Each subject's data in a separate cell
        end
        
        % Run Mann-Whitney U test
        [z, U] = mann_whitney_group(wake_cell, sleep_cell);
        
        % Calculate p-values from z-scores
        pval = erfc(abs(z)/sqrt(2)); % two-tailed significance test
        z_vals = [z_vals, z];
        p_vals = [p_vals, pval];
    end
end

function plot_effect_size_diff(z_vals, sig_vals, title_text, filename)
    figure;
    imagesc(z_vals);  % Plot z-scores (effect sizes)
    colormap(redbluecmap);  % Use red-blue colormap
    colorbar;
    title(title_text);
    set(gca, 'XTick', 1:size(z_vals, 2), 'YTick', 1:size(z_vals, 1));
    xlabel('Regions');
    ylabel('Regions');
    
    % Mark significant differences with white dots
    hold on;
    for i = 1:size(sig_vals, 1)
        for j = 1:size(sig_vals, 2)
            if sig_vals(i, j)
                plot(j, i, 'wo', 'MarkerSize', 5, 'LineWidth', 2);  % White dot for significance
            end
        end
    end
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(results_dir, filename));
end

