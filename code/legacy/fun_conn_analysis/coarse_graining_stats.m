% GroupWise_GWC_GC_Comparison.m
% -------------------------------------------------------------------------
% This script performs a group-wise comparison of Granger-Geweke Causality
% (GWC-GC) F matrices between Wake and Sleep conditions across multiple
% subjects. It defines brain region groups, computes statistical test
% statistics, calculates normalized effect sizes, applies multiple
% comparison corrections, and visualizes the results.
% -------------------------------------------------------------------------

%% ---------------------------- Setup -------------------------------------

% Clear workspace, command window, and close all figures
clear; clc; close all;

% Add necessary paths (modify these paths based on your directory structure)
% Example:
% addpath('/path/to/your/functions');
% addpath('/path/to/your/data');

% Define the path to the data file
dataFile = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/allGWC_GC_results.mat'; % Modify as needed

% Define the directory where results will be saved
resultsDir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/coarse-grained-stats/'; % Modify as needed
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Load the data
if exist(dataFile, 'file')
    load(dataFile, 'allResults');
    fprintf('Data loaded successfully from %s.\n', dataFile);
else
    error('Data file not found: %s', dataFile);
end

% Define the list of subjects
subjectNumbers = sprintfc('%02d', 1:29); % '01', '02', ..., '29'
subjects = strcat('EPCTL', subjectNumbers);

% Define conditions
wakeCondition = 'W';
sleepConditions = {'N1', 'N2', 'N3', 'R'}; % Add more sleep conditions if applicable

% Define group names (ensure this matches your group definitions)
groupNames = {'Frontal Lobe', 'Parietal Lobe', 'Temporal Lobe', 'Occipital Lobe', 'Limbic System'};

% -------------------------------------------------------------------------

%% ---------------- Step 1: Define Groups Based on Labels from a Sample File

fprintf('\n--- Step 1: Define Groups Based on Labels from a Sample File ---\n');

% Define the path to a sample MAT file to extract labels
% Modify 'sampleSubject' and 'sampleCondition' based on your data structure
sampleSubject = subjects{1}; % e.g., 'EPCTL01'
sampleCondition = wakeCondition; % e.g., 'Wake'

% Define base directory where individual subject data is stored
baseDir = '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data'; % Modify as needed

% Construct the path to the sample MAT file
% Assuming the structure: baseDir/subject/mat/subject_condition_source_time_series.mat
sampleMatFile = fullfile(baseDir, sampleSubject, 'mat', sprintf('%s_%s_source_time_series.mat', sampleSubject, sampleCondition));

% Check if the sample MAT file exists and load labels
if exist(sampleMatFile, 'file')
    data = load(sampleMatFile, 'source_ts', 'labels');
    labels = data.labels; % Assuming labels are stored as a character array
    fprintf('Labels loaded from %s.\n', sampleMatFile);
else
    error('Sample MAT file not found: %s', sampleMatFile);
end

% Convert the labels from character matrix to a cell array of strings
labels = cellstr(labels);  % Converts the character matrix into a cell array of strings

% Define group names for reference and plotting (already defined earlier)

% Define groups using exact label matches with hemisphere suffixes
groups = { ...
    % Frontal Lobe
    find(ismember(labels, {'superiorfrontal-lh', 'superiorfrontal-rh', ...
                           'rostralmiddlefrontal-lh', 'rostralmiddlefrontal-rh', ...
                           'caudalmiddlefrontal-lh', 'caudalmiddlefrontal-rh', ...
                           'parsopercularis-lh', 'parsopercularis-rh', ...
                           'parsorbitalis-lh', 'parsorbitalis-rh', ...
                           'parstriangularis-lh', 'parstriangularis-rh', ...
                           'lateralorbitofrontal-lh', 'lateralorbitofrontal-rh', ...
                           'medialorbitofrontal-lh', 'medialorbitofrontal-rh', ...
                           'frontalpole-lh', 'frontalpole-rh'})), ... % Frontal Lobe

    % Parietal Lobe
    find(ismember(labels, {'superiorparietal-lh', 'superiorparietal-rh', ...
                           'inferiorparietal-lh', 'inferiorparietal-rh', ...
                           'supramarginal-lh', 'supramarginal-rh', ...
                           'precuneus-lh', 'precuneus-rh', ...
                           'postcentral-lh', 'postcentral-rh', ...
                           'paracentral-lh', 'paracentral-rh', ...
                           'precentral-lh', 'precentral-rh'})), ... % Parietal Lobe

    % Temporal Lobe
    find(ismember(labels, {'superiortemporal-lh', 'superiortemporal-rh', ...
                           'middletemporal-lh', 'middletemporal-rh', ...
                           'inferiortemporal-lh', 'inferiortemporal-rh', ...
                           'fusiform-lh', 'fusiform-rh', ...
                           'temporalpole-lh', 'temporalpole-rh', ...
                           'transversetemporal-lh', 'transversetemporal-rh', ...
                           'bankssts-lh', 'bankssts-rh'})), ... % Temporal Lobe

    % Occipital Lobe
    find(ismember(labels, {'lateraloccipital-lh', 'lateraloccipital-rh', ...
                           'cuneus-lh', 'cuneus-rh', ...
                           'lingual-lh', 'lingual-rh', ...
                           'pericalcarine-lh', 'pericalcarine-rh'})), ... % Occipital Lobe

    % Limbic System
    find(ismember(labels, {'insula-lh', 'insula-rh', ...
                           'entorhinal-lh', 'entorhinal-rh', ...
                           'parahippocampal-lh', 'parahippocampal-rh', ...
                           'posteriorcingulate-lh', 'posteriorcingulate-rh', ...
                           'isthmuscingulate-lh', 'isthmuscingulate-rh', ...
                           'caudalanteriorcingulate-lh', 'caudalanteriorcingulate-rh', ...
                           'rostralanteriorcingulate-lh', 'rostralanteriorcingulate-rh'})) ... % Limbic System
};

% Optionally, handle unassigned channels
% Initialize a logical array to keep track of assigned channels
assigned = false(length(labels), 1);

% Mark assigned channels
for gIdx = 1:length(groups)
    assigned(groups{gIdx}) = true;
end

% Find unassigned channels
unassignedIdx = find(~assigned);

% If there are unassigned channels, create an additional group
if ~isempty(unassignedIdx)
    groups{end+1} = unassignedIdx;
    groupNames{end+1} = 'Unassigned';
    fprintf('Unassigned Group created with %d channels.\n', length(unassignedIdx));
end

% Verify group assignments to ensure no overlaps
try
    [gnum, gsiz, idx] = check_group(groups, length(labels));
    fprintf('All groups have been defined correctly with no overlapping indices.\n');
catch ME
    warning('Group validation failed: %s', ME.message);
end

% Verify group definitions
for gIdx = 1:length(groups)
    fprintf('Group %d (%s) has %d channels.\n', gIdx, groupNames{gIdx}, length(groups{gIdx}));
end

% Define 'g' as the number of groups
g = length(groups);
fprintf('Number of groups (g): %d\n', g);

%% ---------------- Step 2: Organize F Matrices and Time Series Data for Comparison

fprintf('\n--- Step 2: Organize F Matrices and Time Series Data ---\n');

% Initialize data containers
% Dimensions: g x g x nSubjects
F_Wake = NaN(g, g, length(subjects));
F_Sleep = struct(); % To store each sleep condition separately

for sc = 1:length(sleepConditions)
    F_Sleep.(sleepConditions{sc}) = NaN(g, g, length(subjects));
end

% Populate the F_Wake and F_Sleep structures
for subjIdx = 1:length(subjects)
    subject = subjects{subjIdx};
    fprintf('Processing Subject: %s\n', subject);
    
    % Extract Wake F matrix
    if isfield(allResults, subject) && isfield(allResults.(subject), wakeCondition)
        F_Wake(:, :, subjIdx) = allResults.(subject).(wakeCondition).F;
    else
        warning('Wake condition F matrix missing for %s.', subject);
    end
    
    % Extract Sleep condition F matrices
    for sc = 1:length(sleepConditions)
        sleepCond = sleepConditions{sc};
        if isfield(allResults, subject) && isfield(allResults.(subject), sleepCond)
            F_Sleep.(sleepCond)(:, :, subjIdx) = allResults.(subject).(sleepCond).F;
        else
            warning('%s condition F matrix missing for %s.', sleepCond, subject);
        end
    end
end

%% ---------------- Step 3: Compute Test Statistics Using var_to_gwcgc_tstat

fprintf('\n--- Step 3: Compute Test Statistics Using var_to_gwcgc_tstat ---\n');

% Define VAR parameters
p = 1; % Adjust based on your VAR model
regmode = 'OLS'; % Adjust based on your regression mode
tstat = 'F'; % Choose 'F' for F-test or 'LR' for Likelihood Ratio test

% Initialize a structure to hold test statistics
testStats = struct();

for sc = 1:length(sleepConditions)
    sleepCond = sleepConditions{sc};
    fprintf('Computing Test Statistics for Sleep Condition: %s\n', sleepCond);
    
    % Initialize a matrix to hold statistics for all subjects
    % Dimensions: g x g x nSubjects
    stat_sleep = NaN(g, g, length(subjects));
    
    for subjIdx = 1:length(subjects)
        subject = subjects{subjIdx};
        fprintf('  Subject: %s\n', subject);
        
        % Check if both Wake and Sleep F matrices exist
        if isfield(allResults, subject) && isfield(allResults.(subject), wakeCondition) && isfield(allResults.(subject), sleepCond)
            % Extract time series data for Sleep condition
            if isfield(allResults.(subject).(sleepCond), 'X') && isfield(allResults.(subject).(wakeCondition), 'V')
                X_sleep = allResults.(subject).(sleepCond).X; % n_variables x n_timesteps x n_epochs
                V_wake = allResults.(subject).(wakeCondition).V; % n x n
                
                % Compute test statistics
                try
                    stat = var_to_gwcgc_tstat(X_sleep, V_wake, groups, p, regmode, tstat);
                    stat_sleep(:, :, subjIdx) = stat;
                catch ME
                    warning('Error computing test statistics for %s - %s: %s', subject, sleepCond, ME.message);
                end
            else
                warning('Missing X or V data for subject %s in condition %s.', subject, sleepCond);
            end
        else
            warning('Missing F matrices or data for subject %s in condition %s.', subject, sleepCond);
        end
    end
    
    % Store the test statistics for the current Sleep condition
    testStats.(sleepCond) = stat_sleep;
end

%% ---------------- Step 4: Calculate Effect Sizes Based on Test Statistics

fprintf('\n--- Step 4: Calculate Effect Sizes Based on Test Statistics ---\n');

% Initialize effect size containers
effectSizes = struct();

for sc = 1:length(sleepConditions)
    sleepCond = sleepConditions{sc};
    fprintf('Calculating Effect Sizes for Sleep Condition: %s\n', sleepCond);
    
    % Extract test statistics for the current Sleep condition
    stat_sleep = testStats.(sleepCond); % g x g x nSubjects
    
    % Compute mean across subjects for each group pair
    mean_stat = mean(stat_sleep, 3, 'omitnan');
    
    % Normalize effect sizes between -1 and +1
    max_abs_stat = max(abs(mean_stat(:)));
    if max_abs_stat > 0
        mean_stat_normalized = mean_stat / max_abs_stat;
    else
        mean_stat_normalized = mean_stat; % All zeros
    end
    
    % Store the normalized effect sizes
    effectSizes.(sleepCond) = mean_stat_normalized;
end

%% ---------------- Step 5: Apply Multiple Comparison Correction (FDR)

fprintf('\n--- Step 5: Apply Multiple Comparison Correction (FDR) ---\n');

% Initialize a structure to hold significance masks
significanceMasks = struct();

for sc = 1:length(sleepConditions)
    sleepCond = sleepConditions{sc};
    fprintf('Applying FDR Correction for Sleep Condition: %s\n', sleepCond);
    
    % Extract test statistics for this condition
    stat_sleep = testStats.(sleepCond); % g x g x nSubjects
    
    % Assuming 'F' test statistics, compute p-values
    % Degrees of freedom for F-test depend on your VAR model
    % Define df1 and df2 based on your model
    % Example placeholders (replace with actual degrees of freedom)
    df1 = p; % degrees of freedom numerator
    % df2 is typically (n_timesteps - p * g), but adjust based on your data
    % For illustration, let's assume n_timesteps = 100
    n_timesteps = size(allResults.(subjects{1}).(wakeCondition).X, 2);
    df2 = n_timesteps - p * g; % adjust as necessary
    
    if df2 <= 0
        error('Degrees of freedom df2 must be positive. Check your VAR model parameters.');
    end
    
    % Compute p-values from F-test statistics
    % For each group pair, average F-statistic across subjects and compute p-value
    mean_stat = mean(stat_sleep, 3, 'omitnan');
    p_values = 1 - fcdf(abs(mean_stat), df1, df2); % Two-tailed test if applicable
    p_values(1:g+1:end) = NaN; % Remove self-causality p-values by setting them to NaN
    
    % Flatten the p-values and remove NaNs for FDR
    p_flat = p_values(~isnan(p_values));
    
    % Apply FDR correction using the Benjamini-Hochberg method
    % Ensure you have the fdr_bh function available. Download from MATLAB File Exchange if needed.
    [h_fdr, ~, ~] = fdr_bh(p_flat, 0.05, 'pdep'); % 5% FDR
    
    % Initialize a significance mask matrix
    significance = false(g, g);
    
    % Assign significance values
    significance(~eye(g)) = h_fdr;
    
    % Store the significance mask
    significanceMasks.(sleepCond) = significance;
end

%% ---------------- Step 6: Visualize Effect Sizes with Significance

fprintf('\n--- Step 6: Visualize Effect Sizes with Significance ---\n');

% Check if brewermap is available
if exist('brewermap', 'file') == 2
    cmap = brewermap(256, 'RdBu');
else
    cmap = jet(256);
    warning('Brewermap not found. Using default jet colormap.');
end

% Loop through each Sleep condition for visualization
for sc = 1:length(sleepConditions)
    sleepCond = sleepConditions{sc};
    fprintf('Visualizing Effect Sizes for Sleep Condition: %s\n', sleepCond);
    
    % Extract normalized effect sizes and significance mask
    d_norm = effectSizes.(sleepCond); % g x g matrix
    significance = significanceMasks.(sleepCond); % g x g logical matrix
    
    % Create a figure
    figure('Name', sprintf('Effect Sizes: %s vs Wake', sleepCond), 'NumberTitle', 'off');
    
    % Plot heatmap
    imagesc(d_norm);
    colorbar;
    colormap(cmap);
    caxis([-1 1]); % Ensure color scaling between -1 and +1
    title(sprintf('Effect Sizes (Normalized Test Stat) - %s vs Wake', sleepCond));
    xlabel('Source Groups');
    ylabel('Target Groups');
    
    % Set axis ticks and labels
    set(gca, 'XTick', 1:g, 'XTickLabel', groupNames, ...
             'YTick', 1:g, 'YTickLabel', groupNames, ...
             'TickLabelInterpreter', 'none');
    
    % Improve visualization aesthetics
    axis square;
    set(gca, 'FontSize', 12);
    
    % Overlay significance markers
    hold on;
    [rows, cols] = find(significance);
    for k = 1:length(rows)
        plot(cols(k), rows(k), 'k*', 'MarkerSize', 8); % Black stars for significance
    end
    hold off;
    
    % Add legend or annotations as needed
    annotation('textbox', [0.15, 0.85, 0.1, 0.1], 'String', '* Significant', ...
               'EdgeColor', 'none', 'FontSize', 12, 'BackgroundColor', 'none');
    
    % Save the figure
    saveas(gcf, fullfile(resultsDir, sprintf('EffectSizes_Significant_%s_vs_Wake.png', sleepCond)));
    
    % Close the figure to free memory
    close(gcf);
end

%% ----------------------- Step 7: Final Reporting ----------------------------

fprintf('\n--- Step 7: Final Reporting and Interpretation ---\n');

% (This step involves manual interpretation based on the generated visualizations.)
% Review the saved heatmaps to identify patterns of increased or decreased connectivity.
% Significant changes are marked with black stars.

fprintf('Group-Wise Comparison of GWC-GC F Matrices Completed and Visualizations Saved in %s.\n', resultsDir);

%% ---------------------------- End of Script --------------------------------
