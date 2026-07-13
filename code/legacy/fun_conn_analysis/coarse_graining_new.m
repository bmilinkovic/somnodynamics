%% MATLAB Script: Perform Group-Wise Conditional Granger-Geweke Causality Analysis

% -------------------------------------------------------------------------
% Prerequisites:
% - Ensure that all auxiliary functions (`check_group`, `vardare`, 
%   `sserror`, `logdet`, `findin`, `var_to_gwcgc1`) are implemented and accessible in the MATLAB path.
% - Verify that the directory structure and file naming conventions are as described.
% -------------------------------------------------------------------------clea
clear all;

% Define the base directory for preprocessed data
baseDir = '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data';

% Define the list of subjects
subjectNumbers = sprintfc('%02d', 1:29); % '01', '02', ..., '29'

% Define the list of conditions as they appear in source time series files
sourceConditions = {'N1', 'N2', 'N3', 'R', 'W'};

% Define the corresponding list of conditions as they appear in pwcgc_params files
% For 'R' and 'W', append an underscore '_'; others remain the same
pwcgcConditions = {'N1', 'N2', 'N3', 'R_', 'W_'};

% Define group names (for reference and plotting)
groupNames = {'Frontal Lobe', 'Parietal Lobe', 'Temporal Lobe', 'Occipital Lobe', 'Limbic System'};

% Define the directory where GWC-GC results will be saved
resultsDir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/';
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Initialize a structure to hold all results
allResults = struct();

%%
% -------------------------------------------------------------------------
% Step 1: Define Groups Based on Labels from a Sample File
% -------------------------------------------------------------------------

% Temporary variables to hold labels
labels = [];

% Loop through one subject and one condition to extract labels
% Assuming labels are consistent across all subjects and conditions
sampleSubject = 'EPCTL01';
sampleCondition = 'N1';
sampleMatFile = fullfile(baseDir, sampleSubject, 'mat', sprintf('%s_%s_source_time_series.mat', sampleSubject, sampleCondition));

if exist(sampleMatFile, 'file')
    data = load(sampleMatFile, 'source_ts', 'labels');
    labels = data.labels; % 68 x 27 char array
else
    error('Sample MAT file not found: %s', sampleMatFile);
end

% Convert the labels from character matrix to a cell array of strings
labels = cellstr(labels);  % Converts the character matrix into a cell array of strings

% Define group names for reference and plotting
groupNames = {'Frontal Lobe', 'Parietal Lobe', 'Temporal Lobe', 'Occipital Lobe', 'Limbic System'};

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

%%
% -------------------------------------------------------------------------
% Step 2: Iterate Through Subjects and Conditions to Perform GWC-GC Analysis
% -------------------------------------------------------------------------

% Iterate through each subject
for subjIdx = 1:length(subjectNumbers)
    subject = ['EPCTL' subjectNumbers{subjIdx}];
    fprintf('Processing Subject: %s\n', subject);
    
    % Iterate through each condition
    for condIdx = 1:length(sourceConditions)
        sourceCondition = sourceConditions{condIdx};
        pwcgcCondition = pwcgcConditions{condIdx};
        fprintf('  Condition: %s\n', sourceCondition);
        
        % Define the path to the source time series file
        matFile = fullfile(baseDir, subject, 'mat', sprintf('%s_%s_source_time_series.mat', subject, sourceCondition));
        
        % Check if the MAT file exists
        if ~exist(matFile, 'file')
            warning('MAT file not found: %s. Skipping.', matFile);
            continue;
        end
        
        % Load the source time series data and labels
        data = load(matFile, 'source_ts', 'labels');
        source_ts = data.source_ts; % (n_variables, n_timesteps, n_epochs)
        current_labels = data.labels; % (n_variables, label_length)
        
        % Load VAR parameters (CAK and V) from the pwcgc_params directory
        paramsFile = fullfile(resultsDir, sprintf('pwcgc_params_%s_%s.mat', subject, pwcgcCondition));
        if ~exist(paramsFile, 'file')
            warning('Parameters file not found: %s. Skipping.', paramsFile);
            continue;
        end
        params = load(paramsFile, 'CAK', 'V');
        
        % Validate that VAR parameters are loaded
        if ~isfield(params, 'CAK') || ~isfield(params, 'V')
            warning('VAR parameters CAK and/or V not found in %s. Skipping.', paramsFile);
            continue;
        end
        
        CAK = params.CAK; % (n x n x p)
        V = params.V;      % (n x n)
        
        % Perform GWC-GC analysis using var_to_gwcgc1
        try
            F = var_to_gwcgc1(CAK, V, groups); % g x g matrix
        catch ME
            warning('Error in var_to_gwcgc1 for %s - %s: %s', subject, sourceCondition, ME.message);
            continue;
        end
        
        % Store the result
        allResults.(subject).(sourceCondition).F = F; % g x g matrix
        
        % Save the GWC-GC matrix for this subject-condition pair
        saveFile = fullfile(resultsDir, sprintf('GWC_GC_%s_%s.mat', subject, sourceCondition));
        save(saveFile, 'F');
    end
end

%%
% -------------------------------------------------------------------------
% Step 3: Save All Results
% -------------------------------------------------------------------------

% Save all results into a single MAT file for easy access
save(fullfile(resultsDir, 'allGWC_GC_results.mat'), 'allResults');

%%
% -------------------------------------------------------------------------
% Step 4: Create Summary Statistics Across Subjects for Each Condition
% -------------------------------------------------------------------------

% Initialize a structure to hold summary results
summaryResults = struct();

% Iterate through each condition to compute mean and standard deviation across subjects
for condIdx = 1:length(sourceConditions)
    sourceCondition = sourceConditions{condIdx};
    pwcgcCondition = pwcgcConditions{condIdx};
    fprintf('Computing summary statistics for Condition: %s\n', sourceCondition);
    
    % Initialize a temporary array to hold F matrices from all subjects for this condition
    F_allSubjects = [];
    
    for subjIdx = 1:length(subjectNumbers)
        subject = ['EPCTL' subjectNumbers{subjIdx}];
        if isfield(allResults, subject) && isfield(allResults.(subject), sourceCondition)
            F = allResults.(subject).(sourceCondition).F;
            F_allSubjects(:,:,end+1) = F; %#ok<AGROW>
        end
    end
    
    % Check if any data exists for this condition
    if ~isempty(F_allSubjects)
        % Compute mean and standard deviation across subjects, ignoring NaNs
        summaryResults.(sourceCondition).meanF = mean(F_allSubjects, 3, 'omitnan');
        summaryResults.(sourceCondition).stdF = std(F_allSubjects, 0, 3, 'omitnan');
    else
        warning('No data found for condition %s.', sourceCondition);
    end
end

% Save summary results
save(fullfile(resultsDir, 'summaryGWC_GC_results.mat'), 'summaryResults');

%%
% -------------------------------------------------------------------------
% Step 5: Visualization - Plotting Mean GWC-GC Matrices and F-Statistics with SD
% -------------------------------------------------------------------------

% Iterate through each condition to plot the mean GWC-GC matrix and F-statistics
for condIdx = 1:length(sourceConditions)
    sourceCondition = sourceConditions{condIdx};
    pwcgcCondition = pwcgcConditions{condIdx};
    
    if isfield(summaryResults, sourceCondition)
        meanF = summaryResults.(sourceCondition).meanF;
        stdF = summaryResults.(sourceCondition).stdF;
        
        % Plot Heatmap of Mean GWC-GC Matrix
        figure('Name', sprintf('Average GWC-GC Matrix - Condition: %s', sourceCondition));
        imagesc(meanF);
        colorbar;
        colormap('jet');
        title(sprintf('Average GWC-GC Matrix - Condition: %s', sourceCondition));
        xlabel('Source Groups');
        ylabel('Target Groups');
        set(gca, 'XTick', 1:length(groups), 'XTickLabel', groupNames, ...
                 'YTick', 1:length(groups), 'YTickLabel', groupNames, 'TickLabelInterpreter', 'none');
        axis square;
        
        % Save the heatmap figure
        saveas(gcf, fullfile(resultsDir, sprintf('MeanGWC_GC_%s.png', sourceCondition)));
        
        % Close the figure to free memory
        close(gcf);
        
        % Plot Mean GWC-GC with Standard Deviation (Separate Heatmaps)
        figure('Name', sprintf('Mean and Std Dev GWC-GC - Condition: %s', sourceCondition));
        subplot(1,2,1);
        imagesc(meanF);
        colorbar;
        colormap('jet');
        title(sprintf('Mean GWC-GC - %s', sourceCondition));
        xlabel('Source Groups');
        ylabel('Target Groups');
        set(gca, 'XTick', 1:length(groups), 'XTickLabel', groupNames, ...
                 'YTick', 1:length(groups), 'YTickLabel', groupNames, 'TickLabelInterpreter', 'none');
        axis square;
        
        subplot(1,2,2);
        imagesc(stdF);
        colorbar;
        colormap('hot');
        title(sprintf('Std Dev GWC-GC - %s', sourceCondition));
        xlabel('Source Groups');
        ylabel('Target Groups');
        set(gca, 'XTick', 1:length(groups), 'XTickLabel', groupNames, ...
                 'YTick', 1:length(groups), 'YTickLabel', groupNames, 'TickLabelInterpreter', 'none');
        axis square;
        
        % Save the combined heatmap figure
        saveas(gcf, fullfile(resultsDir, sprintf('MeanStdGWC_GC_%s.png', sourceCondition)));
        
        % Close the figure to free memory
        close(gcf);
        
        % Optional: Plot a network graph representing causality connections
        % Requires defining significant connections based on statistical tests
    else
        fprintf('No summary results found for condition %s.\n', sourceCondition);
    end
end

fprintf('GWC-GC Analysis Completed and Results Saved.\n');
