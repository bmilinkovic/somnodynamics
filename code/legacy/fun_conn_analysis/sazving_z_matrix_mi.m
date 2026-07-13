% MATLAB code to extract z-scores from the mi_time_effect_sizes structure
% and save them in a matrix where each column represents a condition
% and each row represents a region.

% Step 1: Load the .mat file
file_path = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_time_effect_sizes.mat';
data = load(file_path);

% Step 2: Initialize variables
conditions = {'N1', 'N2', 'N3', 'R_'};
num_regions = 68;  % Number of regions (Region_1 to Region_68)
num_conditions = length(conditions);

% Initialize an empty matrix to hold the z-scores
z_matrix = zeros(num_regions, num_conditions);

% Step 3: Extract the z-scores for each condition and store in the matrix
for cond_idx = 1:num_conditions
    condition = conditions{cond_idx};
    
    % Loop through all 68 regions for this condition
    for region_idx = 1:num_regions
        % Create the region field name without leading zero (e.g., 'Region_1', 'Region_2', ..., 'Region_68')
        region_name = sprintf('Region_%d', region_idx);  
        
        % Get the z score for the current region and condition
        z_matrix(region_idx, cond_idx) = data.mi_time_effect_sizes.(condition).(region_name).z;
    end
end

% Step 4: Save the matrix in the same directory
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_time_zscores_matrix.mat', 'z_matrix');
