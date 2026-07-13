% MATLAB code to extract z-scores from the mw_in_time_effect_sizes structure
% and save them in a matrix where each column represents a condition
% and each row represents a region.

% Step 1: Load the .mat file
file_path = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_in_time_mw_effect_sizes.mat';
data = load(file_path);

% Step 2: Initialize variables
conditions = {'N1', 'N2', 'N3', 'R_'};
num_regions = 68;  % Number of regions
num_conditions = length(conditions);

% Initialize an empty matrix to hold the z-scores
z_matrix = zeros(num_regions, num_conditions);

% Step 3: Extract the z-scores for each condition and store in the matrix
for cond_idx = 1:num_conditions
    condition = conditions{cond_idx};
    
    % Get the z field for the current condition and convert it to a vector
    z_scores = data.mw_in_time_effect_sizes.(condition);
    z_vector = [z_scores.z];  % This will concatenate all z-scores into a single array
    
    % Store the z_vector in the appropriate column of z_matrix
    z_matrix(:, cond_idx) = z_vector;
end

% Step 4: Save the matrix in the same directory
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_in_time_zscores_matrix.mat', 'z_matrix');
