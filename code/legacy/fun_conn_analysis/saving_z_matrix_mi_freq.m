% MATLAB code to extract z-scores from the mi_time_effect_sizes structure
% for each frequency band, condition, and region, and save them in separate matrices.

% Step 1: Load the .mat file
file_path = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_freq_effect_sizes.mat';
data = load(file_path);

% Step 2: Initialize variables
conditions = {'N1', 'N2', 'N3', 'R_'};
frequencies = {'delta', 'theta', 'alpha', 'sigma', 'beta'};
num_regions = 68;  % Number of regions (Region_1 to Region_68)
num_conditions = length(conditions);
num_frequencies = length(frequencies);

% Initialize empty matrices to hold the z-scores for each frequency band
z_matrix_delta = zeros(num_regions, num_conditions);
z_matrix_theta = zeros(num_regions, num_conditions);
z_matrix_alpha = zeros(num_regions, num_conditions);
z_matrix_sigma = zeros(num_regions, num_conditions);
z_matrix_beta = zeros(num_regions, num_conditions);

% Step 3: Extract the z-scores for each condition, frequency band, and region, and store in matrices
for freq_idx = 1:num_frequencies
    freq_band = frequencies{freq_idx};  % 'delta', 'theta', 'alpha', 'beta', 'gamma'
    
    for cond_idx = 1:num_conditions
        condition = conditions{cond_idx};  % 'N1', 'N2', 'N3', 'R_'
        
        % Loop through all 68 regions for this condition and frequency band
        for region_idx = 1:num_regions
            region_name = sprintf('Region_%d', region_idx);  % Dynamically generate region name without leading zeros
            
            % Get the z score for the current region, condition, and frequency band
            z_score = data.mi_freq_effect_sizes.(condition).(freq_band).(region_name).z;
            
            % Store the z score in the appropriate matrix for the frequency band
            if strcmp(freq_band, 'delta')
                z_matrix_delta(region_idx, cond_idx) = z_score;
            elseif strcmp(freq_band, 'theta')
                z_matrix_theta(region_idx, cond_idx) = z_score;
            elseif strcmp(freq_band, 'alpha')
                z_matrix_alpha(region_idx, cond_idx) = z_score;
            elseif strcmp(freq_band, 'sigma')
                z_matrix_sigma(region_idx, cond_idx) = z_score;
            elseif strcmp(freq_band, 'beta')
                z_matrix_beta(region_idx, cond_idx) = z_score;
            end
        end
    end
end

% Step 4: Save the matrices in the same directory
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_freq_zscores_delta.mat', 'z_matrix_delta');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_freq_zscores_theta.mat', 'z_matrix_theta');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_freq_zscores_alpha.mat', 'z_matrix_alpha');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_freq_zscores_sigma.mat', 'z_matrix_sigma');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/mi/mi_freq_zscores_beta.mat', 'z_matrix_beta');
