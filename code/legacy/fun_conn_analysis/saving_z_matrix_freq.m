% MATLAB code to extract z-scores from the mw_in_effect_sizes structure
% for each frequency band and condition, and save them in separate matrices.

% Step 1: Load the .mat file
file_path = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_out_freq_mw_effect_sizes.mat';
data = load(file_path);

% Step 2: Initialize variables
conditions = {'N1', 'N2', 'N3', 'R_'};
frequencies = {'delta', 'theta', 'alpha', 'sigma', 'beta'};
num_regions = 68;  % Number of regions
num_conditions = length(conditions);
num_frequencies = length(frequencies);

% Initialize empty matrices to hold the z-scores for each frequency band
z_matrix_delta = zeros(num_regions, num_conditions);
z_matrix_theta = zeros(num_regions, num_conditions);
z_matrix_alpha = zeros(num_regions, num_conditions);
z_matrix_sigma = zeros(num_regions, num_conditions);
z_matrix_beta = zeros(num_regions, num_conditions);

% Step 3: Extract the z-scores for each condition and frequency band, and store in matrices
for freq_idx = 1:num_frequencies
    freq_band = frequencies{freq_idx};  % 'delta', 'theta', 'alpha', 'beta', 'gamma'
    
    for cond_idx = 1:num_conditions
        condition = conditions{cond_idx};  % 'N1', 'N2', 'N3', 'R_'
        
        % Extract the z-scores for the current condition and frequency band
        z_scores = data.mw_out_effect_sizes.(condition).(freq_band);
        z_vector = [z_scores.z];  % This will concatenate all z-scores into a single array
        
        % Store the z_vector in the appropriate matrix for the frequency band
        if strcmp(freq_band, 'delta')
            z_matrix_delta(:, cond_idx) = z_vector;
        elseif strcmp(freq_band, 'theta')
            z_matrix_theta(:, cond_idx) = z_vector;
        elseif strcmp(freq_band, 'alpha')
            z_matrix_alpha(:, cond_idx) = z_vector;
        elseif strcmp(freq_band, 'sigma')
            z_matrix_sigma(:, cond_idx) = z_vector;
        elseif strcmp(freq_band, 'beta')
            z_matrix_beta(:, cond_idx) = z_vector;
        end
    end
end

% Step 4: Save the matrices in the same directory
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_out_freq_zscores_delta.mat', 'z_matrix_delta');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_out_freq_zscores_theta.mat', 'z_matrix_theta');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_out_freq_zscores_alpha.mat', 'z_matrix_alpha');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_out_freq_zscores_sigma.mat', 'z_matrix_sigma');
save('/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/gc_out_freq_zscores_beta.mat', 'z_matrix_beta');
