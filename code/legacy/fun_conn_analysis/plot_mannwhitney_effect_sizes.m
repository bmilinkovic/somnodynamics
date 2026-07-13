% plot_mannwhitney_effect_sizes.m
%
% This script performs the following tasks:
%
% 1. Loads the Desikan-Killiany atlas region labels from the source time series file of one participant.
%    These labels correspond to the indices used in the Granger causality (GC) analysis, where each
%    index represents a brain region.
%
% 2. Loads the cross-subject effect sizes (from the Mann-Whitney U test) for each condition and frequency band.
%    These effect sizes are computed as rank correlations, normalized to lie between -1 and 1, where:
%    - Blue indicates a decrease in GC (compared to wake).
%    - Red indicates an increase in GC (compared to wake).
%
% 3. Plots the effect sizes on a brain surface using the Desikan-Killiany atlas with 68 regions.
%
% Inputs:
%   - Directory to the source time series file of one participant to extract region labels.
%   - Mann-Whitney U test results (effect sizes) previously saved in 'mw_effect_sizes.mat'.
%
% Outputs:
%   - Brain surface plot visualizing the effect sizes for different conditions and frequency bands.
%
% Example Usage:
%   Run the script to load effect sizes and plot them on the brain surface for a specified condition (e.g., N1)
%   and frequency band (e.g., delta). The script uses the Desikan-Killiany atlas to map the effect sizes
%   to the brain regions.
%
% Author: [Your Name]
% Date: [Date]

% Step 1: Load the Desikan-Killiany labels from the source time series file
label_file = '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data/EPCTL01/mat/EPCTL01_N1_source_time_series.mat';
data = load(label_file, 'labels');  % Load the labels variable
labels = data.labels;  % Extract the region labels (Desikan-Killiany atlas)

% Step 2: Load the Mann-Whitney U test effect sizes
mw_effect_sizes_file = 'mw_effect_sizes.mat';  % The file where the cross-subject test results are saved
load(mw_effect_sizes_file, 'mw_effect_sizes');  % Load the Mann-Whitney U test results

% Step 3: Set parameters for the plot
% Choose condition and frequency band to plot
condition = 'N1';   % Example condition (N1)
band = 'delta';     % Example frequency band (delta)

% Extract the effect sizes for the selected condition and frequency band
% mw_effect_sizes.(condition).(band) contains the rank correlation results for each brain region
effect_sizes = mw_effect_sizes.(condition).(band);

% Step 4: Prepare the brain surface plot
% Load a brain surface mesh (this should be a brain model with 68 regions corresponding to the Desikan-Killiany atlas)
% You can use a toolbox like BrainNet Viewer, SPM, or a MATLAB surface plotting method
brain_surface = load_brain_surface_mesh();  % Function to load the brain surface mesh (replace with your own method)

% Assuming 'labels' gives the region names in order, and 'effect_sizes' corresponds to those regions
% For example, labels{1} corresponds to the first region in the atlas, effect_sizes(1) corresponds to that region.

% Plot the brain surface and map effect sizes
figure;
trisurf(brain_surface.faces, brain_surface.vertices(:, 1), brain_surface.vertices(:, 2), brain_surface.vertices(:, 3), ...
    'FaceVertexCData', effect_sizes, 'EdgeColor', 'none', 'FaceColor', 'interp');
colormap('redbluecmap');  % Use a red-blue colormap for visualizing effect sizes
colorbar;  % Add a colorbar to indicate the scale of effect sizes
caxis([-1, 1]);  % Set the color axis limits to match the range of effect sizes (-1 to 1)
title(sprintf('Effect Sizes (Condition: %s, Band: %s)', condition, band));

% Step 5: Adjust plot settings
axis equal;
lighting phong;
camlight headlight;
view([90, 0]);  % Set the view to display the brain from the side

