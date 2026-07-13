% Script to run a Mann-Whitney U rank-sum test on the dynamical dependence
% values across anaesthetic conditions.

% comparison of anaesthetic conditions compared to wake, creating 3
% independent tests. These tests are independent because the same subject
% is not administered a wider variety of the anaesthetic, but only one.

% here we have, wake vs propofol, wake vs ketamine, wake vs xenon. we can
% either do this for the preoptimisation runs in which we have the same
% sized cell arrays, or we can do it for the optimisation condition in
% which case the cell arrays are off different size. 


clear all;
clc;

%% Load in the data for each individual

data_dir = '/Volumes/dataSets/restEEGHealthySubjects/restEEGHealthySubjects/AnesthesiaProjectEmergence/results/ssdiData/';

% define each anaesthetic condition
conditions = {'P', 'X', 'K'};

% initialise a structure that will hold all conditions

condition_data = struct();

% then initialise cell arrays that will hold the data for each condition

for cond = conditions
    condition_data.(cond{1}) = {};
end

% do the same for the wake data
wake_data = struct();

for cond = conditions
    wake_data.(cond{1}) = {};
end

%pcondition_data = {}; % initialise a cell array that will store the column vectors containing dd_values for each subject under a specific condition
%wcondition_data = {};

z_score_final = [];
U_final = [];
pval_final = [];
zcrit_final = [];
pcrit_final = [];
sigs_final = [];

%%
for n = 2:9

    data_files = dir(fullfile(data_dir, sprintf("*_mdim_%d_dynamical_dependence.mat", n)));


    for i = 2:length(data_files)
        file_name = data_files(i).name; % gets the name of file 'i'
        file_dir = fullfile(data_dir, file_name); % gets the full file directory

        for cond = conditions

            if contains(file_name, cond{1})
                fprintf(1, "Now reading %s \n", file_name); % prints a string informing you that it is reading the file into workspace
                dd_data = load(file_dir); % loads in data into a variable called 'dd_data'
                dd_data = dd_data.dopto;

            % Flip the data to make it a (n, 1) column vector
                dd_data_column = dd_data';

            % Concatenate this column vector to the right of the existing matrix
                condition_data.(cond{1}){end+1} = dd_data_column;

            % Look for corresponding wake condition data files
                name_prefix = file_name(1:5);
                w_file_index = find(contains({data_files.name}, [name_prefix, 'W']),1);
                if ~isempty(w_file_index)
                    w_file_name = data_files(w_file_index).name;
                    % now load in the wake condition file
                    fprintf(1, "Now reading %s \n", w_file_name);
                    w_data = load(fullfile(data_dir, w_file_name));
                    w_data = w_data.dopto;
                    w_data_column = w_data';

                    wake_data.(cond{1}){end+1} = w_data_column;
                end
            end
        end
    end


    %% Running Mann-Whitney U rank-sum test

    zscores_all = [];
    U_all = [];

    for cond = conditions
        [z_score, U] = mann_whitney_group(wake_data.(cond{1}), condition_data.(cond{1})); % calculate z score for each condition
        zscores_all = [zscores_all, z_score];
        pval = erfc(abs(zscores_all)/sqrt(2)); % two-tailed significance test
        U_all = [U_all, U];
    end

    z_score_final = [z_score_final; zscores_all]
    pval_final = [pval_final; pval]
    U_final = [U_final; U_all]



    %%

    % multiple comparison with FDR

    [sigs, pcrit] = significance(pval, 0.01, 'FDR');
    pcrit_final = [pcrit_final; pcrit]
    sigs_final = [sigs_final; sigs]

    % critical z-score
    zcrit = abs(sqrt(2)*erfcinv(pcrit));% two-tailed (+/-)
    zcrit_final = [zcrit_final; zcrit]
end


% %% Saving file for each n-macro
% 
% results_dir = '/Users/borjan/code/python/AnesthesiaProjectEmergence/results/stats/';
% results_file = 'nmacro_stats_opt';
% save_file = [fullfile(results_dir, results_file) '.mat'];
% fprintf('saving stats for ''%s''... \n', results_file);
% save(save_file, 'z_score_final', 'pval_final', 'U_final', 'zcrit_final', 'pcrit_final', 'sigs_final');

%% PLOT Z-score figures:

% Pastel colors definitions
pastelGreen = [0.6, 0.9, 0.6]; % Lighter green
pastelRed = [0.9, 0.6, 0.6]; % Lighter red

nFINDs = 2:9; % n-FINDs range from 2 to 9
columnTitles = {'Propofol', 'Xenon', 'Ketamine'}; % Titles for each column

figure; % Create a new figure for all conditions/measurements
for i = 1:size(z_score_final, 2) % Loop through each column
    subplot(1, 3, i); % Create a subplot for each condition/measurement
    for j = 1:length(nFINDs)
        % Determine color based on z-score value
        if z_score_final(j, i) >= 0
            barColor = pastelGreen; % Pastel green for positive
        else
            barColor = pastelRed; % Pastel red for negative
        end
        bar(nFINDs(j), z_score_final(j, i), 'FaceColor', barColor); % Plot each bar individually to color them
        hold on; % Keep the figure open to plot the next bar
    end
    
    % Set x-tick positions and labels
    xticks(nFINDs); % Set x-tick positions to match nFINDs
    xticklabels(arrayfun(@num2str, nFINDs, 'UniformOutput', false)); % Label each x-tick with n-FINDs
    
    title(columnTitles{i}, 'FontSize', 20); % Add title with larger font
    xlabel('n-FINDs', 'FontSize', 14); % x-axis label with larger font
    ylabel('Z-Score', 'FontSize', 14); % y-axis label with larger font
    set(gca, 'FontSize', 12); % Larger font for x-ticks and y-ticks
    xlim([1.5 9.5]); % Adjust x-axis limits to make the plot look nicer
end

% sgtitle('Z-Score Histograms for Each Condition', 'FontSize', 18); % Super title for the entire figure

% Save the figure
filename = 'Z-Score_Histograms_All_Conditions.png'; % Create a filename for the figure
saveas(gcf, filename); % Save the figure as a PNG file
% For a high-resolution version, you might want to use:
print(filename, '-dpng', '-r300'); % Save as high-resolution PNG

%% Plotting for publication

% Pastel colors definitions
pastelGreen = [0.6, 0.9, 0.6]; % Lighter green
pastelRed = [0.9, 0.6, 0.6]; % Lighter red

nFINDs = 2:9; % n-FINDs range from 2 to 9
columnTitles = {'Propofol', 'Xenon', 'Ketamine'}; % Titles for each column

for i = 1:size(z_score_final, 2) % Loop through each column
    figure; % Create a new figure for each condition/measurement
    for j = 1:length(nFINDs)
        % Determine color based on z-score value
        if z_score_final(j, i) >= 0
            barColor = pastelGreen; % Pastel green for positive
        else
            barColor = pastelRed; % Pastel red for negative
        end
        bar(nFINDs(j), z_score_final(j, i), 'FaceColor', barColor); % Plot each bar individually to color them
        hold on; % Keep the figure open to plot the next bar
    end
    hold off; % Release the figure

    % Set x-tick positions and labels
    xticks(nFINDs); % Set x-tick positions to match nFINDs
    xticklabels(arrayfun(@num2str, nFINDs, 'UniformOutput', false)); % Label each x-tick with n-FINDs

    %title(columnTitles{i}, 'FontSize', 20); % Add title with larger font
    xlabel('n-FINDs', 'FontSize', 14); % x-axis label with larger font
    ylabel('Z-Score', 'FontSize', 14); % y-axis label with larger font
    set(gca, 'FontSize', 12); % Larger font for x-ticks and y-ticks
    xlim([1.5 9.5]); % Adjust x-axis limits to make the plot look nicer

    % Save the figure
    filename = sprintf('%s_Z-Score_Histogram.png', columnTitles{i}); % Create a filename for each figure
    saveas(gcf, filename); % Save the figure as a PNG file
    % For a high-resolution version, you might want to use:
    print(filename, '-dpng', '-r300'); % Save as high-resolution PNG
end



