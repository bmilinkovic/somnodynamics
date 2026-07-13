%% 1. Parameters
if ~exist('momax', 'var'), momax = 6; end % maximum model order
if ~exist('moregmode', 'var'), moregmode = 'LWR'; end % model order regression mode 
if ~exist('mosel', 'var'), mosel = 'AIC'; end % which model order to select
if ~exist('plotm', 'var'), plotm = 0; end % for plotting on separate figures
if ~exist('regmode', 'var'), regmode = 'LWR'; end
if ~exist('test_stat', 'var'), test_stat = 'F'; end
if ~exist('alpha', 'var'), alpha = 0.01; end
if ~exist('mhtc', 'var'), mhtc = 'FDRD'; end
if ~exist('fres', 'var'), fres = 256; end

%% 1. Set up paths and directories
% Define the file path to preprocessed data directory
src_data = '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data/';

% Define the file path to the results of the pair-wise condition Granger causality parameters
pwcgc_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc_lowres/';

% Create directory for results if it doesn't exist
if exist(pwcgc_dir, 'dir') == 0
    mkdir(pwcgc_dir)
end

%% 2. Load data and loop through participant folders
% List all participant folders in the src_data directory
participant_folders = dir(src_data);
participant_folders = participant_folders([participant_folders.isdir]); % Only directories

% Loop over each participant folder
for i = 1:length(participant_folders)
    % Skip '.' and '..' directories
    if strcmp(participant_folders(i).name, '.') || strcmp(participant_folders(i).name, '..')
        continue;
    end
    
    % Get current participant folder
    participant_folder = fullfile(src_data, participant_folders(i).name);
    
    % Check if there is a 'mat' subdirectory
    mat_folder = fullfile(participant_folder, 'mat');
    if ~exist(mat_folder, 'dir')
        warning('Mat folder does not exist in %s', participant_folder);
        continue;
    end
    
    % List all .mat files in the 'mat' subfolder
    files = dir(fullfile(mat_folder, '*.mat'));
    
    % Loop over each .mat file in the folder
    for fileNumber = 1:length(files)
        filename = files(fileNumber).name;
        folder = files(fileNumber).folder;

        tic;
        fprintf('..Starting on subject %s (%g / %g) \n', filename, fileNumber, length(files));
        
        % Load the .mat file
        loadedData = load(fullfile(folder, filename));
        
        % Extract 'source_ts' from the loaded structure
        if isfield(loadedData, 'source_ts')
            data = loadedData.source_ts;  % shape (n_regions, n_timesteps, n_epochs)
        else
            warning('No "source_ts" field in file %s', filename);
            continue;
        end
        
        % Perform the Dickey-Fuller test for each region and epoch
        [n_regions, n_timesteps, n_epochs] = size(data);
        fprintf('Dickey-Fuller test results for %s:\n', filename);
        for region = 1:n_regions
            for epoch = 1:n_epochs
                % Extract the time series for the current region and epoch
                ts = squeeze(data(region, :, epoch));
                
                % Perform the Augmented Dickey-Fuller test
                [h, pValue, stat, cValue] = adftest(ts);
                
                % Print the results
                if h == 1
                    fprintf('Region %d, Epoch %d: Stationary (p = %.4f)\n', region, epoch, pValue);
                else
                    fprintf('Region %d, Epoch %d: Non-stationary (p = %.4f)\n', region, epoch, pValue);
                end
            end
        end
        
        % Demean the data if needed (assumed function 'demean' exists)
        if exist('demean', 'file')
            data = demean(data, true);
        end
        
        % Model order selection and VAR model building
        [moaic, mobic] = tsdata_to_varmo(data, momax, moregmode, [], [], []);
        morder = moselect(sprintf('VAR model order selection, max = %d', momax), mosel, 'AIC', moaic, 'BIC', mobic);
        [A, V] = tsdata_to_var(data, morder, regmode);
        info = var_info(A, V);

        % Pairwise Granger Causality and significance testing
        F = var_to_pwcgc(A, V);
        tstat = var_to_pwcgc_tstat(data, V, morder, regmode, test_stat);

        % Get dimensions for significance testing
        nvars = size(A, 1);
        nobs = size(data, 2);
        ntrials = size(data, 3);

        % p-value and significance testing
        pval = mvgc_pval(tstat, test_stat, 1, 1, nvars-2, morder, nobs, ntrials);
        sig = significance(pval, alpha, mhtc);

        % Save eweights matrix
        eweights = F / nanmax(F(:)); 
        eweights_file = sprintf('%seweights_matrix_%s.mat', pwcgc_dir, filename(1:10));
        save(eweights_file, 'eweights');

        % Convert to SS parameters for DI analysis
        [Ass, Vss] = transform_var(A, V);
        [A, C, K] = var_to_ss(Ass);
        
        CAK = Ass;

        H = var2trfun(Ass, fres);

        % Save the PWCGC model
        modfile = sprintf('%spwcgc_params_%s.mat', pwcgc_dir, filename(1:10));
        fprintf('Saving PWCGC model ''%s''. ', modfile);
        save(modfile, 'V', 'CAK', 'H', 'Vss');
        fprintf('Saved..\n');
        toc;
    end
    clear V CAK H Vss A C K eweights sig pval F tstat 
end

