function run_mvgc_atlas(varargin)
%RUN_MVGC_ATLAS Atlas-aware MVGC runner (legacy + yeo7 side-by-side).
%
% This function does not modify legacy outputs. For Yeo7 it writes to:
%   code/results/yeo7/gc/pwcgc/

if ~exist('tsdata_to_varmo','file')
    error('MVGC toolbox functions not found on path. Run startup.m first.');
end

p = inputParser;
addParameter(p, 'atlas', 'yeo7');
addParameter(p, 'src_data', '');
addParameter(p, 'pwcgc_dir', '');
addParameter(p, 'momax', 6);
addParameter(p, 'moregmode', 'LWR');
addParameter(p, 'mosel', 'AIC');
addParameter(p, 'regmode', 'LWR');
addParameter(p, 'test_stat', 'F');
addParameter(p, 'alpha', 0.01);
addParameter(p, 'mhtc', 'FDRD');
addParameter(p, 'fres', 256);
addParameter(p, 'run_adf', true);
parse(p, varargin{:});
cfg = p.Results;

atlas = lower(string(cfg.atlas));
if strlength(cfg.src_data) == 0
    if atlas == "yeo7"
        cfg.src_data = '/Users/borjan/code/healthy-sleep-project/code/results/source_recon/yeo7/';
    else
        cfg.src_data = '/Users/borjan/code/healthy-sleep-project/code/results/preprocessed_data/';
    end
end

if strlength(cfg.pwcgc_dir) == 0
    if atlas == "yeo7"
        cfg.pwcgc_dir = '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/pwcgc/';
    else
        cfg.pwcgc_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc_lowres/';
    end
end

if exist(cfg.pwcgc_dir, 'dir') == 0
    mkdir(cfg.pwcgc_dir);
end

source_files = dir(fullfile(cfg.src_data, 'EPCTL*', 'mat', '*_source_time_series.mat'));
if isempty(source_files)
    error('No source files found in %s', cfg.src_data);
end

fprintf('Running MVGC for atlas="%s" on %d files.\n', atlas, numel(source_files));

for file_idx = 1:numel(source_files)
    filename = source_files(file_idx).name;
    folder = source_files(file_idx).folder;
    fprintf('.. (%d/%d) %s\n', file_idx, numel(source_files), filename);

    parsed = parse_source_filename(filename);
    if parsed.valid == false
        warning('Skipping file with unsupported name: %s', filename);
        continue;
    end

    loaded_data = load(fullfile(folder, filename));
    if isfield(loaded_data, 'source_ts')
        data = loaded_data.source_ts;
    elseif isfield(loaded_data, 'data')
        data = loaded_data.data;
    else
        warning('No source_ts/data field in %s', filename);
        continue;
    end

    if isfield(loaded_data, 'labels')
        labels = loaded_data.labels;
    else
        labels = {};
    end

    if cfg.run_adf
        [n_regions, ~, n_epochs] = size(data);
        non_stationary = 0;
        total_traces = n_regions * n_epochs;
        for region = 1:n_regions
            for epoch = 1:n_epochs
                ts = squeeze(data(region, :, epoch));
                if ~adftest(ts)
                    non_stationary = non_stationary + 1;
                end
            end
        end
        fprintf('   ADF non-stationary traces: %d/%d\n', non_stationary, total_traces);
    end

    if exist('demean', 'file')
        data = demean(data, true);
    else
        data = data - mean(data, 2, 'omitnan');
    end

    [moaic, mobic] = tsdata_to_varmo(data, cfg.momax, cfg.moregmode, [], [], []);
    morder = moselect(sprintf('VAR model order selection, max = %d', cfg.momax), ...
        cfg.mosel, 'AIC', moaic, 'BIC', mobic);
    [A, V] = tsdata_to_var(data, morder, cfg.regmode);
    var_info(A, V);

    F = var_to_pwcgc(A, V);
    tstat = var_to_pwcgc_tstat(data, V, morder, cfg.regmode, cfg.test_stat);

    nvars = size(A, 1);
    nobs = size(data, 2);
    ntrials = size(data, 3);
    pval = mvgc_pval(tstat, cfg.test_stat, 1, 1, nvars-2, morder, nobs, ntrials);
    sig = significance(pval, cfg.alpha, cfg.mhtc);

    eweights = F / nanmax(F(:));
    file_tag = sprintf('%s_%s', parsed.subject_id, normalize_condition(parsed.condition));
    eweights_file = fullfile(cfg.pwcgc_dir, sprintf('eweights_matrix_%s.mat', file_tag));
    save(eweights_file, 'eweights', 'atlas');

    [Ass, Vss] = transform_var(A, V);
    [A_ss, C_ss, K_ss] = var_to_ss(Ass);
    CAK = Ass;
    H = var2trfun(Ass, cfg.fres);

    modfile = fullfile(cfg.pwcgc_dir, sprintf('pwcgc_params_%s.mat', file_tag));
    save(modfile, ...
        'V', 'CAK', 'H', 'Vss', 'A', 'A_ss', 'C_ss', 'K_ss', ...
        'morder', 'pval', 'sig', 'labels', 'atlas');
end

fprintf('Atlas "%s" MVGC complete. Output: %s\n', atlas, cfg.pwcgc_dir);
end

function parsed = parse_source_filename(filename)
parsed.valid = false;
parsed.subject_id = '';
parsed.condition = '';
tokens = regexp(filename, '^(EPCTL\d{2})_([A-Za-z0-9_]+)_source_time_series\.mat$', ...
    'tokens', 'once');
if isempty(tokens)
    return;
end
parsed.subject_id = tokens{1};
parsed.condition = tokens{2};
parsed.valid = true;
end

function out_condition = normalize_condition(in_condition)
if strcmpi(in_condition, 'W')
    out_condition = 'W_';
elseif strcmpi(in_condition, 'R')
    out_condition = 'R_';
else
    out_condition = in_condition;
end
end
