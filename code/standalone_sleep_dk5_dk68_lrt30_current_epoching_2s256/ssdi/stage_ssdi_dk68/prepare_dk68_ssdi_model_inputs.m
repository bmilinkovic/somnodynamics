function prepare_dk68_ssdi_model_inputs(varargin)
%PREPARE_DK68_SSDI_MODEL_INPUTS Fit DK68 VAR models for SSDI without full DK68 GC.
%
% This writes subject_measure-compatible files containing model.A/model.V/model.morder
% and meta labels. The DK68 SSDI loader only consumes those model fields, so this
% avoids computing/storing full 68x68 GC/MI products that are unnecessary for SSDI.

p = inputParser;
results_root = sd_external_results_root(project_root());
addParameter(p, 'source_root', fullfile(results_root, 'source_recon', 'dk5_dk68_lcmv_nai_current_epoching_2s_256hz', 'dk68'));
addParameter(p, 'output_root', fullfile(results_root, 'ssdi', 'stage_ssdi_dk68_lrt30_current_epoching_2s_256hz', 'model_inputs', 'subject_measures'));
addParameter(p, 'fs_original', 256);
addParameter(p, 'fs_target', 256);
addParameter(p, 'momax', 30);
addParameter(p, 'mosel', 'LRT');
addParameter(p, 'moregmode', 'LWR');
addParameter(p, 'regmode', 'LWR');
addParameter(p, 'normalization_eps', 1e-12);
addParameter(p, 'parallel_enable', true);
addParameter(p, 'parallel_workers', 4);
addParameter(p, 'recompute', false);
addParameter(p, 'max_files', inf);
parse(p, varargin{:});
cfg = p.Results;

repo = project_root();
addpath(fullfile(repo, 'code', 'run_all'));
setup_matlab_toolboxes(repo);

stages = {'W_', 'N1', 'N2', 'N3', 'R_'};
for i = 1:numel(stages)
    od = fullfile(cfg.output_root, stages{i});
    if exist(od, 'dir') ~= 7
        mkdir(od);
    end
end
log_dir = fullfile(fileparts(fileparts(cfg.output_root)), 'logs');
if exist(log_dir, 'dir') ~= 7
    mkdir(log_dir);
end
status_file = fullfile(log_dir, 'dk68_ssdi_model_input_status.csv');

files = dir(fullfile(cfg.source_root, 'EPCTL*', 'mat', '*_source_time_series_68.mat'));
tasks = struct('subject', {}, 'stage', {}, 'input_file', {}, 'output_file', {});
for k = 1:numel(files)
    parsed = parse_source_filename(files(k).name);
    if ~parsed.valid || ~ismember(parsed.stage, stages)
        continue;
    end
    out_file = fullfile(cfg.output_root, parsed.stage, sprintf('%s_subject_measures.mat', parsed.subject));
    if ~cfg.recompute && exist(out_file, 'file') == 2
        continue;
    end
    t = struct();
    t.subject = parsed.subject;
    t.stage = parsed.stage;
    t.input_file = fullfile(files(k).folder, files(k).name);
    t.output_file = out_file;
    tasks(end+1) = t; %#ok<AGROW>
end
if ~isempty(tasks)
    [~, idx] = sortrows([{tasks.stage}', {tasks.subject}']);
    tasks = tasks(idx);
end
if isfinite(cfg.max_files)
    tasks = tasks(1:min(numel(tasks), cfg.max_files));
end

fprintf('[DK68 SSDI model inputs] source: %s\n', cfg.source_root);
fprintf('[DK68 SSDI model inputs] output: %s\n', cfg.output_root);
fprintf('[DK68 SSDI model inputs] pending tasks: %d\n', numel(tasks));
if isempty(tasks)
    return;
end

if cfg.parallel_enable && numel(tasks) > 1
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', cfg.parallel_workers);
    elseif pool.NumWorkers ~= cfg.parallel_workers
        delete(pool);
        parpool('local', cfg.parallel_workers);
    end
    records = repmat(init_record(), numel(tasks), 1);
    parfor i = 1:numel(tasks)
        records(i) = process_task(tasks(i), cfg);
    end
else
    records = repmat(init_record(), numel(tasks), 1);
    for i = 1:numel(tasks)
        records(i) = process_task(tasks(i), cfg);
    end
end

T = struct2table(records, 'AsArray', true);
writetable(T, status_file);
fprintf('[DK68 SSDI model inputs] status: %s\n', status_file);
if any(~[records.success])
    warning('%d/%d model-input tasks failed.', sum(~[records.success]), numel(records));
end
end

function rec = process_task(task, cfg)
rec = init_record();
rec.subject = task.subject;
rec.stage = task.stage;
rec.input_file = task.input_file;
rec.output_file = task.output_file;
try
    fprintf('  [%s | %s] fitting VAR...\n', task.subject, task.stage);
    S = load(task.input_file);
    if isfield(S, 'source_ts')
        X = double(S.source_ts);
    elseif isfield(S, 'data')
        X = double(S.data);
    else
        error('Missing source_ts/data in %s', task.input_file);
    end
    if isfield(S, 'labels')
        labels = normalize_labels(S.labels);
    else
        labels = arrayfun(@(i) sprintf('DK%02d', i), 1:size(X,1), 'UniformOutput', false);
    end
    if size(X,1) ~= 68
        error('Expected 68 regions, got %d', size(X,1));
    end
    dsfac = round(cfg.fs_original / cfg.fs_target);
    if abs(cfg.fs_original / cfg.fs_target - dsfac) > 1e-9
        error('Non-integer downsample ratio.');
    end
    X = X(:, 1:dsfac:end, :);
    [X, norm_meta] = pooled_zscore(X, cfg.normalization_eps);
    [nreg, ntime, nep] = size(X);

    [aic, bic, hqc, lrt] = tsdata_to_varmo(X, cfg.momax, cfg.moregmode, [], [], []);
    morder = moselect('dk68-ssdi-model-input morder', cfg.mosel, 'AIC', aic, 'BIC', bic, 'HQC', hqc, 'LRT', lrt);
    morder_hqc = moselect('dk68-ssdi-model-input hqc-check', 'HQC', 'AIC', aic, 'BIC', bic, 'HQC', hqc, 'LRT', lrt);
    [A, V] = tsdata_to_var(X, morder, cfg.regmode);
    if isempty(A)
        error('VAR estimation failed');
    end

    subject_measure = struct();
    subject_measure.meta = struct( ...
        'input_file', task.input_file, ...
        'subject_id', task.subject, ...
        'stage_id', task.stage, ...
        'nreg', nreg, ...
        'ntime', ntime, ...
        'nep', nep, ...
        'labels', {labels}, ...
        'group_names', {labels}, ...
        'fs_original', cfg.fs_original, ...
        'fs_target', cfg.fs_target, ...
        'normalization_mode', 'pooled_zscore', ...
        'normalization', norm_meta, ...
        'source_type', 'dk68_source_recon_var_model_only_for_ssdi');
    subject_measure.model = struct( ...
        'A', A, ...
        'V', V, ...
        'morder', morder, ...
        'morder_hqc', morder_hqc, ...
        'momax', cfg.momax, ...
        'mosel', cfg.mosel, ...
        'moregmode', cfg.moregmode, ...
        'regmode', cfg.regmode, ...
        'aic', aic, ...
        'bic', bic, ...
        'hqc', hqc, ...
        'lrt', lrt);
    save(task.output_file, 'subject_measure', '-v7.3');

    rec.success = true;
    rec.nreg = nreg;
    rec.ntime = ntime;
    rec.nep = nep;
    rec.morder = morder;
    fprintf('    OK [%s | %s] morder=%d ntime=%d nep=%d\n', task.subject, task.stage, morder, ntime, nep);
catch ME
    rec.message = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(2, '    FAILED [%s | %s]: %s\n', task.subject, task.stage, ME.message);
end
end

function rec = init_record()
rec = struct('subject','','stage','','input_file','','output_file','','success',false,'message','','nreg',NaN,'ntime',NaN,'nep',NaN,'morder',NaN);
end

function [X, meta] = pooled_zscore(X, epsv)
mu = mean(X, [2 3], 'omitnan');
sd = std(reshape(permute(X, [1 2 3]), size(X,1), []), 0, 2, 'omitnan');
sd(sd < epsv | ~isfinite(sd)) = 1;
X = (X - mu) ./ sd;
meta = struct('mean', mu, 'std', sd, 'eps', epsv);
end

function parsed = parse_source_filename(filename)
parsed = struct('valid', false, 'subject', '', 'stage', '');
tok = regexp(filename, '^(EPCTL\d{2})_([A-Za-z0-9_]+)_source_time_series_68\.mat$', 'tokens', 'once');
if isempty(tok), return; end
stage = normalize_stage_token(tok{2});
if isempty(stage), return; end
parsed.valid = true;
parsed.subject = tok{1};
parsed.stage = stage;
end

function stage = normalize_stage_token(raw)
t = upper(string(raw));
if startsWith(t, "W")
    stage = 'W_';
elseif startsWith(t, "R")
    stage = 'R_';
elseif t == "N1"
    stage = 'N1';
elseif t == "N2"
    stage = 'N2';
elseif t == "N3"
    stage = 'N3';
else
    stage = '';
end
end

function labels = normalize_labels(raw)
if isstring(raw)
    labels = cellstr(raw(:));
elseif ischar(raw)
    labels = {raw};
elseif iscell(raw)
    labels = cell(size(raw));
    for i = 1:numel(raw)
        labels{i} = char(string(raw{i}));
    end
    labels = labels(:)';
else
    labels = cellstr(string(raw));
    labels = labels(:)';
end
end

function r = project_root()
r = '/Users/borjan/code/somnodynamics';
end
