function sd_gc01_compute_subject_measures_by_stage(cfg)
% Stage 1: estimate subject/stage VAR models and GC/MI measures using MVGC2.

if nargin < 1
    cfg = sd_gc_config_anatomical7();
end

fprintf('\n[Stage 1] Subject-level model fitting and measure computation\n');

tasks = discover_tasks(cfg);
if isempty(tasks)
    if ~cfg.recompute_subject_measures
        fprintf('  No pending files (recompute_subject_measures=false). Stage 1 skipped.\n');
        return;
    end
    error('Stage 1 found no input files under %s', cfg.paths.source_root);
end

fprintf('  Candidate files: %d\n', numel(tasks));
if isfinite(cfg.max_files)
    tasks = tasks(1:min(numel(tasks), cfg.max_files));
    fprintf('  Debug max_files limit applied: %d\n', numel(tasks));
end

if cfg.parallel.enable && numel(tasks) > 1
    ensure_parallel_pool(cfg.parallel.n_workers);
    records = run_parallel(tasks, cfg);
else
    records = run_serial(tasks, cfg);
end

T = struct2table(records, 'AsArray', true);
out_csv = fullfile(cfg.paths.logs, 'stage1_subject_measure_status.csv');
writetable(T, out_csv);
fprintf('  Stage 1 status table saved: %s\n', out_csv);

n_ok = sum([records.success]);
if n_ok == 0
    error('Stage 1 completed with zero successful subject measure files.');
end
if n_ok < numel(records)
    fprintf(2, '  WARNING: Stage 1 had %d failures out of %d files.\n', numel(records)-n_ok, numel(records));
end

end

function tasks = discover_tasks(cfg)
files = dir(fullfile(cfg.paths.source_root, 'EPCTL*', 'mat', '*_source_time_series.mat'));

tasks = struct('subject', {}, 'stage', {}, 'stage_label', {}, 'input_file', {}, 'output_file', {});
for k = 1:numel(files)
    parsed = parse_source_filename(files(k).name);
    if ~parsed.valid || ~ismember(parsed.stage, cfg.stages)
        continue;
    end

    out_dir = fullfile(cfg.paths.subject_measures, parsed.stage);
    out_file = fullfile(out_dir, sprintf('%s_subject_measures.mat', parsed.subject));

    if ~cfg.recompute_subject_measures && exist(out_file, 'file') == 2
        continue;
    end

    t = struct();
    t.subject = parsed.subject;
    t.stage = parsed.stage;
    t.stage_label = cfg.stage_display.(parsed.stage);
    t.input_file = fullfile(files(k).folder, files(k).name);
    t.output_file = out_file;
    tasks(end+1) = t; %#ok<AGROW>
end

% Stable deterministic ordering by stage then subject.
if ~isempty(tasks)
    [~, idx] = sortrows([{tasks.stage}', {tasks.subject}']);
    tasks = tasks(idx);
end

end

function ensure_parallel_pool(n_workers)
try
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(n_workers)
            parpool('local');
        else
            parpool('local', n_workers);
        end
    else
        if ~isempty(n_workers) && pool.NumWorkers ~= n_workers
            delete(pool);
            parpool('local', n_workers);
        end
    end
catch ME
    warning('Parallel pool setup failed, falling back to serial execution: %s', ME.message);
end
end

function records = run_parallel(tasks, cfg)
n = numel(tasks);
records = repmat(init_record(), n, 1);
parfor i = 1:n
    records(i) = process_single_task(tasks(i), cfg);
end
end

function records = run_serial(tasks, cfg)
n = numel(tasks);
records = repmat(init_record(), n, 1);
for i = 1:n
    records(i) = process_single_task(tasks(i), cfg);
end
end

function rec = init_record()
rec = struct(...
    'subject', '', ...
    'stage', '', ...
    'stage_label', '', ...
    'input_file', '', ...
    'output_file', '', ...
    'success', false, ...
    'message', '', ...
    'nreg', NaN, ...
    'ntime', NaN, ...
    'nep', NaN, ...
    'morder', NaN, ...
    'fs_target', NaN ...
);
end

function rec = process_single_task(task, cfg)
rec = init_record();
rec.subject = task.subject;
rec.stage = task.stage;
rec.stage_label = task.stage_label;
rec.input_file = task.input_file;
rec.output_file = task.output_file;
rec.fs_target = cfg.fs_target;

fprintf('    [%s | %s] computing...\n', task.subject, task.stage);

try
    subject_measure = compute_subject_measures(task.input_file, task.subject, task.stage, task.stage_label, cfg);
    save(task.output_file, 'subject_measure', '-v7.3');

    rec.success = true;
    rec.nreg = subject_measure.meta.nreg;
    rec.ntime = subject_measure.meta.ntime;
    rec.nep = subject_measure.meta.nep;
    rec.morder = subject_measure.model.morder;

    fprintf('      OK [%s | %s] (morder=%d, ntime=%d, nep=%d)\n', ...
        task.subject, task.stage, rec.morder, rec.ntime, rec.nep);
catch ME
    rec.message = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(2, '      FAILED [%s | %s]: %s\n', task.subject, task.stage, ME.message);
end
end

function subject_measure = compute_subject_measures(in_file, subject_id, stage_id, stage_label, cfg)
S = load(in_file);
if isfield(S, 'source_ts')
    X = double(S.source_ts);
elseif isfield(S, 'data')
    X = double(S.data);
else
    error('Missing source_ts/data in %s', in_file);
end

if isfield(S, 'labels')
    labels = normalize_labels(S.labels);
else
    error('Missing labels in %s', in_file);
end

[nreg, ~, nep0] = size(X);
if nreg < 3
    error('Need at least 3 variables for conditional GC; got %d', nreg);
end

% Downsample then normalize in multi-trial/panel style.
dsfac = round(cfg.fs_original / cfg.fs_target);
if abs(cfg.fs_original / cfg.fs_target - dsfac) > 1e-9
    error('Non-integer downsample ratio: %.6f', cfg.fs_original / cfg.fs_target);
end
X = X(:, 1:dsfac:end, :);
[X, norm_meta] = apply_normalization(X, cfg);

[nreg, ntime, nep] = size(X);
groups = num2cell(1:nreg);
group_names = labels;

[aic, bic, hqc, lrt] = tsdata_to_varmo(X, cfg.momax, cfg.moregmode, [], [], []);
morder = moselect('somnodynamics-stage morder', cfg.mosel, ...
    'AIC', aic, 'BIC', bic, 'HQC', hqc, 'LRT', lrt);
morder_hqc = moselect('somnodynamics-stage morder hqc-check', 'HQC', ...
    'AIC', aic, 'BIC', bic, 'HQC', hqc, 'LRT', lrt);
[A, V] = tsdata_to_var(X, morder, cfg.regmode);
if isempty(A)
    error('VAR estimation failed for %s', in_file);
end

freqs = sfreqs(cfg.fres, cfg.fs_target);
band_info = build_band_info(freqs, cfg.band_defs);

meas_dir_t.src_in = var_to_gwiogc(A, V, groups, 'in');
meas_dir_t.src_out = var_to_gwiogc(A, V, groups, 'out');
meas_dir_t.roi = var_to_gwcgc(A, V, groups);

meas_dir_f.src_in = var_to_sgwiogc(A, V, groups, 'in', cfg.fres);
meas_dir_f.src_out = var_to_sgwiogc(A, V, groups, 'out', cfg.fres);
meas_dir_f.roi = var_to_sgwcgc(A, V, groups, cfg.fres);

[G, ~] = var_to_autocov(A, V, 0);
if isempty(G)
    error('var_to_autocov returned empty sequence for %s', in_file);
end
C0 = G(:, :, 1);

meas_undir_t.src = cov_to_gwiomi(C0, groups);
meas_undir_t.roi = cov_to_gwcmi(C0, groups);

[Sxx, ~] = var_to_cpsd(A, V, cfg.fres);
meas_undir_f.src = cpsd_to_sgwiomi(Sxx, groups);
meas_undir_f.roi = cpsd_to_sgwcmi(Sxx, groups);

meas_dir_t.roi = set_diag_nan(meas_dir_t.roi);
meas_undir_t.roi = set_diag_nan(meas_undir_t.roi);
meas_dir_f.roi = set_diag_nan_3d(meas_dir_f.roi);
meas_undir_f.roi = set_diag_nan_3d(meas_undir_f.roi);

estimator_gc = struct();
if cfg.compute_estimator_level_gc_sig
    tstat = var_to_pwcgc_tstat(X, V, morder, cfg.regmode, cfg.estimator_test_stat);
    pval = mvgc_pval(tstat, cfg.estimator_test_stat, 1, 1, nreg - 2, morder, ntime, nep);
    sig = significance(pval, cfg.alpha_estimator, cfg.estimator_mhtc);
    sig(1:nreg+1:end) = false;

    estimator_gc.test_stat = cfg.estimator_test_stat;
    estimator_gc.alpha = cfg.alpha_estimator;
    estimator_gc.correction = cfg.estimator_mhtc;
    estimator_gc.tstat = tstat;
    estimator_gc.pval = pval;
    estimator_gc.sig = sig;
end

subject_measure = struct();
subject_measure.meta = struct(...
    'input_file', in_file, ...
    'subject_id', subject_id, ...
    'stage_id', stage_id, ...
    'stage_label', stage_label, ...
    'nreg', nreg, ...
    'ntime', ntime, ...
    'nep', nep, ...
    'nep_original', nep0, ...
    'labels', {labels}, ...
    'group_names', {group_names}, ...
    'groups', {groups}, ...
    'directed_convention', cfg.directed_convention, ...
    'fs_original', cfg.fs_original, ...
    'fs_target', cfg.fs_target, ...
    'fres', cfg.fres, ...
    'freqs', freqs, ...
    'band_info', band_info, ...
    'normalization_mode', cfg.normalization_mode, ...
    'normalization', norm_meta ...
);

subject_measure.model = struct(...
    'morder', morder, ...
    'morder_hqc', morder_hqc, ...
    'aic', aic, ...
    'bic', bic, ...
    'hqc', hqc, ...
    'lrt', lrt, ...
    'mosel_used', cfg.mosel, ...
    'A', A, ...
    'V', V ...
);

subject_measure.measures = struct();
subject_measure.measures.directed = struct('time', meas_dir_t, 'spectral', meas_dir_f);
subject_measure.measures.undirected = struct('time', meas_undir_t, 'spectral', meas_undir_f);
subject_measure.estimator_level_pairwise_gc = estimator_gc;

end

function parsed = parse_source_filename(filename)
parsed = struct('valid', false, 'subject', '', 'stage', '');

tok = regexp(filename, '^(EPCTL\d{2})_([A-Za-z0-9_]+)_source_time_series(?:_46)?\.mat$', 'tokens', 'once');
if isempty(tok)
    return;
end

stage = normalize_stage_token(tok{2});
if isempty(stage)
    return;
end

parsed.valid = true;
parsed.subject = tok{1};
parsed.stage = stage;
end

function stage = normalize_stage_token(raw)
t = upper(string(raw));
if startsWith(t, "W")
    stage = 'W_';
    return;
end
if startsWith(t, "R")
    stage = 'R_';
    return;
end
if t == "N1"
    stage = 'N1';
elseif t == "N2"
    stage = 'N2';
elseif t == "N3"
    stage = 'N3';
else
    stage = '';
end
end

function labels = normalize_labels(raw_labels)
if isstring(raw_labels)
    labels = cellstr(raw_labels(:));
elseif ischar(raw_labels)
    labels = cellstr(raw_labels);
elseif iscell(raw_labels)
    labels = raw_labels(:);
else
    labels = cellstr(string(raw_labels(:)));
end
labels = cellfun(@(x) strtrim(char(x)), labels, 'UniformOutput', false);
end

function [X, norm_meta] = apply_normalization(X, cfg)
mode = lower(string(cfg.normalization_mode));
[nreg, ~, ~] = size(X);

norm_meta = struct();
norm_meta.mode = char(mode);
norm_meta.eps = cfg.normalization_eps;
norm_meta.n_zero_std = 0;
norm_meta.pooled_mean = nan(nreg, 1);
norm_meta.pooled_std = nan(nreg, 1);

switch mode
    case "pooled_zscore"
        % Barnett-style panel normalization: one mean/std per region and subject-stage.
        for r = 1:nreg
            v = reshape(X(r, :, :), 1, []);
            mu = mean(v, 'omitnan');
            sd = std(v, 0, 'omitnan');

            if ~isfinite(mu)
                mu = 0;
            end
            if ~isfinite(sd) || sd < cfg.normalization_eps
                sd = 1;
                norm_meta.n_zero_std = norm_meta.n_zero_std + 1;
            end

            X(r, :, :) = (X(r, :, :) - mu) ./ sd;
            norm_meta.pooled_mean(r) = mu;
            norm_meta.pooled_std(r) = sd;
        end

    case "demean_epoch"
        X = X - mean(X, 2, 'omitnan');
        norm_meta.pooled_mean = [];
        norm_meta.pooled_std = [];

    case "none"
        norm_meta.pooled_mean = [];
        norm_meta.pooled_std = [];

    otherwise
        error('Unsupported normalization_mode: %s', cfg.normalization_mode);
end
end

function band_info = build_band_info(freqs, band_defs)
nyq = max(freqs);
n_def = size(band_defs, 1);

tmp = struct('name', {}, 'lo', {}, 'hi', {}, 'mask', {}, 'n_bins', {});
for b = 1:n_def
    name = band_defs{b, 1};
    lo = band_defs{b, 2};
    hi_raw = band_defs{b, 3};
    hi = min(hi_raw, nyq);
    if hi <= lo
        continue;
    end

    mask = (freqs >= lo) & (freqs <= hi);
    if ~any(mask)
        continue;
    end

    s = struct();
    s.name = name;
    s.lo = lo;
    s.hi = hi;
    s.mask = mask;
    s.n_bins = sum(mask);
    tmp(end+1) = s; %#ok<AGROW>
end

if isempty(tmp)
    error('No valid band definitions after Nyquist clipping.');
end
band_info = tmp;
end

function M = set_diag_nan(M)
n = size(M, 1);
M(1:n+1:end) = NaN;
end

function M = set_diag_nan_3d(M)
n = size(M, 1);
for k = 1:size(M, 3)
    Mk = M(:, :, k);
    Mk(1:n+1:end) = NaN;
    M(:, :, k) = Mk;
end
end
