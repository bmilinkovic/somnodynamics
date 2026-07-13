function sd_ssdi01_compute_stage_spectral(cfg)
% Stage 1: run stage-wise SSDI spectral optimisation from GC subject measures.

if nargin < 1
    cfg = sd_ssdi_config_anatomical7();
end

fprintf('[Stage 1] Stage-wise SSDI spectral optimisation\n');

band_info = build_band_info(sfreqs(cfg.fres, cfg.fs), cfg.band_defs);
write_band_report(cfg, band_info);

tasks = discover_tasks(cfg);
if isempty(tasks)
    if ~cfg.recompute
        fprintf('  No pending files (recompute=false). Stage 1 skipped.\n');
        return;
    end
    error('No input tasks found under %s', cfg.paths.gc_subject_measures);
end

fprintf('  Candidate files: %d\n', numel(tasks));
if isfinite(cfg.max_files)
    tasks = tasks(1:min(numel(tasks), cfg.max_files));
    fprintf('  Debug max_files limit applied: %d\n', numel(tasks));
end

if cfg.parallel.enable && numel(tasks) > 1
    ensure_parallel_pool(cfg.parallel.n_workers);
    records = run_parallel(tasks, cfg, band_info);
else
    records = run_serial(tasks, cfg, band_info);
end

T = struct2table(records, 'AsArray', true);
out_csv = fullfile(cfg.paths.logs, 'stage1_ssdi_status.csv');
writetable(T, out_csv);
fprintf('  Stage 1 status table saved: %s\n', out_csv);

n_ok = sum([records.success]);
if n_ok == 0
    error('Stage 1 completed with zero successful SSDI files.');
end
if n_ok < numel(records)
    fprintf(2, '  WARNING: Stage 1 had %d failures out of %d files.\n', numel(records)-n_ok, numel(records));
end

end

function tasks = discover_tasks(cfg)
files = dir(fullfile(cfg.paths.gc_subject_measures, '*', '*_subject_measures.mat'));
tasks = struct('subject', {}, 'stage', {}, 'stage_label', {}, 'input_file', {}, 'output_file', {});
for k = 1:numel(files)
    subject = erase(files(k).name, '_subject_measures.mat');
    stage = string(split(string(files(k).folder), filesep));
    stage = char(stage(end));
    if ~ismember(stage, cfg.stages)
        continue;
    end

    out_file = fullfile(cfg.paths.spectral_out, sprintf('%s_%s_ssdi_spectral_domain_results.mat', subject, stage));
    if ~cfg.recompute && exist(out_file, 'file') == 2
        continue;
    end

    t = struct();
    t.subject = subject;
    t.stage = stage;
    t.stage_label = cfg.stage_display.(stage);
    t.input_file = fullfile(files(k).folder, files(k).name);
    t.output_file = out_file;
    tasks(end+1) = t; %#ok<AGROW>
end

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

function records = run_parallel(tasks, cfg, band_info)
n = numel(tasks);
records = repmat(init_record(), n, 1);
parfor i = 1:n
    records(i) = process_single_task(tasks(i), cfg, band_info);
end
end

function records = run_serial(tasks, cfg, band_info)
n = numel(tasks);
records = repmat(init_record(), n, 1);
for i = 1:n
    records(i) = process_single_task(tasks(i), cfg, band_info);
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
    'n_regions', NaN, ...
    'n_bands', NaN, ...
    'n_m', NaN ...
);
end

function rec = process_single_task(task, cfg, band_info)
rec = init_record();
rec.subject = task.subject;
rec.stage = task.stage;
rec.stage_label = task.stage_label;
rec.input_file = task.input_file;
rec.output_file = task.output_file;

fprintf('    [%s | %s] computing SSDI...\n', task.subject, task.stage);

try
    results = compute_ssdi_results(task, cfg, band_info);
    save(task.output_file, 'results');

    rec.success = true;
    rec.n_regions = results.meta.n_regions;
    rec.n_bands = numel(band_info);
    rec.n_m = numel(cfg.m_values);

    fprintf('      OK [%s | %s] (n=%d bands=%d m=%d)\n', ...
        task.subject, task.stage, rec.n_regions, rec.n_bands, rec.n_m);
catch ME
    rec.message = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(2, '      FAILED [%s | %s]: %s\n', task.subject, task.stage, ME.message);
end
end

function results = compute_ssdi_results(task, cfg, band_info)
S = load(task.input_file, 'subject_measure');
if ~isfield(S, 'subject_measure')
    error('Missing subject_measure in %s', task.input_file);
end
M = S.subject_measure;

if ~isfield(M, 'model') || ~isfield(M.model, 'A') || ~isfield(M.model, 'V')
    error('Missing model.A/model.V in %s', task.input_file);
end

A = double(M.model.A);
V = double(M.model.V);
labels = normalize_labels(M.meta.group_names);

n = size(A, 1);
if n ~= cfg.required_regions
    error('Expected %d regions for anatomical7, got %d in %s', cfg.required_regions, n, task.input_file);
end

[Ass, ~] = transform_var(A, V);
[A_ss, C_ss, K_ss] = var_to_ss(Ass);
H = ss2trfun(A_ss, C_ss, K_ss, cfg.fres);
freqs = sfreqs(cfg.fres, cfg.fs);

results = struct();
results.meta = struct();
results.meta.subject_id = task.subject;
results.meta.stage = task.stage;
results.meta.stage_label = task.stage_label;
results.meta.labels = labels;
results.meta.n_regions = n;
results.meta.m_values = cfg.m_values;
results.meta.fs = cfg.fs;
results.meta.fres = cfg.fres;
results.meta.freqs = freqs(:)';
results.meta.input_file = task.input_file;
results.meta.timestamp = datestr(now, 31);
results.meta.band_info = band_info;
results.meta.optimiser = 'opt_gd_dds_mruns';
results.meta.optimiser_objective = 'band_limited_spectral_dd';
results.meta.nrunsp = cfg.nrunsp;
results.meta.nitersp = cfg.nitersp;

for m = cfg.m_values
    if m >= n
        continue;
    end

    for b = 1:numel(band_info)
        bname = band_info(b).name;
        idx = band_info(b).indices;
        if isempty(idx)
            continue;
        end

        H_band = H(:,:,idx);

        seed = cfg.iseed + sum(double(task.subject)) + 1000*m + 100*b + sum(double(task.stage));
        rs = rng_seed(seed);
        L0 = rand_orthonormal(n, m, cfg.nrunsp);
        rng_restore(rs);

        [dopto, Lo, conv, iopt, sopt, cput, ohist] = opt_gd_dds_mruns( ...
            H_band, L0, cfg.nitersp, cfg.gdesp, cfg.sig0p, cfg.gdlsp, cfg.gdtolp, cfg.histp, cfg.ppp);

        Lopto = itransform_subspace(Lo, V);
        gopto = gmetrics(Lopto);
        [node_weights, beta_stat] = compute_node_metrics(Lopto);

        results.(sprintf('%s_Lopto_m%d', bname, m)) = Lopto;
        results.(sprintf('%s_dopto_m%d', bname, m)) = dopto;
        results.(sprintf('%s_gopto_m%d', bname, m)) = gopto;
        results.(sprintf('%s_node_weights_m%d', bname, m)) = node_weights;
        results.(sprintf('%s_beta_m%d', bname, m)) = beta_stat;
        results.(sprintf('%s_opt_m%d', bname, m)) = struct( ...
            'conv', conv, 'iopt', iopt, 'sopt', sopt, 'cput', cput, 'ohist', ohist, ...
            'freq_idx', idx(:)', ...
            'freq_lo', band_info(b).lo, ...
            'freq_hi', band_info(b).hi);

        cell_key = sprintf('%s_Lopto_m', bname);
        if ~isfield(results, cell_key)
            results.(cell_key) = cell(1, max(cfg.m_values)-1);
        end
        tmp = results.(cell_key);
        tmp{m-1} = Lopto;
        results.(cell_key) = tmp;
    end
end
end

function [node_weights, beta_stat] = compute_node_metrics(Lopto)
n = size(Lopto,1);
r = size(Lopto,3);
node_weights = zeros(n, r);
beta_stat = zeros(n, r);
for k = 1:r
    node_weights(:,k) = 1 - gmetricsx(Lopto(:,:,k));
    if exist('habeta', 'file') == 2
        bvec = double(real(habeta(Lopto(:,:,k))));
        bvec(~isfinite(bvec)) = 0;
        beta_stat(:,k) = min(max(bvec, 0), 1);
    else
        beta_stat(:,k) = nan(n,1);
    end
end
end

function labels = normalize_labels(raw)
if isstring(raw)
    labels = cellstr(raw(:));
elseif ischar(raw)
    labels = cellstr(raw);
elseif iscell(raw)
    labels = raw(:);
else
    labels = cellstr(string(raw(:)));
end
labels = cellfun(@(x) strtrim(char(x)), labels, 'UniformOutput', false);
end

function band_info = build_band_info(freqs, band_defs)
band_info = struct('name', {}, 'lo', {}, 'hi', {}, 'indices', {}, 'n_bins', {});
for b = 1:numel(band_defs)
    lo = band_defs(b).lo;
    hi = band_defs(b).hi;
    idx = find(freqs >= lo & freqs <= hi);
    if isempty(idx)
        continue;
    end

    s = struct();
    s.name = band_defs(b).name;
    s.lo = lo;
    s.hi = hi;
    s.indices = idx(:)';
    s.n_bins = numel(idx);
    band_info(end+1) = s; %#ok<AGROW>
end
if isempty(band_info)
    error('No valid frequency bands after clipping to available freqs.');
end
end

function write_band_report(cfg, band_info)
f = fullfile(cfg.paths.diag, 'frequency_band_summary.txt');
fid = fopen(f, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'SomnoDynamics SSDI Stage Bands\n');
fprintf(fid, '=============================\n');
fprintf(fid, 'fs=%g, fres=%g\n\n', cfg.fs, cfg.fres);
for i = 1:numel(band_info)
    b = band_info(i);
    fprintf(fid, '%s: %g-%g Hz (%d bins)\n', b.name, b.lo, b.hi, b.n_bins);
end
end
