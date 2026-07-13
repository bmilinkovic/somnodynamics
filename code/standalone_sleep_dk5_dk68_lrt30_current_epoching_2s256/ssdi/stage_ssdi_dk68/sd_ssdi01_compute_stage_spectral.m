function sd_ssdi01_compute_stage_spectral(cfg)
% Stage 1: DMT-faithful SSDI spectral optimisation for sleep stages.

if nargin < 1
    cfg = sd_ssdi_config_anatomical7();
end

fprintf('[Stage 1] DMT-faithful stage-wise SSDI spectral optimisation\n');

freqs = sfreqs(cfg.fres, cfg.fs);
band_info = build_band_info(freqs, cfg.band_defs);
write_frequency_diagnostics(cfg, freqs, band_info);

jobs = discover_jobs(cfg);
if isempty(jobs)
    if ~cfg.recompute
        fprintf('  No pending files (recompute=false). Stage 1 skipped.\n');
        return;
    end
    error('No input tasks found under %s', cfg.subject_measure_root);
end

fprintf('  Candidate jobs: %d\n', size(jobs, 1));
if isfinite(cfg.max_files)
    jobs = jobs(1:min(size(jobs, 1), cfg.max_files), :);
    fprintf('  max_files limit applied: %d\n', size(jobs, 1));
end

maybe_start_parpool(cfg);

records = repmat(init_record(), size(jobs, 1), 1);
for i = 1:size(jobs, 1)
    subject_id = jobs{i, 1};
    stage = jobs{i, 2};
    records(i) = process_job(subject_id, stage, cfg, band_info);
end

T = struct2table(records, 'AsArray', true);
out_csv = fullfile(cfg.logs_dir, 'stage1_ssdi_status.csv');
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

function rec = process_job(subject_id, stage, cfg, band_info)
rec = init_record();
rec.subject = subject_id;
rec.stage = stage;
rec.stage_label = cfg.stage_display.(stage);
rec.input_file = fullfile(cfg.subject_measure_root, stage, sprintf('%s_subject_measures.mat', subject_id));
rec.output_file = fullfile(stage_out_dir(cfg, stage), sprintf('%s_ssdi_spectral_domain_results.mat', subject_id));
rec.nrunsp = cfg.nrunsp;
rec.nitersp = cfg.nitersp;
rec.niterso = cfg.niterso;

fprintf('    [%s | %s] computing SSDI...\n', subject_id, stage);
t0 = tic;

try
    [results, debug_rows] = compute_ssdi_results(subject_id, stage, cfg, band_info);

    if cfg.save_v73
        save(rec.output_file, 'results', '-v7.3');
    else
        save(rec.output_file, 'results');
    end

    append_dd_debug(cfg, debug_rows);
    append_status(cfg, subject_id, stage, true, 'ok');

    rec.success = true;
    rec.runtime_sec = toc(t0);
    rec.n_regions = results.meta.n_regions;
    rec.n_bands = numel(band_info);
    rec.n_m = numel(cfg.m_values);

    fprintf('      OK [%s | %s] (n=%d bands=%d m=%d, %.1fs)\n', ...
        subject_id, stage, rec.n_regions, rec.n_bands, rec.n_m, rec.runtime_sec);
catch ME
    rec.runtime_sec = toc(t0);
    rec.message = getReport(ME, 'extended', 'hyperlinks', 'off');
    append_status(cfg, subject_id, stage, false, ME.message);
    fprintf(2, '      FAILED [%s | %s]: %s\n', subject_id, stage, ME.message);
end

end

function [results, debug_rows] = compute_ssdi_results(subject_id, stage, cfg, band_info)
state = load_subject_state(subject_id, stage, cfg);

n = size(state.A, 1);
if cfg.strict_n7 && n ~= cfg.required_regions
    error('Expected n=%d regions, got n=%d for %s %s', cfg.required_regions, n, subject_id, stage);
end

fprintf('      Loaded model source: %s\n', state.source_type);
fprintf('      Input source file: %s\n', state.input_file);
fprintf('      Model size n=%d, morder=%d\n', n, state.morder);

[Ass, ~] = transform_var(state.A, state.V);
[A_ss, C_ss, K_ss] = var_to_ss(Ass);
CAK = iss2cak(A_ss, C_ss, K_ss);
spec_radius = max(abs(eig(A_ss)));
if spec_radius >= 1
    warning('State-space model appears unstable (spectral radius >= 1).');
end

freqs = sfreqs(cfg.fres, cfg.fs);
H = ss2trfun(A_ss, C_ss, K_ss, cfg.fres);

results = struct();
results.meta = struct();
results.meta.subject_id = subject_id;
results.meta.stage = stage;
results.meta.stage_label = cfg.stage_display.(stage);
results.meta.labels = state.labels;
results.meta.n_regions = n;
results.meta.morder = state.morder;
results.meta.fres = cfg.fres;
results.meta.fs = cfg.fs;
results.meta.freqs = freqs(:)';
results.meta.input_file = state.input_file;
results.meta.model_source = state.source_type;
results.meta.timestamp = datestr(now, 31);
results.meta.band_info = band_info;
results.meta.stability = struct('spectral_radius_A_ss', spec_radius, 'is_stable', spec_radius < 1);
results.meta.preopt_is_band_specific = true;
results.meta.preopt_note = [ ...
    "Pre-optimisation uses CAK-based proxy DD (canonical_opt_ddx_mruns); " ...
    "final optimisation uses canonical spectral DD Eq25/26 (canonical_opt_dds_mruns)."];
results.meta.pre_optimiser = 'canonical_opt_ddx_mruns';
results.meta.pre_optimiser_objective = 'proxy_ddx_cak';
results.meta.final_optimiser = 'canonical_opt_dds_mruns';
results.meta.final_optimiser_objective = 'canonical_band_dd_eq25_26';
results.meta.optimiser = 'canonical_opt_dds_mruns';
results.meta.optimiser_objective = 'canonical_band_dd_eq25_26';
results.meta.nrunsp = cfg.nrunsp;
results.meta.nitersp = cfg.nitersp;
results.meta.niterso = cfg.niterso;

debug_rows = struct([]);
tiny_tol = 1e-10;

for m = cfg.m_values
    if m >= n
        fprintf('      Skipping m=%d (must be < n=%d).\n', m, n);
        continue;
    end

    for b = 1:numel(band_info)
        bname = band_info(b).name;
        fidx = band_info(b).indices;

        if isempty(fidx)
            fprintf('      Skipping %s for m=%d (empty band).\n', bname, m);
            continue;
        end

        H_band = H(:,:,fidx);
        ranks = zeros(1, numel(fidx));
        for kk = 1:numel(fidx)
            ranks(kk) = rank(H_band(:,:,kk));
        end

        fprintf('      Optimising %s m=%d (%d bins, rank min/max=%d/%d)\n', ...
            bname, m, numel(fidx), min(ranks), max(ranks));

        % Band-specific pre-optimisation.
        seed = cfg.iseed + sum(double(subject_id)) + 1000*m + 100*b + sum(double(stage));
        rstate = rng_seed(seed);
        L0p = rand_orthonormal(n, m, cfg.nrunsp);
        rng_restore(rstate);

        [doptp, Lp, convp, ioptp, soptp, cputp, ohistp] = canonical_opt_ddx_mruns( ...
            CAK, L0p, cfg.nitersp, cfg.gdesp, cfg.sig0p, cfg.gdlsp, cfg.gdtolp, cfg.histp, cfg.ppp, 'pre');

        Loptp = itransform_subspace(Lp, state.V);
        goptp = gmetrics(Loptp);

        [uidx, usiz, nrunso] = Lcluster(goptp, cfg.ctol, doptp);
        if isempty(uidx)
            warning('No clusters returned at %s m=%d; falling back to all runs.', bname, m);
            uidx = 1:size(Lp,3);
            usiz = ones(size(uidx));
            nrunso = numel(uidx);
        end
        L0o = Lp(:,:,uidx);

        pre = struct();
        pre.method = 'canonical_opt_ddx_mruns';
        pre.objective = 'proxy_ddx_cak';
        pre.band = bname;
        pre.band_freq_idx = fidx(:)';
        pre.band_freq_vals = freqs(fidx)';
        pre.Loptp = Loptp;
        pre.Lp = Lp;
        pre.doptp = doptp;
        pre.goptp = goptp;
        pre.uidx = uidx;
        pre.usiz = usiz;
        pre.nrunso = nrunso;
        pre.convp = convp;
        pre.ioptp = ioptp;
        pre.soptp = soptp;
        pre.cputp = cputp;
        if cfg.histp
            pre.ohistp = ohistp;
        end

        results.(sprintf('%s_preopt_m%d', bname, m)) = pre;
        results.(sprintf('%s_Loptp_m%d', bname, m)) = Loptp;
        results.(sprintf('%s_doptp_m%d', bname, m)) = doptp;
        results.(sprintf('%s_goptp_m%d', bname, m)) = goptp;

        if b == 1 && ~isfield(results, sprintf('preopt_m%d', m))
            results.(sprintf('preopt_m%d', m)) = pre;
            results.(sprintf('Loptp_m%d', m)) = Loptp;
            results.(sprintf('doptp_m%d', m)) = doptp;
            results.(sprintf('goptp_m%d', m)) = goptp;
        end

        node_weights_pre = zeros(n, size(Loptp,3));
        beta_stat_pre = zeros(n, size(Loptp,3));
        for kp = 1:size(Loptp,3)
            node_weights_pre(:,kp) = 1 - gmetricsx(Loptp(:,:,kp));
            bvecp = double(real(habeta(Loptp(:,:,kp))));
            bvecp(~isfinite(bvecp)) = 0;
            beta_stat_pre(:,kp) = min(max(bvecp, 0), 1);
        end
        results.(sprintf('%s_pre_node_weights_m%d', bname, m)) = node_weights_pre;
        results.(sprintf('%s_pre_beta_m%d', bname, m)) = beta_stat_pre;
        if b == 1 && ~isfield(results, sprintf('pre_node_weights_m%d', m))
            results.(sprintf('pre_node_weights_m%d', m)) = node_weights_pre;
            results.(sprintf('pre_beta_m%d', m)) = beta_stat_pre;
        end

        [dopto_raw, Lo, convo, iopto, sopto, cputo, ohisto] = canonical_opt_dds_mruns( ...
            H, freqs, fidx, L0o, cfg.niterso, cfg.gdeso, cfg.sig0o, cfg.gdlso, cfg.gdtolo, cfg.histo, cfg.ppo, 'final');

        Lopto = itransform_subspace(Lo, state.V);
        gopto = gmetrics(Lopto);

        [dd_full, dd_full_raw] = trfun_dd_batch(Lo, H);
        dd_band = dopto_raw;
        dd_band_raw = dopto_raw;
        [dd_full, n_clip_full] = clip_tiny_negative(dd_full, tiny_tol);
        n_clip_band = 0;

        Lo_roundtrip = transform_subspace(Lopto, state.V);
        dd_roundtrip_band = canonical_band_dd_batch(Lo_roundtrip, H, freqs, fidx);

        node_weights = zeros(n, size(Lopto,3));
        beta_stat = zeros(n, size(Lopto,3));
        for k = 1:size(Lopto,3)
            node_weights(:,k) = 1 - gmetricsx(Lopto(:,:,k));
            bvec = double(real(habeta(Lopto(:,:,k))));
            bvec(~isfinite(bvec)) = 0;
            beta_stat(:,k) = min(max(bvec, 0), 1);
        end

        beta_best = beta_stat(:,1);
        [beta_cval, beta_pvals, beta_sigs] = habeta_statinf( ...
            beta_best, n, m, cfg.alpha, cfg.habeta_tails, cfg.habeta_mhtc);

        opt = struct();
        opt.method = 'canonical_opt_dds_mruns';
        opt.objective = 'canonical_band_dd_eq25_26';
        opt.band = bname;
        opt.band_freq_idx = fidx(:)';
        opt.band_freq_vals = freqs(fidx)';
        opt.dopto = dopto_raw;
        opt.convo = convo;
        opt.iopto = iopto;
        opt.sopto = sopto;
        opt.cputo = cputo;
        opt.niterso = cfg.niterso;
        opt.gdtolo = cfg.gdtolo;
        if cfg.histo
            opt.ohisto = ohisto;
        end

        results.(sprintf('%s_opt_m%d', bname, m)) = opt;
        results.(sprintf('%s_Lopto_m%d', bname, m)) = Lopto;
        results.(sprintf('%s_dopto_m%d', bname, m)) = dopto_raw;
        results.(sprintf('%s_dopto_raw_band_m%d', bname, m)) = dopto_raw;
        results.(sprintf('%s_dd_full_raw_m%d', bname, m)) = dd_full_raw;
        results.(sprintf('%s_dd_full_m%d', bname, m)) = dd_full;
        results.(sprintf('%s_dd_band_m%d', bname, m)) = dd_band;
        results.(sprintf('%s_dd_band_raw_m%d', bname, m)) = dd_band_raw;
        results.(sprintf('%s_gopto_m%d', bname, m)) = gopto;
        results.(sprintf('%s_node_weights_m%d', bname, m)) = node_weights;
        results.(sprintf('%s_beta_m%d', bname, m)) = beta_stat;
        results.(sprintf('%s_beta_cvals_m%d', bname, m)) = beta_cval;
        results.(sprintf('%s_beta_pvals_m%d', bname, m)) = beta_pvals;
        results.(sprintf('%s_beta_sigs_m%d', bname, m)) = beta_sigs;

        legacy_key = sprintf('%s_Lopto_m', bname);
        if ~isfield(results, legacy_key)
            results.(legacy_key) = cell(1, max(cfg.m_values));
        end
        tmp = results.(legacy_key);
        if numel(tmp) < m
            tmp{m} = [];
        end
        tmp{m} = Lopto;
        results.(legacy_key) = tmp;

        stage_mins = struct();
        stage_maxs = struct();
        stage_mins.preopt = min(doptp);
        stage_maxs.preopt = max(doptp);
        stage_mins.cluster = min(doptp(uidx));
        stage_maxs.cluster = max(doptp(uidx));
        stage_mins.second_stage_raw = min(dopto_raw);
        stage_maxs.second_stage_raw = max(dopto_raw);
        stage_mins.roundtrip_band = min(dd_roundtrip_band);
        stage_maxs.roundtrip_band = max(dd_roundtrip_band);
        stage_mins.full_dd_check = min(dd_full);
        stage_maxs.full_dd_check = max(dd_full);

        stage_flags = struct();
        stage_flags.preopt = any(doptp < -tiny_tol);
        stage_flags.cluster = any(doptp(uidx) < -tiny_tol);
        stage_flags.second_stage_raw = any(dopto_raw < -tiny_tol);
        stage_flags.roundtrip_band = any(dd_roundtrip_band < -tiny_tol);
        stage_flags.full_dd_check = any(dd_full < -tiny_tol);
        first_stage = first_negative_stage(stage_flags);
        if ~strcmp(first_stage, 'none')
            fprintf(2, '[NEG-DD] %s %s %s m=%d first negative stage: %s\n', ...
                subject_id, stage, bname, m, first_stage);
        end

        diag_struct = struct();
        diag_struct.freq_idx = fidx(:)';
        diag_struct.freq_vals = freqs(fidx)';
        diag_struct.rank_vals = ranks;
        diag_struct.rank_min = min(ranks);
        diag_struct.rank_max = max(ranks);
        diag_struct.n_bins = numel(fidx);
        diag_struct.convo = convo;
        diag_struct.iopto = iopto;
        diag_struct.sopto = sopto;
        diag_struct.cputo = cputo;
        diag_struct.stability_spectral_radius = spec_radius;
        diag_struct.dd_stage_mins = stage_mins;
        diag_struct.dd_stage_maxs = stage_maxs;
        diag_struct.dd_negative_flags = stage_flags;
        diag_struct.dd_first_negative_stage = first_stage;
        diag_struct.dd_clip_counts = struct('full_dd', n_clip_full, 'band_dd', n_clip_band);
        diag_struct.notes = [ ...
            "second_stage_raw and roundtrip_band are band-limited spectral objectives and may be negative; " ...
            "full_dd_check is tracked separately as a positivity diagnostic on full-spectrum DD."];
        if cfg.histo
            diag_struct.ohisto = ohisto;
        end
        results.(sprintf('%s_diag_m%d', bname, m)) = diag_struct;

        dbg = struct();
        dbg.timestamp = datestr(now, 31);
        dbg.subject = subject_id;
        dbg.stage = stage;
        dbg.band = bname;
        dbg.m = m;
        dbg.n_bins = numel(fidx);
        dbg.freq_min = min(freqs(fidx));
        dbg.freq_max = max(freqs(fidx));
        dbg.stage_mins = stage_mins;
        dbg.stage_maxs = stage_maxs;
        dbg.stage_negative = stage_flags;
        dbg.first_negative_stage = first_stage;
        dbg.clip_full = n_clip_full;
        dbg.clip_band = n_clip_band;
        dbg.raw_second_stage_note = 'Band-limited spectral objective can be negative by construction.';
        if isempty(debug_rows)
            debug_rows = dbg;
        else
            debug_rows(end+1,1) = dbg; %#ok<AGROW>
        end
    end
end

results.model = struct('A', state.A, 'V', state.V, 'A_ss', A_ss, 'C_ss', C_ss, 'K_ss', K_ss, 'CAK', CAK);

end

function rec = init_record()
rec = struct( ...
    'subject', '', ...
    'stage', '', ...
    'stage_label', '', ...
    'input_file', '', ...
    'output_file', '', ...
    'success', false, ...
    'message', '', ...
    'runtime_sec', NaN, ...
    'n_regions', NaN, ...
    'n_bands', NaN, ...
    'n_m', NaN, ...
    'nrunsp', NaN, ...
    'nitersp', NaN, ...
    'niterso', NaN ...
);
end

function out = stage_out_dir(cfg, stage)
out = fullfile(cfg.results_root, stage);
if ~exist(out, 'dir')
    mkdir(out);
end
end

function jobs = discover_jobs(cfg)
jobs = {};
for s = 1:numel(cfg.stages)
    stage = cfg.stages{s};
    sm_dir = fullfile(cfg.subject_measure_root, stage);
    files = dir(fullfile(sm_dir, '*_subject_measures.mat'));
    for i = 1:numel(files)
        fn = files(i).name;
        subject_id = erase(fn, '_subject_measures.mat');
        out_file = fullfile(stage_out_dir(cfg, stage), sprintf('%s_ssdi_spectral_domain_results.mat', subject_id));
        if ~cfg.recompute && exist(out_file, 'file') == 2
            continue;
        end
        jobs(end+1,:) = {subject_id, stage}; %#ok<AGROW>
    end
end
if isempty(jobs)
    return;
end
jobs = sortrows(jobs, [2 1]);
end

function state = load_subject_state(subject_id, stage, cfg)
sm_file = fullfile(cfg.subject_measure_root, stage, sprintf('%s_subject_measures.mat', subject_id));
if ~exist(sm_file, 'file')
    error('Missing subject measure file: %s', sm_file);
end

S = load(sm_file, 'subject_measure');
if ~isfield(S, 'subject_measure')
    error('Missing subject_measure struct in %s', sm_file);
end
sm = S.subject_measure;

if ~isfield(sm, 'model') || ~isfield(sm.model, 'A') || ~isfield(sm.model, 'V')
    error('Missing model.A/model.V in %s', sm_file);
end

state = struct();
state.A = double(sm.model.A);
state.V = double(sm.model.V);
state.morder = scalar_default(getfield_safe(sm.model, 'morder'), size(state.A,3)); %#ok<GFLD>
labels = [];
if isfield(sm, 'meta') && isfield(sm.meta, 'labels')
    labels = normalize_labels(sm.meta.labels);
end
if isempty(labels) && isfield(sm, 'meta') && isfield(sm.meta, 'group_names')
    labels = normalize_labels(sm.meta.group_names);
end
if isempty(labels)
    labels = cfg.required_labels;
end
state.labels = labels;
state.input_file = sm_file;
state.source_type = 'subject_measure_model';
end

function labels = normalize_labels(raw)
if isempty(raw)
    labels = {};
    return;
end
if isstring(raw)
    labels = cellstr(raw(:));
elseif ischar(raw)
    labels = {raw};
elseif iscell(raw)
    labels = cell(size(raw));
    for i = 1:numel(raw)
        if isstring(raw{i})
            labels{i} = char(raw{i});
        elseif ischar(raw{i})
            labels{i} = raw{i};
        else
            labels{i} = strtrim(char(string(raw{i})));
        end
    end
    labels = labels(:)';
else
    labels = cellstr(string(raw));
end
labels = labels(:)';
end

function x = scalar_default(v, fallback)
if isempty(v)
    x = fallback;
    return;
end
if isnumeric(v)
    x = double(v(1));
else
    x = fallback;
end
end

function v = getfield_safe(s, field)
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = [];
end
end

function band_info = build_band_info(freqs, band_defs)
n = numel(band_defs);
band_info = repmat(struct('name','','lo',0,'hi',0,'indices',[],'freq_values',[],'n_bins',0,'is_empty',true), n, 1);
for i = 1:n
    lo = band_defs(i).lo;
    hi = band_defs(i).hi;
    mask = freqs >= lo & freqs <= hi;
    idx = find(mask);
    band_info(i).name = band_defs(i).name;
    band_info(i).lo = lo;
    band_info(i).hi = hi;
    band_info(i).indices = idx(:)';
    band_info(i).freq_values = freqs(idx)';
    band_info(i).n_bins = numel(idx);
    band_info(i).is_empty = isempty(idx);
end
end

function write_frequency_diagnostics(cfg, freqs, band_info)
out_mat = fullfile(cfg.diag_dir, 'frequency_vector.mat');
out_txt = fullfile(cfg.diag_dir, 'frequency_band_summary.txt');

save(out_mat, 'freqs', 'band_info');

fid = fopen(out_txt, 'w');
assert(fid > 0, 'Could not write %s', out_txt);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Frequency diagnostics for SomnoDynamics stage SSDI anatomical7\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'fs = %.6f Hz\n', cfg.fs);
fprintf(fid, 'fres = %d\n', cfg.fres);
fprintf(fid, 'n_freq = %d\n', numel(freqs));
fprintf(fid, 'range = [%.6f, %.6f] Hz\n', min(freqs), max(freqs));
if numel(freqs) > 1
    fprintf(fid, 'delta_f = %.6f Hz\n\n', freqs(2)-freqs(1));
else
    fprintf(fid, 'delta_f = n/a\n\n');
end

for i = 1:numel(band_info)
    bi = band_info(i);
    fprintf(fid, '%s: [%.3f, %.3f] (inclusive)\n', bi.name, bi.lo, bi.hi);
    fprintf(fid, '  n_bins: %d\n', bi.n_bins);
    if bi.n_bins > 0
        fprintf(fid, '  idx: %s\n', mat2str(bi.indices));
        fprintf(fid, '  freqs: %s\n\n', mat2str(bi.freq_values, 6));
    else
        fprintf(fid, '  idx: []\n  freqs: []\n\n');
    end
end
end

function maybe_start_parpool(cfg)
want_parallel = logical(cfg.ppp) || logical(cfg.ppo);
if ~want_parallel
    return;
end
p = gcp('nocreate');
if ~isempty(p)
    fprintf('Using existing parpool with %d workers.\n', p.NumWorkers);
    return;
end

ncores = feature('numcores');
workers = min([cfg.max_parpool_workers, ncores]);
workers = max(workers, 1);

try
    parpool('local', workers);
    fprintf('Started parpool with %d workers.\n', workers);
catch ME
    warning('Could not start parpool (%s). Falling back to serial mode in SSDI routines.', ME.message);
end
end

function append_status(cfg, subject_id, stage, ok, msg)
status_file = fullfile(cfg.logs_dir, 'ssdi_anatomical7_stage_status.csv');
need_header = ~exist(status_file, 'file');
fid = fopen(status_file, 'a');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
if need_header
    fprintf(fid, 'timestamp,subject,stage,stage_label,success,message\n');
end
msg = strrep(msg, '"', '''');
msg = strrep(msg, sprintf('\n'), ' ');
fprintf(fid, '%s,%s,%s,%s,%d,"%s"\n', datestr(now, 31), subject_id, stage, cfg.stage_display.(stage), ok, msg);
end

function [dd, dd_raw] = trfun_dd_batch(Lset, H)
nr = size(Lset,3);
dd_raw = nan(1, nr);
for k = 1:nr
    dd_raw(k) = trfun2dd(Lset(:,:,k), H);
end
dd = dd_raw;
end

function dd = canonical_band_dd_batch(Lset, H, freqs, fidx)
nr = size(Lset, 3);
dd = nan(1, nr);
for k = 1:nr
    [d, ~, ~, ok] = canonical_band_objective_and_grad(Lset(:,:,k), H, freqs, fidx);
    if ok && isfinite(d)
        dd(k) = d;
    end
end
end

function [x, nclip] = clip_tiny_negative(x, tol)
if nargin < 2
    tol = 1e-10;
end
mask_tiny = isfinite(x) & x < 0 & x >= -tol;
nclip = nnz(mask_tiny);
x(mask_tiny) = 0;
end

function stage = first_negative_stage(flags)
order = {'preopt','cluster','second_stage_raw','roundtrip_band','full_dd_check'};
stage = 'none';
for i = 1:numel(order)
    key = order{i};
    if isfield(flags, key) && logical(flags.(key))
        stage = key;
        return;
    end
end
end

function append_dd_debug(cfg, debug_rows)
if isempty(debug_rows)
    return;
end

mat_file = fullfile(cfg.diag_dir, 'dd_negativity_debug.mat');
if exist(mat_file, 'file')
    S = load(mat_file, 'debug_entries');
    if isfield(S, 'debug_entries')
        debug_entries = S.debug_entries; %#ok<NASGU>
    else
        debug_entries = struct([]);
    end
else
    debug_entries = struct([]);
end

if isempty(debug_entries)
    debug_entries = debug_rows(:);
else
    debug_entries = [debug_entries(:); debug_rows(:)]; %#ok<AGROW>
end
save(mat_file, 'debug_entries');

write_dd_debug_summary(cfg, debug_entries);
end

function write_dd_debug_summary(cfg, entries)
txt_file = fullfile(cfg.diag_dir, 'dd_negativity_debug_summary.txt');
fid = fopen(txt_file, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'DD negativity debug summary (SomnoDynamics stage SSDI)\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Rows: %d\n\n', numel(entries));

fprintf(fid, 'Interpretation\n');
fprintf(fid, '- second_stage_raw / roundtrip_band refer to band-limited spectral objective values.\n');
fprintf(fid, '- full_dd_check is the true DD evaluated on full-spectrum transfer function.\n');
fprintf(fid, '- full_dd_check should be non-negative except tiny floating-point error near zero.\n\n');

n_neg_pre = 0;
n_neg_clu = 0;
n_neg_raw = 0;
n_neg_round = 0;
n_neg_full = 0;

for i = 1:numel(entries)
    e = entries(i);
    n_neg_pre = n_neg_pre + double(logical(stage_field(e, 'stage_negative', 'preopt', false)));
    n_neg_clu = n_neg_clu + double(logical(stage_field(e, 'stage_negative', 'cluster', false)));
    n_neg_raw = n_neg_raw + double(logical(stage_field(e, 'stage_negative', 'second_stage_raw', false)));
    n_neg_round = n_neg_round + double(logical(stage_field(e, 'stage_negative', 'roundtrip_band', false)));
    n_neg_full = n_neg_full + double(logical(stage_field(e, 'stage_negative', 'full_dd_check', false)));
end

fprintf(fid, 'Negative-flag counts by stage\n');
fprintf(fid, '- preopt: %d\n', n_neg_pre);
fprintf(fid, '- cluster: %d\n', n_neg_clu);
fprintf(fid, '- second_stage_raw: %d\n', n_neg_raw);
fprintf(fid, '- roundtrip_band: %d\n', n_neg_round);
fprintf(fid, '- full_dd_check: %d\n\n', n_neg_full);

fprintf(fid, 'Detailed rows\n');
for i = 1:numel(entries)
    e = entries(i);
    min_pre = stage_field(e, 'stage_mins', 'preopt', nan);
    min_clu = stage_field(e, 'stage_mins', 'cluster', nan);
    min_raw = stage_field(e, 'stage_mins', 'second_stage_raw', nan);
    min_round = stage_field(e, 'stage_mins', 'roundtrip_band', nan);
    min_full = stage_field(e, 'stage_mins', 'full_dd_check', nan);
    max_pre = stage_field(e, 'stage_maxs', 'preopt', nan);
    max_clu = stage_field(e, 'stage_maxs', 'cluster', nan);
    max_raw = stage_field(e, 'stage_maxs', 'second_stage_raw', nan);
    max_round = stage_field(e, 'stage_maxs', 'roundtrip_band', nan);
    max_full = stage_field(e, 'stage_maxs', 'full_dd_check', nan);
    neg_pre = stage_field(e, 'stage_negative', 'preopt', 0);
    neg_clu = stage_field(e, 'stage_negative', 'cluster', 0);
    neg_raw = stage_field(e, 'stage_negative', 'second_stage_raw', 0);
    neg_round = stage_field(e, 'stage_negative', 'roundtrip_band', 0);
    neg_full = stage_field(e, 'stage_negative', 'full_dd_check', 0);
    clip_full = scalar_field(e, 'clip_full', 0);
    clip_band = scalar_field(e, 'clip_band', 0);
    fprintf(fid, '%s | %s %s %s | m=%d | bins=%d [%.3f-%.3f] | first_negative=%s\n', ...
        e.timestamp, e.subject, e.stage, e.band, e.m, e.n_bins, e.freq_min, e.freq_max, e.first_negative_stage);
    fprintf(fid, '  mins: pre=%.6g cluster=%.6g raw=%.6g round=%.6g fullcheck=%.6g\n', ...
        min_pre, min_clu, min_raw, min_round, min_full);
    fprintf(fid, '  maxs: pre=%.6g cluster=%.6g raw=%.6g round=%.6g fullcheck=%.6g\n', ...
        max_pre, max_clu, max_raw, max_round, max_full);
    fprintf(fid, '  neg: pre=%d cluster=%d raw=%d round=%d fullcheck=%d\n', ...
        neg_pre, neg_clu, neg_raw, neg_round, neg_full);
    fprintf(fid, '  clips: full=%d band=%d\n\n', clip_full, clip_band);
end
end

function v = stage_field(s, stage_name, field_name, fallback)
v = fallback;
if ~isstruct(s) || ~isfield(s, stage_name)
    return;
end
stage = s.(stage_name);
if ~isstruct(stage) || ~isfield(stage, field_name)
    return;
end
v = stage.(field_name);
if isempty(v)
    v = fallback;
end
end

function v = scalar_field(s, field_name, fallback)
v = fallback;
if ~isstruct(s) || ~isfield(s, field_name)
    return;
end
val = s.(field_name);
if isempty(val)
    return;
end
v = val;
end
