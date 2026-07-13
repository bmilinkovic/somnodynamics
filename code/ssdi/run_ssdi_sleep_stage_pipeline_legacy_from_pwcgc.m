function run_ssdi_sleep_stage_pipeline_legacy_from_pwcgc(varargin)
%RUN_SSDI_SLEEP_STAGE_PIPELINE Stage-based SSDI from GC pwcgc parameter files.
%
% Input:
%   results/gc/pwcgc/pwcgc_params_EPCTLxx_<stage>.mat
%
% Output:
%   results/ssdi/time/*_ssdi_time_domain_results.mat
%   results/ssdi/spectral/*_ssdi_spectral_domain_results.mat

p = inputParser;
addParameter(p, 'pwcgc_dir', '');
addParameter(p, 'out_root', '');
addParameter(p, 'fres', 128);
addParameter(p, 'fs', 128);
addParameter(p, 'nrunsp', 30);
addParameter(p, 'niters', 3000);
addParameter(p, 'm_values', []);
addParameter(p, 'bands', { ...
    'delta', 1, 4; ...
    'theta', 5, 8; ...
    'alpha', 9, 12; ...
    'sigma', 12, 16; ...
    'beta', 13, 30; ...
    'gamma', 31, 44 ...
});
parse(p, varargin{:});
cfg = p.Results;

this_file = mfilename('fullpath');
code_dir = fileparts(this_file);
repo_root = fileparts(fileparts(code_dir));
legacy_root = fullfile(repo_root, 'code', 'legacy');
run_all_root = fullfile(repo_root, 'code', 'run_all');
addpath(legacy_root);
addpath(run_all_root);
setup_matlab_toolboxes(repo_root);

if isempty(cfg.pwcgc_dir)
    cfg.pwcgc_dir = fullfile(sd_external_results_root(repo_root), 'gc', 'pwcgc');
end
if isempty(cfg.out_root)
    cfg.out_root = fullfile(sd_external_results_root(repo_root), 'ssdi');
end

out_time = fullfile(cfg.out_root, 'time');
out_spec = fullfile(cfg.out_root, 'spectral');
out_logs = fullfile(cfg.out_root, 'logs');
mkdir_if_missing(cfg.out_root);
mkdir_if_missing(out_time);
mkdir_if_missing(out_spec);
mkdir_if_missing(out_logs);

if exist('opt_gd_dds_mruns', 'file') ~= 2
    error(['SSDI optimisation function opt_gd_dds_mruns not found. ' ...
        'Add the SSDI toolbox to the MATLAB path before running.']);
end
if exist('var_to_ss', 'file') ~= 2
    error('MVGC/SSDI dependencies missing on path (var_to_ss unresolved).');
end

files = dir(fullfile(cfg.pwcgc_dir, 'pwcgc_params_EPCTL*_*.mat'));
if isempty(files)
    error('No pwcgc files found in %s', cfg.pwcgc_dir);
end

freqs = sfreqs(cfg.fres, cfg.fs);

fprintf('\n[SomnoDynamics SSDI] input pwcgc: %s\n', cfg.pwcgc_dir);
fprintf('[SomnoDynamics SSDI] output root: %s\n', cfg.out_root);
fprintf('[SomnoDynamics SSDI] files: %d\n', numel(files));

for i = 1:numel(files)
    fname = files(i).name;
    fpath = fullfile(files(i).folder, fname);

    parsed = regexp(fname, '^pwcgc_params_(EPCTL\d{2})_([A-Za-z0-9_]+)\.mat$', 'tokens', 'once');
    if isempty(parsed)
        warning('Skipping unmatched file: %s', fname);
        continue;
    end
    subject_id = parsed{1};
    stage = parsed{2};

    S = load(fpath);
    if ~isfield(S, 'CAK') || ~isfield(S, 'V')
        warning('Skipping %s (missing CAK/V)', fname);
        continue;
    end

    Ass = double(S.CAK);
    V = double(S.V);
    [A_ss, C_ss, K_ss] = var_to_ss(Ass);
    H = ss2trfun(A_ss, C_ss, K_ss, cfg.fres);

    n = size(Ass, 1);
    if isempty(cfg.m_values)
        m_values = 2:max(2, n-1);
    else
        m_values = cfg.m_values;
    end

    labels = {};
    if isfield(S, 'labels')
        labels = normalize_labels(S.labels);
    end
    if isempty(labels)
        labels = default_labels(n);
    end

    fprintf('[%d/%d] %s %s | n=%d\n', i, numel(files), subject_id, stage, n);

    res_time = struct();
    res_time.meta = struct('subject_id', subject_id, 'stage', stage, 'labels', {labels}, ...
        'fres', cfg.fres, 'fs', cfg.fs, 'n_regions', n, 'source_file', fpath, ...
        'note', 'Broadband objective over full spectral transfer function (sleep-stage window).');

    res_spec = struct();
    res_spec.meta = struct('subject_id', subject_id, 'stage', stage, 'labels', {labels}, ...
        'fres', cfg.fres, 'fs', cfg.fs, 'n_regions', n, 'source_file', fpath);

    for m = m_values
        if m >= n
            continue;
        end

        seed = 1000 + i*100 + m;
        rs = rng_seed(seed);
        L0 = rand_orthonormal(n, m, cfg.nrunsp);
        rng_restore(rs);

        [dopto_t, Lo_t, conv_t, iopt_t, sopt_t, cput_t, ohist_t] = opt_gd_dds_mruns( ...
            H, L0, cfg.niters, 2, 0.1, 2, 1e-10, true, true);
        Lopto_t = itransform_subspace(Lo_t, V);
        res_time.(sprintf('Lopto_m%d', m)) = Lopto_t;
        res_time.(sprintf('dopto_m%d', m)) = dopto_t;
        res_time.(sprintf('gopto_m%d', m)) = gmetrics(Lopto_t);
        res_time.(sprintf('opt_m%d', m)) = struct( ...
            'conv', conv_t, 'iopt', iopt_t, 'sopt', sopt_t, 'cput', cput_t, 'ohist', ohist_t);
        [nw_t, beta_t] = compute_node_metrics(Lopto_t);
        res_time.(sprintf('node_weights_m%d', m)) = nw_t;
        res_time.(sprintf('beta_m%d', m)) = beta_t;

        for b = 1:size(cfg.bands,1)
            bname = cfg.bands{b,1};
            lo = cfg.bands{b,2};
            hi = cfg.bands{b,3};
            idx = find(freqs >= lo & freqs <= hi);
            if isempty(idx)
                continue;
            end

            H_band = H(:,:,idx);
            rsb = rng_seed(seed + 100*b);
            L0b = rand_orthonormal(n, m, cfg.nrunsp);
            rng_restore(rsb);

            [dopto, Lo, conv, iopt, sopt, cput, ohist] = opt_gd_dds_mruns( ...
                H_band, L0b, cfg.niters, 2, 0.1, 2, 1e-10, true, true);

            Lopto = itransform_subspace(Lo, V);
            keyL = sprintf('%s_Lopto_m', bname);
            if ~isfield(res_spec, keyL)
                res_spec.(keyL) = cell(1, max(m_values)-1);
            end
            temp = res_spec.(keyL);
            temp{m-1} = Lopto;
            res_spec.(keyL) = temp;

            res_spec.(sprintf('%s_Lopto_m%d', bname, m)) = Lopto;
            res_spec.(sprintf('%s_dopto_m%d', bname, m)) = dopto;
            res_spec.(sprintf('%s_gopto_m%d', bname, m)) = gmetrics(Lopto);
            res_spec.(sprintf('%s_opt_m%d', bname, m)) = struct( ...
                'conv', conv, 'iopt', iopt, 'sopt', sopt, 'cput', cput, 'ohist', ohist, ...
                'freq_idx', idx(:)', 'freq_lo', lo, 'freq_hi', hi);

            [nw, beta] = compute_node_metrics(Lopto);
            res_spec.(sprintf('%s_node_weights_m%d', bname, m)) = nw;
            res_spec.(sprintf('%s_beta_m%d', bname, m)) = beta;
        end
    end

    out_time_file = fullfile(out_time, sprintf('%s_%s_ssdi_time_domain_results.mat', subject_id, stage));
    out_spec_file = fullfile(out_spec, sprintf('%s_%s_ssdi_spectral_domain_results.mat', subject_id, stage));

    results = res_time; %#ok<NASGU>
    save(out_time_file, 'results', '-v7.3');
    results = res_spec; %#ok<NASGU>
    save(out_spec_file, 'results', '-v7.3');

    fprintf('  saved time    -> %s\n', out_time_file);
    fprintf('  saved spectral-> %s\n', out_spec_file);
end

fprintf('[SomnoDynamics SSDI] complete.\n');
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
    labels = {raw};
elseif iscell(raw)
    labels = cellfun(@(x) char(string(x)), raw, 'UniformOutput', false);
else
    labels = cellstr(string(raw));
end
labels = labels(:)';
end

function labels = default_labels(n)
labels = cell(1,n);
for i = 1:n
    labels{i} = sprintf('Region_%d', i);
end
end

function mkdir_if_missing(path)
if exist(path, 'dir') ~= 7
    mkdir(path);
end
end
