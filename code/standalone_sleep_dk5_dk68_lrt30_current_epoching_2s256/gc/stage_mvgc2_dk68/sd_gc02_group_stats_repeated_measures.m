function sd_gc02_group_stats_repeated_measures(cfg)
% Stage 2: repeated-measures GC/MI statistics across sleep stages.

if nargin < 1
    cfg = sd_gc_config_anatomical7();
end

out_mat = fullfile(cfg.paths.stats, 'repeated_measures_stats.mat');
out_omnibus_csv = fullfile(cfg.paths.stats, 'omnibus_long.csv');
out_pairwise_csv = fullfile(cfg.paths.stats, 'pairwise_long.csv');
out_summary_txt = fullfile(cfg.paths.reports, 'repeated_measures_summary.txt');

if ~cfg.recompute_group_stats && exist(out_mat, 'file') == 2
    fprintf('\n[Stage 2] Existing repeated-measures stats found, skipping: %s\n', out_mat);
    return;
end

fprintf('  Loading stage-wise subject measures...\n');
D = load_stage_measure_stacks(cfg);

fprintf('  Building repeated-measures tests (subjects=%d, stages=%d)\n', D.nsub, D.nstage);

[omni_rows, pair_rows] = build_all_stats(D, cfg);

T_omni = struct2table(omni_rows, 'AsArray', true);
T_pair = struct2table(pair_rows, 'AsArray', true);

T_omni = apply_family_fdr(T_omni, cfg.fdr_primary, cfg.fdr_sensitivity, cfg.alpha_group);
T_pair = apply_family_fdr(T_pair, cfg.fdr_primary, cfg.fdr_sensitivity, cfg.alpha_group);

writetable(T_omni, out_omnibus_csv);
writetable(T_pair, out_pairwise_csv);

stats = struct();
stats.meta = struct(...
    'subjects', {D.subjects}, ...
    'stages', {cfg.stages}, ...
    'stage_display', cfg.stage_display, ...
    'labels', {D.labels}, ...
    'freqs', D.freqs, ...
    'band_info', D.band_info, ...
    'alpha_group', cfg.alpha_group, ...
    'fdr_primary', cfg.fdr_primary, ...
    'fdr_sensitivity', cfg.fdr_sensitivity, ...
    'directed_convention', cfg.directed_convention ...
);
stats.omnibus = T_omni;
stats.pairwise = T_pair;
save(out_mat, 'stats', '-v7.3');

write_summary_report(out_summary_txt, D, T_omni, T_pair, cfg);

fprintf('  Saved repeated-measures stats MAT: %s\n', out_mat);
fprintf('  Saved omnibus table            : %s\n', out_omnibus_csv);
fprintf('  Saved pairwise table           : %s\n', out_pairwise_csv);
fprintf('  Saved summary report           : %s\n', out_summary_txt);

end

function D = load_stage_measure_stacks(cfg)
% Discover subjects present in all required stages.
subject_sets = cell(numel(cfg.stages), 1);
for s = 1:numel(cfg.stages)
    stage = cfg.stages{s};
    files = dir(fullfile(cfg.paths.subject_measures, stage, '*_subject_measures.mat'));
    ids = erase({files.name}, '_subject_measures.mat');
    subject_sets{s} = ids(:);
end

subjects = subject_sets{1};
for s = 2:numel(subject_sets)
    subjects = intersect(subjects, subject_sets{s}, 'stable');
end
if isempty(subjects)
    error('No subjects found with all required stages in %s', cfg.paths.subject_measures);
end

nsub = numel(subjects);
nstage = numel(cfg.stages);

% Reference metadata.
ref_file = fullfile(cfg.paths.subject_measures, cfg.stages{1}, sprintf('%s_subject_measures.mat', subjects{1}));
R0 = load(ref_file, 'subject_measure');
M0 = R0.subject_measure;
labels = M0.meta.group_names(:);
nreg = numel(labels);
freqs = M0.meta.freqs(:)';
band_info = M0.meta.band_info;
h = numel(freqs);

% Preallocate stacks.
dir_src_in_t = nan(nsub, nstage, nreg);
dir_src_out_t = nan(nsub, nstage, nreg);
dir_roi_t = nan(nsub, nstage, nreg, nreg);

undir_src_t = nan(nsub, nstage, nreg);
undir_roi_t = nan(nsub, nstage, nreg, nreg);

dir_src_in_f = nan(nsub, nstage, nreg, h);
dir_src_out_f = nan(nsub, nstage, nreg, h);
dir_roi_f = nan(nsub, nstage, nreg, nreg, h);

undir_src_f = nan(nsub, nstage, nreg, h);
undir_roi_f = nan(nsub, nstage, nreg, nreg, h);

for si = 1:nsub
    subj = subjects{si};
    for st = 1:nstage
        stage = cfg.stages{st};
        f = fullfile(cfg.paths.subject_measures, stage, sprintf('%s_subject_measures.mat', subj));
        S = load(f, 'subject_measure');
        M = S.subject_measure;

        if ~isequal(labels, M.meta.group_names(:))
            error('Label mismatch for subject %s stage %s', subj, stage);
        end

        dir_src_in_t(si, st, :) = reshape(M.measures.directed.time.src_in, 1, 1, nreg);
        dir_src_out_t(si, st, :) = reshape(M.measures.directed.time.src_out, 1, 1, nreg);
        dir_roi_t(si, st, :, :) = reshape(M.measures.directed.time.roi, 1, 1, nreg, nreg);

        undir_src_t(si, st, :) = reshape(M.measures.undirected.time.src, 1, 1, nreg);
        undir_roi_t(si, st, :, :) = reshape(M.measures.undirected.time.roi, 1, 1, nreg, nreg);

        dir_src_in_f(si, st, :, :) = reshape(M.measures.directed.spectral.src_in, 1, 1, nreg, h);
        dir_src_out_f(si, st, :, :) = reshape(M.measures.directed.spectral.src_out, 1, 1, nreg, h);
        dir_roi_f(si, st, :, :, :) = reshape(M.measures.directed.spectral.roi, 1, 1, nreg, nreg, h);

        undir_src_f(si, st, :, :) = reshape(M.measures.undirected.spectral.src, 1, 1, nreg, h);
        undir_roi_f(si, st, :, :, :) = reshape(M.measures.undirected.spectral.roi, 1, 1, nreg, nreg, h);
    end
end

D = struct();
D.subjects = subjects;
D.nsub = nsub;
D.nstage = nstage;
D.labels = labels;
D.nreg = nreg;
D.freqs = freqs;
D.band_info = band_info;

D.time = struct();
D.time.directed = struct('src_in', dir_src_in_t, 'src_out', dir_src_out_t, 'roi', dir_roi_t);
D.time.undirected = struct('src', undir_src_t, 'roi', undir_roi_t);

D.spectral = struct();
D.spectral.directed = struct('src_in', dir_src_in_f, 'src_out', dir_src_out_f, 'roi', dir_roi_f);
D.spectral.undirected = struct('src', undir_src_f, 'roi', undir_roi_f);
end

function [omni_rows, pair_rows] = build_all_stats(D, cfg)
omni_rows = repmat(init_omni_row(), 0, 1);
pair_rows = repmat(init_pair_row(), 0, 1);

% Time-domain families.
[omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, ...
    squeeze(D.time.directed.src_in), D.labels, cfg, 'time', 'all', 'directed_src_in', 'nodewise', 'in');
[omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, ...
    squeeze(D.time.directed.src_out), D.labels, cfg, 'time', 'all', 'directed_src_out', 'nodewise', 'out');
[omni_rows, pair_rows] = append_matrix_family(omni_rows, pair_rows, ...
    D.time.directed.roi, D.labels, cfg, 'time', 'all', 'directed_roi', 'edgewise', 'directed');

[omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, ...
    squeeze(D.time.undirected.src), D.labels, cfg, 'time', 'all', 'undirected_src', 'nodewise', 'undirected');
[omni_rows, pair_rows] = append_matrix_family(omni_rows, pair_rows, ...
    D.time.undirected.roi, D.labels, cfg, 'time', 'all', 'undirected_roi', 'edgewise', 'undirected');

% Spectral families collapsed by configured bands.
for b = 1:numel(D.band_info)
    bname = D.band_info(b).name;
    idx = D.band_info(b).mask;

    x = mean(D.spectral.directed.src_in(:, :, :, idx), 4, 'omitnan');
    [omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, squeeze(x), D.labels, ...
        cfg, 'spectral', bname, 'directed_src_in', 'nodewise', 'in');

    x = mean(D.spectral.directed.src_out(:, :, :, idx), 4, 'omitnan');
    [omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, squeeze(x), D.labels, ...
        cfg, 'spectral', bname, 'directed_src_out', 'nodewise', 'out');

    x = mean(D.spectral.directed.roi(:, :, :, :, idx), 5, 'omitnan');
    [omni_rows, pair_rows] = append_matrix_family(omni_rows, pair_rows, x, D.labels, ...
        cfg, 'spectral', bname, 'directed_roi', 'edgewise', 'directed');

    x = mean(D.spectral.undirected.src(:, :, :, idx), 4, 'omitnan');
    [omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, squeeze(x), D.labels, ...
        cfg, 'spectral', bname, 'undirected_src', 'nodewise', 'undirected');

    x = mean(D.spectral.undirected.roi(:, :, :, :, idx), 5, 'omitnan');
    [omni_rows, pair_rows] = append_matrix_family(omni_rows, pair_rows, x, D.labels, ...
        cfg, 'spectral', bname, 'undirected_roi', 'edgewise', 'undirected');
end

end

function [omni_rows, pair_rows] = append_vector_family(omni_rows, pair_rows, X, labels, cfg, domain, band, measure_family, feature_type, direction)
% X: nsub x nstage x nfeat
if ndims(X) ~= 3
    error('Vector family expects 3D array [nsub x nstage x nfeat].');
end

nfeat = size(X, 3);
for f = 1:nfeat
    feat = squeeze(X(:, :, f));
    target = labels{f};
    source = '';

    [orow, prows] = compute_feature_stats(feat, cfg, domain, band, measure_family, feature_type, direction, target, source);

    omni_rows(end+1, 1) = orow; %#ok<AGROW>
    pair_rows = [pair_rows; prows]; %#ok<AGROW>
end
end

function [omni_rows, pair_rows] = append_matrix_family(omni_rows, pair_rows, X, labels, cfg, domain, band, measure_family, feature_type, direction)
% X: nsub x nstage x nto x nfrom
if ndims(X) ~= 4
    error('Matrix family expects 4D array [nsub x nstage x nto x nfrom].');
end

nreg = numel(labels);
for to = 1:nreg
    for from = 1:nreg
        if to == from
            continue;
        end
        feat = squeeze(X(:, :, to, from));
        target = labels{to};
        source = labels{from};

        [orow, prows] = compute_feature_stats(feat, cfg, domain, band, measure_family, feature_type, direction, target, source);

        omni_rows(end+1, 1) = orow; %#ok<AGROW>
        pair_rows = [pair_rows; prows]; %#ok<AGROW>
    end
end
end

function [orow, prows] = compute_feature_stats(feat, cfg, domain, band, measure_family, feature_type, direction, target, source)
% feat: nsub x nstage
if size(feat, 2) ~= numel(cfg.stages)
    error('Feature matrix second dimension must match configured stages.');
end

valid = all(isfinite(feat), 2);
X = feat(valid, :);

orow = init_omni_row();
orow.domain = domain;
orow.band = band;
orow.measure_family = measure_family;
orow.feature_type = feature_type;
orow.direction = direction;
orow.target = target;
orow.source = source;
orow.n_subjects_total = size(feat, 1);
orow.n_subjects_effective = size(X, 1);
orow.correction_family = sprintf('%s|%s|%s', domain, band, measure_family);

if size(X, 1) < 3
    orow.friedman_chi2 = NaN;
    orow.p = NaN;
else
    [p_friedman, chi2] = friedman_stat(X);
    orow.friedman_chi2 = chi2;
    orow.p = p_friedman;
end

pairs = nchoosek(1:numel(cfg.stages), 2);
prows = repmat(init_pair_row(), size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    a = pairs(i, 1);
    b = pairs(i, 2);

    pr = init_pair_row();
    pr.domain = domain;
    pr.band = band;
    pr.measure_family = measure_family;
    pr.feature_type = feature_type;
    pr.direction = direction;
    pr.target = target;
    pr.source = source;

    pr.stage_a = cfg.stages{a};
    pr.stage_b = cfg.stages{b};
    pr.stage_a_label = cfg.stage_display.(pr.stage_a);
    pr.stage_b_label = cfg.stage_display.(pr.stage_b);

    pr.correction_family = sprintf('%s|%s|%s|%s_vs_%s', domain, band, measure_family, pr.stage_a, pr.stage_b);

    if size(X, 1) < 3
        pr.n_subjects_effective = size(X, 1);
        pr.effect_mean = NaN;
        pr.effect_median = NaN;
        pr.rank_corr = NaN;
        pr.W = NaN;
        pr.z = NaN;
        pr.p = NaN;
    else
        xa = X(:, a);
        xb = X(:, b);
        st = paired_signedrank_stats(xa, xb);
        pr.n_subjects_effective = st.n_eff;
        pr.effect_mean = st.effect_mean;
        pr.effect_median = st.effect_median;
        pr.rank_corr = st.rank_corr;
        pr.W = st.W;
        pr.z = st.z;
        pr.p = st.p;
    end

    prows(i, 1) = pr;
end
end

function st = paired_signedrank_stats(a, b)
mask = isfinite(a) & isfinite(b);
a = a(mask);
b = b(mask);
diffv = b - a;

st = struct('n_eff', 0, 'effect_mean', NaN, 'effect_median', NaN, 'rank_corr', NaN, 'W', NaN, 'z', NaN, 'p', NaN);
if numel(diffv) < 3
    st.n_eff = numel(diffv);
    return;
end

st.n_eff = sum(diffv ~= 0);
st.effect_mean = mean(diffv, 'omitnan');
st.effect_median = median(diffv, 'omitnan');

try
    [p, ~, stats] = signrank(b, a, 'method', 'approximate');
catch
    [p, ~, stats] = signrank(b, a);
end

st.p = p;
if isfield(stats, 'signedrank')
    st.W = stats.signedrank;
end
if isfield(stats, 'zval')
    st.z = stats.zval;
elseif isfield(stats, 'z')
    st.z = stats.z;
end

if st.n_eff > 0 && isfinite(st.W)
    T = st.n_eff * (st.n_eff + 1) / 2;
    if T > 0
        st.rank_corr = (2 * st.W / T) - 1;
    end
end
end

function [p, chi2] = friedman_stat(X)
chi2 = NaN;
p = NaN;
try
    [p, tbl] = friedman(X, 1, 'off');
    if iscell(tbl) && size(tbl, 1) >= 2 && size(tbl, 2) >= 5
        v = tbl{2, 5};
        if isnumeric(v) && isscalar(v)
            chi2 = v;
        end
    end
catch
    p = NaN;
    chi2 = NaN;
end
end

function T = apply_family_fdr(T, method_primary, method_sensitivity, alpha)
T.q_primary = nan(height(T), 1);
T.sig_primary = false(height(T), 1);
T.q_sensitivity = nan(height(T), 1);
T.sig_sensitivity = false(height(T), 1);
T.n_tests_family = zeros(height(T), 1);

families = unique(T.correction_family, 'stable');
for i = 1:numel(families)
    fam = families{i};
    idx = strcmp(T.correction_family, fam);
    pvals = T.p(idx);

    q_primary = fdr_adjust(pvals, method_primary);
    q_sens = fdr_adjust(pvals, method_sensitivity);

    T.q_primary(idx) = q_primary;
    T.sig_primary(idx) = q_primary <= alpha;
    T.q_sensitivity(idx) = q_sens;
    T.sig_sensitivity(idx) = q_sens <= alpha;
    T.n_tests_family(idx) = sum(isfinite(pvals));
end

T.correction_method_primary = repmat({method_primary}, height(T), 1);
T.correction_method_sensitivity = repmat({method_sensitivity}, height(T), 1);
end

function q = fdr_adjust(pvals, method)
q = nan(size(pvals));
idx = find(isfinite(pvals));
if isempty(idx)
    return;
end

p = pvals(idx);
[m_sorted, ord] = sort(p(:));
m = numel(m_sorted);

if strcmpi(method, 'BY')
    c_m = sum(1 ./ (1:m));
else
    c_m = 1;
end

ranks = (1:m)';
q_sorted = (m_sorted .* m .* c_m) ./ ranks;
q_sorted = flipud(cummin(flipud(q_sorted)));
q_sorted = min(q_sorted, 1);

q_local = nan(size(p));
q_local(ord) = q_sorted;
q(idx) = q_local;
end

function write_summary_report(out_file, D, T_omni, T_pair, cfg)
fid = fopen(out_file, 'w');
if fid < 0
    warning('Could not open summary report file for writing: %s', out_file);
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'SomnoDynamics stage-based MVGC2 repeated-measures report\n');
fprintf(fid, '=======================================================\n\n');
fprintf(fid, 'Timestamp: %s\n', datestr(now, 31));
fprintf(fid, 'Subjects with all stages: %d\n', D.nsub);
fprintf(fid, 'Stages: %s\n', strjoin(cfg.stages, ', '));
fprintf(fid, 'Stage labels: W_=Wake, N1=N1, N2=N2, N3=N3, R_=REM\n');
fprintf(fid, 'Primary FDR: %s\n', cfg.fdr_primary);
fprintf(fid, 'Sensitivity FDR: %s\n\n', cfg.fdr_sensitivity);

fprintf(fid, 'Omnibus rows: %d\n', height(T_omni));
fprintf(fid, 'Pairwise rows: %d\n\n', height(T_pair));

if ~isempty(T_omni)
    min_q = min(T_omni.q_primary, [], 'omitnan');
    n_sig = sum(T_omni.sig_primary > 0);
    fprintf(fid, 'Omnibus min q_primary: %.6g\n', min_q);
    fprintf(fid, 'Omnibus significant rows (q<=%.3f): %d\n\n', cfg.alpha_group, n_sig);
end

if ~isempty(T_pair)
    min_q = min(T_pair.q_primary, [], 'omitnan');
    n_sig = sum(T_pair.sig_primary > 0);
    fprintf(fid, 'Pairwise min q_primary: %.6g\n', min_q);
    fprintf(fid, 'Pairwise significant rows (q<=%.3f): %d\n\n', cfg.alpha_group, n_sig);
end

fprintf(fid, 'Top omnibus hits (up to 20 rows by q_primary then p):\n');
if isempty(T_omni)
    fprintf(fid, '  none\n');
else
    Ts = sortrows(T_omni, {'q_primary','p'}, {'ascend','ascend'});
    nshow = min(20, height(Ts));
    for i = 1:nshow
        r = Ts(i, :);
        fprintf(fid, '  %s | %s | %s | %s <- %s | p=%.4g q=%.4g\n', ...
            r.domain{1}, r.band{1}, r.measure_family{1}, r.target{1}, r.source{1}, r.p, r.q_primary);
    end
end

end

function r = init_omni_row()
r = struct(...
    'domain', '', ...
    'band', '', ...
    'measure_family', '', ...
    'feature_type', '', ...
    'direction', '', ...
    'target', '', ...
    'source', '', ...
    'n_subjects_total', NaN, ...
    'n_subjects_effective', NaN, ...
    'friedman_chi2', NaN, ...
    'p', NaN, ...
    'q_primary', NaN, ...
    'sig_primary', false, ...
    'q_sensitivity', NaN, ...
    'sig_sensitivity', false, ...
    'correction_family', '', ...
    'n_tests_family', NaN, ...
    'correction_method_primary', '', ...
    'correction_method_sensitivity', '' ...
);
end

function r = init_pair_row()
r = struct(...
    'domain', '', ...
    'band', '', ...
    'measure_family', '', ...
    'feature_type', '', ...
    'direction', '', ...
    'target', '', ...
    'source', '', ...
    'stage_a', '', ...
    'stage_a_label', '', ...
    'stage_b', '', ...
    'stage_b_label', '', ...
    'n_subjects_effective', NaN, ...
    'effect_mean', NaN, ...
    'effect_median', NaN, ...
    'rank_corr', NaN, ...
    'W', NaN, ...
    'z', NaN, ...
    'p', NaN, ...
    'q_primary', NaN, ...
    'sig_primary', false, ...
    'q_sensitivity', NaN, ...
    'sig_sensitivity', false, ...
    'correction_family', '', ...
    'n_tests_family', NaN, ...
    'correction_method_primary', '', ...
    'correction_method_sensitivity', '' ...
);
end
