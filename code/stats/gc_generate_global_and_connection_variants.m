function gc_generate_global_and_connection_variants(results_root, varargin)
%GC_GENERATE_GLOBAL_AND_CONNECTION_VARIANTS
% Build two post-hoc GC statistics variants:
% 1) Stage-contrast stats with one global correction pool across all tests.
% 2) One-sample connection-presence stats per stage/band (H1: median > 0),
%    corrected globally across all such tests.
%
% Usage:
%   gc_generate_global_and_connection_variants('/abs/path/to/gc_root')

p = inputParser;
addRequired(p, 'results_root', @(x) ischar(x) || isstring(x));
addParameter(p, 'alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'method_primary', 'BH', @(x) ischar(x) || isstring(x));
addParameter(p, 'method_sensitivity', 'BY', @(x) ischar(x) || isstring(x));
parse(p, results_root, varargin{:});
opt = p.Results;

results_root = char(results_root);
stats_dir = fullfile(results_root, 'stats');
subject_root = fullfile(results_root, 'data_intermediate', 'subject_measures');

pairwise_csv = fullfile(stats_dir, 'pairwise_long.csv');
omnibus_csv = fullfile(stats_dir, 'omnibus_long.csv');
if exist(pairwise_csv, 'file') ~= 2 || exist(omnibus_csv, 'file') ~= 2
    error('Missing stats CSVs in %s', stats_dir);
end
if exist(subject_root, 'dir') ~= 7
    error('Missing subject measures dir: %s', subject_root);
end

fprintf('\n[GC variants] results_root: %s\n', results_root);

% -------------------------------------------------------------------------
% Variant A: stage-contrast stats with global correction across all tests.
% -------------------------------------------------------------------------
T_pair = readtable(pairwise_csv, 'TextType', 'string');
T_omni = readtable(omnibus_csv, 'TextType', 'string');

if ~ismember('p', T_pair.Properties.VariableNames) || ~ismember('p', T_omni.Properties.VariableNames)
    error('Expected column ''p'' in pairwise/omnibus tables.');
end

T_pair_global = T_pair;
T_omni_global = T_omni;

if ismember('q_primary', T_pair_global.Properties.VariableNames)
    T_pair_global.q_family_primary = T_pair_global.q_primary;
end
if ismember('sig_primary', T_pair_global.Properties.VariableNames)
    T_pair_global.sig_family_primary = T_pair_global.sig_primary;
end
if ismember('q_sensitivity', T_pair_global.Properties.VariableNames)
    T_pair_global.q_family_sensitivity = T_pair_global.q_sensitivity;
end
if ismember('sig_sensitivity', T_pair_global.Properties.VariableNames)
    T_pair_global.sig_family_sensitivity = T_pair_global.sig_sensitivity;
end

if ismember('q_primary', T_omni_global.Properties.VariableNames)
    T_omni_global.q_family_primary = T_omni_global.q_primary;
end
if ismember('sig_primary', T_omni_global.Properties.VariableNames)
    T_omni_global.sig_family_primary = T_omni_global.sig_primary;
end
if ismember('q_sensitivity', T_omni_global.Properties.VariableNames)
    T_omni_global.q_family_sensitivity = T_omni_global.q_sensitivity;
end
if ismember('sig_sensitivity', T_omni_global.Properties.VariableNames)
    T_omni_global.sig_family_sensitivity = T_omni_global.sig_sensitivity;
end

T_pair_global = apply_global_fdr_table(T_pair_global, opt.method_primary, opt.method_sensitivity, opt.alpha, 'GLOBAL_ALL_PAIRWISE');
T_omni_global = apply_global_fdr_table(T_omni_global, opt.method_primary, opt.method_sensitivity, opt.alpha, 'GLOBAL_ALL_OMNIBUS');

out_global = fullfile(results_root, 'stats_global_all_tests');
if exist(out_global, 'dir') ~= 7
    mkdir(out_global);
end

writetable(T_pair_global, fullfile(out_global, 'pairwise_long.csv'));
writetable(T_omni_global, fullfile(out_global, 'omnibus_long.csv'));

fid = fopen(fullfile(out_global, 'summary.txt'), 'w');
if fid >= 0
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'GC stage-contrast global-correction variant\n');
    fprintf(fid, '==========================================\n');
    fprintf(fid, 'Input stats dir: %s\n', stats_dir);
    fprintf(fid, 'Global pool (pairwise): %d tests\n', sum(isfinite(T_pair.p)));
    fprintf(fid, 'Global pool (omnibus) : %d tests\n', sum(isfinite(T_omni.p)));
    fprintf(fid, 'Primary correction: %s\n', upper(string(opt.method_primary)));
    fprintf(fid, 'Sensitivity correction: %s\n', upper(string(opt.method_sensitivity)));
    fprintf(fid, 'Alpha: %.3f\n\n', opt.alpha);
    fprintf(fid, 'Pairwise significant (primary): %d\n', sum(T_pair_global.sig_primary > 0));
    fprintf(fid, 'Pairwise significant (sensitivity): %d\n', sum(T_pair_global.sig_sensitivity > 0));
    fprintf(fid, 'Omnibus significant (primary): %d\n', sum(T_omni_global.sig_primary > 0));
    fprintf(fid, 'Omnibus significant (sensitivity): %d\n', sum(T_omni_global.sig_sensitivity > 0));
end

% -------------------------------------------------------------------------
% Variant B: per-stage/per-band connection presence (one-sample > 0).
% -------------------------------------------------------------------------
T_conn = build_connection_presence_table(subject_root);
T_conn = apply_global_fdr_table(T_conn, opt.method_primary, opt.method_sensitivity, opt.alpha, 'GLOBAL_ALL_STAGE_BAND_CONNECTION');

out_conn = fullfile(results_root, 'stats_connection_presence_global_all_tests');
if exist(out_conn, 'dir') ~= 7
    mkdir(out_conn);
end

writetable(T_conn, fullfile(out_conn, 'connection_presence_long.csv'));
writetable(T_conn(T_conn.sig_primary > 0, :), fullfile(out_conn, 'connection_presence_significant_primary.csv'));
writetable(T_conn(T_conn.sig_sensitivity > 0, :), fullfile(out_conn, 'connection_presence_significant_sensitivity.csv'));

fid = fopen(fullfile(out_conn, 'summary.txt'), 'w');
if fid >= 0
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'GC connection-presence variant (one-sample signed-rank)\n');
    fprintf(fid, '======================================================\n');
    fprintf(fid, 'Input subject-measures root: %s\n', subject_root);
    fprintf(fid, 'Total tests in global pool: %d\n', sum(isfinite(T_conn.p)));
    fprintf(fid, 'Primary correction: %s\n', upper(string(opt.method_primary)));
    fprintf(fid, 'Sensitivity correction: %s\n', upper(string(opt.method_sensitivity)));
    fprintf(fid, 'Alpha: %.3f\n\n', opt.alpha);
    fprintf(fid, 'Significant (primary): %d\n', sum(T_conn.sig_primary > 0));
    fprintf(fid, 'Significant (sensitivity): %d\n', sum(T_conn.sig_sensitivity > 0));
    fprintf(fid, '\nNull hypothesis: median GC connection strength <= 0.\n');
    fprintf(fid, 'Alternative hypothesis: median GC connection strength > 0.\n');
end

fprintf('[GC variants] wrote: %s\n', out_global);
fprintf('[GC variants] wrote: %s\n', out_conn);
end

function T = build_connection_presence_table(subject_root)
stages = {'W_', 'N1', 'N2', 'N3', 'R_'};
stage_labels = containers.Map(...
    {'W_', 'N1', 'N2', 'N3', 'R_'}, ...
    {'Wake', 'N1', 'N2', 'N3', 'REM'});

subject_sets = cell(numel(stages), 1);
for s = 1:numel(stages)
    files = dir(fullfile(subject_root, stages{s}, '*_subject_measures.mat'));
    subject_sets{s} = erase({files.name}, '_subject_measures.mat');
end
subjects = subject_sets{1};
for s = 2:numel(subject_sets)
    subjects = intersect(subjects, subject_sets{s}, 'stable');
end
if isempty(subjects)
    error('No subjects with all required stages in %s', subject_root);
end

ref_file = fullfile(subject_root, stages{1}, sprintf('%s_subject_measures.mat', subjects{1}));
Sref = load(ref_file, 'subject_measure');
Mref = Sref.subject_measure;
labels = Mref.meta.group_names(:);
nreg = numel(labels);
band_info = Mref.meta.band_info;

rows = repmat(init_conn_row(), 0, 1);

nsub = numel(subjects);
for st = 1:numel(stages)
    stage = stages{st};
    Mstage = cell(nsub, 1);
    for si = 1:nsub
        f = fullfile(subject_root, stage, sprintf('%s_subject_measures.mat', subjects{si}));
        S = load(f, 'subject_measure');
        Mstage{si} = S.subject_measure;
    end

    rows = append_time_rows(rows, Mstage, labels, stage, stage_labels(stage));
    rows = append_spectral_rows(rows, Mstage, labels, band_info, stage, stage_labels(stage));
end

T = struct2table(rows, 'AsArray', true);
end

function rows = append_time_rows(rows, Mstage, labels, stage, stage_label)
nsub = numel(Mstage);
nreg = numel(labels);

for r = 1:nreg
    x = nan(nsub, 1);
    for si = 1:nsub
        x(si) = Mstage{si}.measures.directed.time.src_in(r);
    end
    rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'time', 'all', 'directed_src_in', 'nodewise', 'in', labels{r}, ''); %#ok<AGROW>

    x = nan(nsub, 1);
    for si = 1:nsub
        x(si) = Mstage{si}.measures.directed.time.src_out(r);
    end
    rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'time', 'all', 'directed_src_out', 'nodewise', 'out', labels{r}, ''); %#ok<AGROW>

    x = nan(nsub, 1);
    for si = 1:nsub
        x(si) = Mstage{si}.measures.undirected.time.src(r);
    end
    rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'time', 'all', 'undirected_src', 'nodewise', 'undirected', labels{r}, ''); %#ok<AGROW>
end

for to = 1:nreg
    for from = 1:nreg
        if to == from
            continue;
        end

        x = nan(nsub, 1);
        for si = 1:nsub
            x(si) = Mstage{si}.measures.directed.time.roi(to, from);
        end
        rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'time', 'all', 'directed_roi', 'edgewise', 'directed', labels{to}, labels{from}); %#ok<AGROW>

        x = nan(nsub, 1);
        for si = 1:nsub
            x(si) = Mstage{si}.measures.undirected.time.roi(to, from);
        end
        rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'time', 'all', 'undirected_roi', 'edgewise', 'undirected', labels{to}, labels{from}); %#ok<AGROW>
    end
end
end

function rows = append_spectral_rows(rows, Mstage, labels, band_info, stage, stage_label)
nsub = numel(Mstage);
nreg = numel(labels);

for b = 1:numel(band_info)
    bname = band_info(b).name;
    mask = logical(band_info(b).mask(:));

    for r = 1:nreg
        x = nan(nsub, 1);
        for si = 1:nsub
            v = squeeze(Mstage{si}.measures.directed.spectral.src_in(r, :));
            x(si) = mean(v(mask), 'omitnan');
        end
        rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'spectral', bname, 'directed_src_in', 'nodewise', 'in', labels{r}, ''); %#ok<AGROW>

        x = nan(nsub, 1);
        for si = 1:nsub
            v = squeeze(Mstage{si}.measures.directed.spectral.src_out(r, :));
            x(si) = mean(v(mask), 'omitnan');
        end
        rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'spectral', bname, 'directed_src_out', 'nodewise', 'out', labels{r}, ''); %#ok<AGROW>

        x = nan(nsub, 1);
        for si = 1:nsub
            v = squeeze(Mstage{si}.measures.undirected.spectral.src(r, :));
            x(si) = mean(v(mask), 'omitnan');
        end
        rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'spectral', bname, 'undirected_src', 'nodewise', 'undirected', labels{r}, ''); %#ok<AGROW>
    end

    for to = 1:nreg
        for from = 1:nreg
            if to == from
                continue;
            end

            x = nan(nsub, 1);
            for si = 1:nsub
                v = squeeze(Mstage{si}.measures.directed.spectral.roi(to, from, :));
                x(si) = mean(v(mask), 'omitnan');
            end
            rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'spectral', bname, 'directed_roi', 'edgewise', 'directed', labels{to}, labels{from}); %#ok<AGROW>

            x = nan(nsub, 1);
            for si = 1:nsub
                v = squeeze(Mstage{si}.measures.undirected.spectral.roi(to, from, :));
                x(si) = mean(v(mask), 'omitnan');
            end
            rows(end+1, 1) = build_conn_row(x, stage, stage_label, 'spectral', bname, 'undirected_roi', 'edgewise', 'undirected', labels{to}, labels{from}); %#ok<AGROW>
        end
    end
end
end

function row = build_conn_row(x, stage, stage_label, domain, band, measure_family, feature_type, direction, target, source)
row = init_conn_row();
row.stage = string(stage);
row.stage_label = string(stage_label);
row.domain = string(domain);
row.band = string(band);
row.measure_family = string(measure_family);
row.feature_type = string(feature_type);
row.direction = string(direction);
row.target = string(target);
row.source = string(source);
row.n_subjects_total = numel(x);

st = one_sample_signedrank_stats(x);
row.n_subjects_effective = st.n_eff;
row.effect_mean = st.effect_mean;
row.effect_median = st.effect_median;
row.rank_corr = st.rank_corr;
row.W = st.W;
row.z = st.z;
row.p = st.p;

row.test_type = "one_sample_signrank_right_tailed";
row.hypothesis_null = "median_connection_strength <= 0";
row.hypothesis_alt = "median_connection_strength > 0";
end

function st = one_sample_signedrank_stats(x)
x = x(isfinite(x));
st = struct('n_eff', 0, 'effect_mean', NaN, 'effect_median', NaN, 'rank_corr', NaN, 'W', NaN, 'z', NaN, 'p', NaN);
if numel(x) < 3
    st.n_eff = sum(x ~= 0);
    st.effect_mean = mean(x, 'omitnan');
    st.effect_median = median(x, 'omitnan');
    return;
end

st.n_eff = sum(x ~= 0);
st.effect_mean = mean(x, 'omitnan');
st.effect_median = median(x, 'omitnan');

try
    [p, ~, stats] = signrank(x, 0, 'tail', 'right', 'method', 'approximate');
catch
    try
        [p, ~, stats] = signrank(x, 0, 'tail', 'right');
    catch
        p = NaN;
        stats = struct();
    end
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

function T = apply_global_fdr_table(T, method_primary, method_sensitivity, alpha, scope_name)
if ~ismember('p', T.Properties.VariableNames)
    error('Expected column ''p'' for FDR correction.');
end

pvals = T.p;
q_primary = fdr_adjust(pvals, method_primary);
q_sens = fdr_adjust(pvals, method_sensitivity);

T.q_primary = q_primary;
T.sig_primary = q_primary <= alpha;
T.q_sensitivity = q_sens;
T.sig_sensitivity = q_sens <= alpha;

n_tests = sum(isfinite(pvals));
T.n_tests_family = repmat(n_tests, height(T), 1);
T.correction_family = repmat(string(scope_name), height(T), 1);
T.correction_method_primary = repmat(string(upper(method_primary)), height(T), 1);
T.correction_method_sensitivity = repmat(string(upper(method_sensitivity)), height(T), 1);
end

function q = fdr_adjust(pvals, method)
q = nan(size(pvals));
idx = find(isfinite(pvals));
if isempty(idx)
    return;
end

p = pvals(idx);
[p_sorted, ord] = sort(p(:));
m = numel(p_sorted);

if strcmpi(method, 'BY')
    c_m = sum(1 ./ (1:m));
else
    c_m = 1;
end

ranks = (1:m)';
q_sorted = (p_sorted .* m .* c_m) ./ ranks;
q_sorted = flipud(cummin(flipud(q_sorted)));
q_sorted = min(q_sorted, 1);

q_local = nan(size(p));
q_local(ord) = q_sorted;
q(idx) = q_local;
end

function r = init_conn_row()
r = struct(...
    'stage', "", ...
    'stage_label', "", ...
    'domain', "", ...
    'band', "", ...
    'measure_family', "", ...
    'feature_type', "", ...
    'direction', "", ...
    'target', "", ...
    'source', "", ...
    'n_subjects_total', NaN, ...
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
    'correction_family', "", ...
    'n_tests_family', NaN, ...
    'correction_method_primary', "", ...
    'correction_method_sensitivity', "", ...
    'test_type', "", ...
    'hypothesis_null', "", ...
    'hypothesis_alt', "" ...
);
end
