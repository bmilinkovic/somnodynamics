function run_gc_mi_stats_atlas(varargin)
%RUN_GC_MI_STATS_ATLAS Run GC/MI downstream statistics for a selected atlas.
%
% Legacy pipeline remains unchanged; Yeo7 outputs are written separately.

if ~exist('var_to_siogc','file')
    error('Required MVGC functions not found on path. Run startup.m first.');
end

p = inputParser;
addParameter(p, 'atlas', 'yeo7');
addParameter(p, 'pwcgc_dir', '');
addParameter(p, 'output_dir', '');
addParameter(p, 'fs', 256);
parse(p, varargin{:});
cfg = p.Results;

atlas = lower(string(cfg.atlas));
if strlength(cfg.pwcgc_dir) == 0
    if atlas == "yeo7"
        cfg.pwcgc_dir = '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/pwcgc/';
    else
        cfg.pwcgc_dir = '/Users/borjan/code/healthy-sleep-project/code/results/pwcgc/';
    end
end

if strlength(cfg.output_dir) == 0
    if atlas == "yeo7"
        cfg.output_dir = '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/';
    else
        cfg.output_dir = '/Users/borjan/code/healthy-sleep-project/code/results/new_results_GC/';
    end
end

mi_dir = fullfile(cfg.output_dir, 'mi');
if ~exist(cfg.output_dir, 'dir'), mkdir(cfg.output_dir); end
if ~exist(mi_dir, 'dir'), mkdir(mi_dir); end

bands = {'delta', 'theta', 'alpha', 'sigma', 'beta'};
freq_bands = [1 4; 4 8; 8 12; 12 16; 16 30];
wake_condition = 'W_';

gc_in_data = struct();
gc_out_data = struct();
gc_in_time_data = struct();
gc_out_time_data = struct();
mi_data = struct();
cmi_data = struct();
freq_cmi_data = struct();
labels_order = {};

files = dir(fullfile(cfg.pwcgc_dir, 'pwcgc_params_EPCTL*_*.mat'));
if isempty(files)
    error('No pwcgc parameter files found in %s', cfg.pwcgc_dir);
end

for f = 1:length(files)
    file_name = files(f).name;
    parsed = parse_pwcgc_filename(file_name);
    if parsed.valid == false
        warning('Skipping non-matching file: %s', file_name);
        continue;
    end
    participant_id = parsed.subject_id;
    condition = parsed.condition;

    data = load(fullfile(cfg.pwcgc_dir, file_name));
    if ~isfield(data, 'CAK') || ~isfield(data, 'V')
        warning('Missing CAK/V in %s', file_name);
        continue;
    end

    if isempty(labels_order) && isfield(data, 'labels')
        labels_order = normalize_labels(data.labels);
    end

    fres = 256;
    if isfield(data, 'H')
        H = data.H;
        if isfield(data, 'fres')
            fres = data.fres;
        end
    else
        H = var2trfun(data.CAK, fres);
    end
    freqs = sfreqs(fres, cfg.fs);

    gc_in_all_freqs = var_to_siogc(data.CAK, data.V, 'in', fres, H);
    gc_out_all_freqs = var_to_siogc(data.CAK, data.V, 'out', fres, H);
    gc_in_time_data.(participant_id).(condition) = mean(gc_in_all_freqs, 2);
    gc_out_time_data.(participant_id).(condition) = mean(gc_out_all_freqs, 2);

    for b = 1:length(bands)
        band_name = bands{b};
        freq_range = freq_bands(b, :);
        freq_idx = find(freqs >= freq_range(1) & freqs <= freq_range(2));
        gc_in_band = mean(gc_in_all_freqs(:, freq_idx), 2);
        gc_out_band = mean(gc_out_all_freqs(:, freq_idx), 2);
        gc_in_data.(participant_id).(condition).(band_name) = num2cell(gc_in_band);
        gc_out_data.(participant_id).(condition).(band_name) = num2cell(gc_out_band);
    end

    mi_values = cov_to_iomi(data.V);
    cmi_values = cpsd_to_scwiomi(var_to_cpsd(H, data.V, fres));
    mi_data.(participant_id).(condition) = mi_values;
    cmi_data.(participant_id).(condition) = cmi_values;
    for b = 1:length(bands)
        band_name = bands{b};
        freq_range = freq_bands(b, :);
        freq_idx = find(freqs >= freq_range(1) & freqs <= freq_range(2));
        freq_cmi_data.(participant_id).(condition).(band_name) = mean(cmi_values(:, freq_idx), 2);
    end
end

save(fullfile(cfg.output_dir, 'gcin_freq_subject_data_struct.mat'), 'gc_in_data');
save(fullfile(cfg.output_dir, 'gcout_freq_subject_data_struct.mat'), 'gc_out_data');
save(fullfile(cfg.output_dir, 'gcin_time_subject_data_struct.mat'), 'gc_in_time_data');
save(fullfile(cfg.output_dir, 'gcout_time_subject_data_struct.mat'), 'gc_out_time_data');
save(fullfile(mi_dir, 'mi_data.mat'), 'mi_data');
save(fullfile(mi_dir, 'cmi_data.mat'), 'cmi_data');
save(fullfile(mi_dir, 'freq_cmi_data.mat'), 'freq_cmi_data');
save(fullfile(cfg.output_dir, 'labels_order.mat'), 'labels_order');

participants = participants_with_conditions(gc_in_data, wake_condition);
if isempty(participants)
    error('No participants found with wake condition "%s".', wake_condition);
end

conditions = available_conditions(gc_in_data, participants{1});
conditions = conditions(~strcmp(conditions, wake_condition));

sample_band = bands{1};
num_regions = length(gc_in_data.(participants{1}).(wake_condition).(sample_band));

mw_in_effect_sizes = struct();
mw_out_effect_sizes = struct();
mw_in_time_effect_sizes = struct();
mw_out_time_effect_sizes = struct();

for c = 1:length(conditions)
    condition = conditions{c};
    valid_participants = participants_with_conditions(gc_in_data, wake_condition, condition);

    for b = 1:length(bands)
        band_name = bands{b};
        for r = 1:num_regions
            group1_in = cellfun(@(s) gc_in_data.(s).(wake_condition).(band_name){r}, ...
                valid_participants, 'UniformOutput', false);
            group2_in = cellfun(@(s) gc_in_data.(s).(condition).(band_name){r}, ...
                valid_participants, 'UniformOutput', false);
            [z_in, U_in, p_in] = run_mwu(group1_in, group2_in);
            mw_in_effect_sizes.(condition).(band_name)(r).z = z_in;
            mw_in_effect_sizes.(condition).(band_name)(r).U = U_in;
            mw_in_effect_sizes.(condition).(band_name)(r).p_value = p_in;

            group1_out = cellfun(@(s) gc_out_data.(s).(wake_condition).(band_name){r}, ...
                valid_participants, 'UniformOutput', false);
            group2_out = cellfun(@(s) gc_out_data.(s).(condition).(band_name){r}, ...
                valid_participants, 'UniformOutput', false);
            [z_out, U_out, p_out] = run_mwu(group1_out, group2_out);
            mw_out_effect_sizes.(condition).(band_name)(r).z = z_out;
            mw_out_effect_sizes.(condition).(band_name)(r).U = U_out;
            mw_out_effect_sizes.(condition).(band_name)(r).p_value = p_out;
        end
    end

    for r = 1:num_regions
        group1_in_t = cellfun(@(s) gc_in_time_data.(s).(wake_condition)(r), ...
            valid_participants, 'UniformOutput', false);
        group2_in_t = cellfun(@(s) gc_in_time_data.(s).(condition)(r), ...
            valid_participants, 'UniformOutput', false);
        [z_in_t, U_in_t, p_in_t] = run_mwu(group1_in_t, group2_in_t);
        mw_in_time_effect_sizes.(condition)(r).z = z_in_t;
        mw_in_time_effect_sizes.(condition)(r).U = U_in_t;
        mw_in_time_effect_sizes.(condition)(r).p_value = p_in_t;

        group1_out_t = cellfun(@(s) gc_out_time_data.(s).(wake_condition)(r), ...
            valid_participants, 'UniformOutput', false);
        group2_out_t = cellfun(@(s) gc_out_time_data.(s).(condition)(r), ...
            valid_participants, 'UniformOutput', false);
        [z_out_t, U_out_t, p_out_t] = run_mwu(group1_out_t, group2_out_t);
        mw_out_time_effect_sizes.(condition)(r).z = z_out_t;
        mw_out_time_effect_sizes.(condition)(r).U = U_out_t;
        mw_out_time_effect_sizes.(condition)(r).p_value = p_out_t;
    end
end

save(fullfile(cfg.output_dir, 'gc_in_freq_mw_effect_sizes.mat'), 'mw_in_effect_sizes');
save(fullfile(cfg.output_dir, 'gc_out_freq_mw_effect_sizes.mat'), 'mw_out_effect_sizes');
save(fullfile(cfg.output_dir, 'gc_in_time_mw_effect_sizes.mat'), 'mw_in_time_effect_sizes');
save(fullfile(cfg.output_dir, 'gc_out_time_mw_effect_sizes.mat'), 'mw_out_time_effect_sizes');

mi_time_effect_sizes = struct();
mi_freq_effect_sizes = struct();

for c = 1:length(conditions)
    condition = conditions{c};
    valid_participants = participants_with_conditions(mi_data, wake_condition, condition);
    if isempty(valid_participants)
        continue;
    end

    for r = 1:num_regions
        group1_mi = cellfun(@(s) mi_data.(s).(wake_condition)(r), ...
            valid_participants, 'UniformOutput', false);
        group2_mi = cellfun(@(s) mi_data.(s).(condition)(r), ...
            valid_participants, 'UniformOutput', false);
        [z_mi, U_mi, p_mi] = run_mwu(group1_mi, group2_mi);
        region_field = sprintf('Region_%d', r);
        mi_time_effect_sizes.(condition).(region_field).z = z_mi;
        mi_time_effect_sizes.(condition).(region_field).U = U_mi;
        mi_time_effect_sizes.(condition).(region_field).p_value = p_mi;
    end

    for b = 1:length(bands)
        band_name = bands{b};
        for r = 1:num_regions
            group1_cmi = cellfun(@(s) freq_cmi_data.(s).(wake_condition).(band_name)(r), ...
                valid_participants, 'UniformOutput', false);
            group2_cmi = cellfun(@(s) freq_cmi_data.(s).(condition).(band_name)(r), ...
                valid_participants, 'UniformOutput', false);
            [z_cmi, U_cmi, p_cmi] = run_mwu(group1_cmi, group2_cmi);
            region_field = sprintf('Region_%d', r);
            mi_freq_effect_sizes.(condition).(band_name).(region_field).z = z_cmi;
            mi_freq_effect_sizes.(condition).(band_name).(region_field).U = U_cmi;
            mi_freq_effect_sizes.(condition).(band_name).(region_field).p_value = p_cmi;
        end
    end
end

save(fullfile(mi_dir, 'mi_time_effect_sizes.mat'), 'mi_time_effect_sizes');
save(fullfile(mi_dir, 'mi_freq_effect_sizes.mat'), 'mi_freq_effect_sizes');

save_gc_z_matrices(cfg.output_dir, mw_in_effect_sizes, mw_out_effect_sizes, ...
    mw_in_time_effect_sizes, mw_out_time_effect_sizes, conditions, bands, num_regions);
save_mi_z_matrices(mi_dir, mi_time_effect_sizes, mi_freq_effect_sizes, conditions, bands, num_regions);

fprintf('Atlas "%s" GC/MI downstream analysis complete. Output: %s\n', atlas, cfg.output_dir);
end

function parsed = parse_pwcgc_filename(file_name)
parsed.valid = false;
parsed.subject_id = '';
parsed.condition = '';
tokens = regexp(file_name, 'pwcgc_params_(EPCTL\d{2})_([A-Za-z0-9_]+)\.mat', 'tokens', 'once');
if isempty(tokens)
    return;
end
parsed.subject_id = tokens{1};
parsed.condition = tokens{2};
parsed.valid = true;
end

function labels = normalize_labels(raw_labels)
if iscell(raw_labels)
    labels = raw_labels(:)';
elseif isstring(raw_labels)
    labels = cellstr(raw_labels);
elseif ischar(raw_labels)
    labels = cellstr(raw_labels);
else
    labels = {};
end
end

function participants = participants_with_conditions(data_struct, varargin)
participants = fieldnames(data_struct);
if isempty(varargin)
    return;
end

required_conditions = varargin;
keep = false(size(participants));
for i = 1:length(participants)
    pid = participants{i};
    has_all = true;
    for c = 1:length(required_conditions)
        if ~isfield(data_struct.(pid), required_conditions{c})
            has_all = false;
            break;
        end
    end
    keep(i) = has_all;
end
participants = participants(keep);
end

function conds = available_conditions(data_struct, participant_id)
conds = fieldnames(data_struct.(participant_id));
end

function [z, U, p] = run_mwu(group1, group2)
if isempty(group1) || isempty(group2) || length(group1) ~= length(group2)
    z = NaN; U = NaN; p = NaN;
    return;
end
[z, U] = mann_whitney_group(group1, group2);
p = erfc(abs(z) / sqrt(2));
end

function save_gc_z_matrices(output_dir, mw_in_effect_sizes, mw_out_effect_sizes, ...
    mw_in_time_effect_sizes, mw_out_time_effect_sizes, conditions, bands, num_regions)

num_conditions = length(conditions);
z_matrix_in_time = nan(num_regions, num_conditions);
z_matrix_out_time = nan(num_regions, num_conditions);

for cond_idx = 1:num_conditions
    condition = conditions{cond_idx};
    if isfield(mw_in_time_effect_sizes, condition)
        z_vector = [mw_in_time_effect_sizes.(condition).z];
        z_matrix_in_time(:, cond_idx) = z_vector(:);
    end
    if isfield(mw_out_time_effect_sizes, condition)
        z_vector = [mw_out_time_effect_sizes.(condition).z];
        z_matrix_out_time(:, cond_idx) = z_vector(:);
    end
end
save(fullfile(output_dir, 'gc_in_time_zscores_matrix.mat'), 'z_matrix_in_time');
save(fullfile(output_dir, 'gc_out_time_zscores_matrix.mat'), 'z_matrix_out_time');

for b = 1:length(bands)
    band_name = bands{b};
    z_matrix_in_freq = nan(num_regions, num_conditions);
    z_matrix_out_freq = nan(num_regions, num_conditions);
    for cond_idx = 1:num_conditions
        condition = conditions{cond_idx};
        if isfield(mw_in_effect_sizes, condition) && isfield(mw_in_effect_sizes.(condition), band_name)
            z_vector = [mw_in_effect_sizes.(condition).(band_name).z];
            z_matrix_in_freq(:, cond_idx) = z_vector(:);
        end
        if isfield(mw_out_effect_sizes, condition) && isfield(mw_out_effect_sizes.(condition), band_name)
            z_vector = [mw_out_effect_sizes.(condition).(band_name).z];
            z_matrix_out_freq(:, cond_idx) = z_vector(:);
        end
    end
    save(fullfile(output_dir, sprintf('gc_in_freq_zscores_%s.mat', band_name)), 'z_matrix_in_freq');
    save(fullfile(output_dir, sprintf('gc_out_freq_zscores_%s.mat', band_name)), 'z_matrix_out_freq');
end
end

function save_mi_z_matrices(mi_dir, mi_time_effect_sizes, mi_freq_effect_sizes, ...
    conditions, bands, num_regions)

num_conditions = length(conditions);
z_matrix_mi_time = nan(num_regions, num_conditions);
for cond_idx = 1:num_conditions
    condition = conditions{cond_idx};
    for region_idx = 1:num_regions
        region_field = sprintf('Region_%d', region_idx);
        if isfield(mi_time_effect_sizes, condition) && isfield(mi_time_effect_sizes.(condition), region_field)
            z_matrix_mi_time(region_idx, cond_idx) = mi_time_effect_sizes.(condition).(region_field).z;
        end
    end
end
save(fullfile(mi_dir, 'mi_time_zscores_matrix.mat'), 'z_matrix_mi_time');

for b = 1:length(bands)
    band_name = bands{b};
    z_matrix_mi_freq = nan(num_regions, num_conditions);
    for cond_idx = 1:num_conditions
        condition = conditions{cond_idx};
        for region_idx = 1:num_regions
            region_field = sprintf('Region_%d', region_idx);
            if isfield(mi_freq_effect_sizes, condition) && ...
                    isfield(mi_freq_effect_sizes.(condition), band_name) && ...
                    isfield(mi_freq_effect_sizes.(condition).(band_name), region_field)
                z_matrix_mi_freq(region_idx, cond_idx) = ...
                    mi_freq_effect_sizes.(condition).(band_name).(region_field).z;
            end
        end
    end
    save(fullfile(mi_dir, sprintf('mi_freq_zscores_%s.mat', band_name)), 'z_matrix_mi_freq');
end
end
