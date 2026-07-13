function run_stage_mvgc2_anatomical7_pipeline(varargin)
% Master runner for SomnoDynamics stage-based MVGC2 (anatomical7) pipeline.

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);

cfg = sd_gc_config_anatomical7();
cfg = apply_overrides(cfg, varargin{:});

prepare_dirs(cfg);
initialize_mvgc2(cfg);
verify_source_inputs(cfg);

log_file = fullfile(cfg.paths.logs, sprintf('stage_mvgc2_pipeline_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
diary(log_file);

fprintf('\n=== SomnoDynamics stage-based MVGC2 anatomical7 pipeline ===\n');
fprintf('Project root : %s\n', cfg.project_root);
fprintf('Source root  : %s\n', cfg.paths.source_root);
fprintf('Results root : %s\n', cfg.paths.results_root);
fprintf('Start time   : %s\n\n', datestr(now, 31));

fprintf('[1/2] Subject measures by stage\n');
sd_gc01_compute_subject_measures_by_stage(cfg);

fprintf('\n[2/2] Repeated-measures group statistics\n');
sd_gc02_group_stats_repeated_measures(cfg);

fprintf('\nPipeline completed successfully at %s\n', datestr(now, 31));
fprintf('Log file: %s\n', log_file);

diary off;

end

function cfg = apply_overrides(cfg, varargin)
if isempty(varargin)
    return;
end

p = inputParser;
addParameter(p, 'source_root', cfg.paths.source_root);
addParameter(p, 'results_root', cfg.paths.results_root);
addParameter(p, 'recompute_subject_measures', cfg.recompute_subject_measures);
addParameter(p, 'recompute_group_stats', cfg.recompute_group_stats);
addParameter(p, 'parallel_enable', cfg.parallel.enable);
addParameter(p, 'parallel_workers', cfg.parallel.n_workers);
addParameter(p, 'max_files', cfg.max_files);
addParameter(p, 'normalization_mode', cfg.normalization_mode);
addParameter(p, 'normalization_eps', cfg.normalization_eps);
addParameter(p, 'fs_original', cfg.fs_original);
addParameter(p, 'fs_target', cfg.fs_target);
addParameter(p, 'fres', cfg.fres);
addParameter(p, 'momax', cfg.momax);
addParameter(p, 'mosel', cfg.mosel);
parse(p, varargin{:});

cfg.paths.source_root = p.Results.source_root;
cfg.paths.results_root = p.Results.results_root;
cfg.recompute_subject_measures = p.Results.recompute_subject_measures;
cfg.recompute_group_stats = p.Results.recompute_group_stats;
cfg.parallel.enable = p.Results.parallel_enable;
cfg.parallel.n_workers = p.Results.parallel_workers;
cfg.max_files = p.Results.max_files;
cfg.normalization_mode = p.Results.normalization_mode;
cfg.normalization_eps = p.Results.normalization_eps;
cfg.fs_original = p.Results.fs_original;
cfg.fs_target = p.Results.fs_target;
cfg.fres = p.Results.fres;
cfg.momax = p.Results.momax;
cfg.mosel = char(p.Results.mosel);

cfg.paths.data_intermediate = fullfile(cfg.paths.results_root, 'data_intermediate');
cfg.paths.subject_measures = fullfile(cfg.paths.data_intermediate, 'subject_measures');
cfg.paths.stats = fullfile(cfg.paths.results_root, 'stats');
cfg.paths.figures = fullfile(cfg.paths.results_root, 'figures');
cfg.paths.logs = fullfile(cfg.paths.results_root, 'logs');
cfg.paths.reports = fullfile(cfg.paths.results_root, 'reports');

end

function prepare_dirs(cfg)
outdirs = {
    cfg.paths.results_root,
    cfg.paths.data_intermediate,
    cfg.paths.subject_measures,
    cfg.paths.stats,
    cfg.paths.figures,
    cfg.paths.logs,
    cfg.paths.reports
};
for s = 1:numel(cfg.stages)
    outdirs{end+1} = fullfile(cfg.paths.subject_measures, cfg.stages{s}); %#ok<AGROW>
end

for k = 1:numel(outdirs)
    if exist(outdirs{k}, 'dir') ~= 7
        mkdir(outdirs{k});
    end
end
end

function initialize_mvgc2(cfg)
if exist('tsdata_to_varmo', 'file') == 2
    return;
end

mvgc_root = cfg.paths.mvgc_root;
if exist(mvgc_root, 'dir') ~= 7
    error('Vendored MVGC2 not found: %s', mvgc_root);
end

startup_file = fullfile(mvgc_root, 'startup.m');
if exist(startup_file, 'file') ~= 2
    error('MVGC2 startup.m not found at %s', startup_file);
end

old_pwd = pwd;
cleanup = onCleanup(@() cd(old_pwd)); %#ok<NASGU>
cd(mvgc_root);
run(startup_file);

if exist('tsdata_to_varmo', 'file') ~= 2
    error('MVGC2 initialization failed (tsdata_to_varmo unresolved).');
end

fprintf('[SomnoDynamics] MVGC2 initialized from %s\n', mvgc_root);
end

function verify_source_inputs(cfg)
files = dir(fullfile(cfg.paths.source_root, 'EPCTL*', 'mat', '*_source_time_series.mat'));
if isempty(files)
    error('No source files found under %s', cfg.paths.source_root);
end

counts = zeros(numel(cfg.stages), 1);
for k = 1:numel(files)
    parsed = parse_source_filename(files(k).name);
    if ~parsed.valid
        continue;
    end
    idx = find(strcmp(cfg.stages, parsed.stage), 1);
    if ~isempty(idx)
        counts(idx) = counts(idx) + 1;
    end
end

fprintf('[SomnoDynamics] Source input files found: %d\n', numel(files));
for s = 1:numel(cfg.stages)
    fprintf('  Stage %-2s (%-4s): %d files\n', cfg.stages{s}, cfg.stage_display.(cfg.stages{s}), counts(s));
end

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
