function run_stage_ssdi_anatomical7_pipeline(varargin)
% Master runner for DMT-faithful SomnoDynamics stage SSDI (anatomical7).

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);

cfg = sd_ssdi_config_anatomical7();
cfg = apply_overrides(cfg, varargin{:});

prepare_dirs(cfg);
initialize_toolboxes(cfg);
verify_inputs(cfg);

log_file = fullfile(cfg.logs_dir, sprintf('stage_ssdi_pipeline_%s.log', datestr(now, 'yyyymmdd_HHMMSS')));
diary(log_file);

fprintf('\n=== SomnoDynamics stage-based SSDI anatomical7 pipeline ===\n');
fprintf('Project root : %s\n', cfg.project_root);
fprintf('Input root   : %s\n', cfg.subject_measure_root);
fprintf('Results root : %s\n', cfg.results_root);
fprintf('Start time   : %s\n\n', datestr(now, 31));

sd_ssdi01_compute_stage_spectral(cfg);

fprintf('\nPipeline completed successfully at %s\n', datestr(now, 31));
fprintf('Log file: %s\n', log_file);
diary off;

end

function cfg = apply_overrides(cfg, varargin)
if isempty(varargin)
    return;
end

p = inputParser;
addParameter(p, 'input_root', cfg.subject_measure_root);
addParameter(p, 'results_root', cfg.results_root);
addParameter(p, 'recompute', cfg.recompute);
addParameter(p, 'parallel_enable', cfg.parallel.enable);
addParameter(p, 'parallel_workers', cfg.parallel.n_workers);
addParameter(p, 'max_files', cfg.max_files);
addParameter(p, 'nrunsp', cfg.nrunsp);
addParameter(p, 'nitersp', cfg.nitersp);
addParameter(p, 'niterso', cfg.niterso);
parse(p, varargin{:});

cfg.subject_measure_root = p.Results.input_root;
cfg.results_root = p.Results.results_root;
cfg.recompute = p.Results.recompute;
cfg.parallel.enable = p.Results.parallel_enable;
cfg.parallel.n_workers = p.Results.parallel_workers;
cfg.max_files = p.Results.max_files;
cfg.nrunsp = p.Results.nrunsp;
cfg.nitersp = p.Results.nitersp;
cfg.niterso = p.Results.niterso;

cfg.logs_dir = fullfile(cfg.results_root, 'logs');
cfg.figures_dir = fullfile(cfg.results_root, 'figures');
cfg.diag_dir = fullfile(cfg.results_root, 'diagnostics');
cfg.summary_dir = fullfile(cfg.results_root, 'summaries');
cfg.spectral_dir = fullfile(cfg.results_root, 'spectral');
end

function prepare_dirs(cfg)
outdirs = {cfg.results_root, cfg.spectral_dir, cfg.logs_dir, cfg.figures_dir, cfg.summary_dir, cfg.diag_dir};
for s = 1:numel(cfg.stages)
    outdirs{end+1} = fullfile(cfg.results_root, cfg.stages{s}); %#ok<AGROW>
end
for i = 1:numel(outdirs)
    if exist(outdirs{i}, 'dir') ~= 7
        mkdir(outdirs{i});
    end
end
end

function initialize_toolboxes(cfg)
if exist('opt_gd_dds_mruns', 'file') == 2 && exist('var_to_ss', 'file') == 2
    return;
end

run_all_root = fullfile(cfg.project_root, 'code', 'run_all');
addpath(run_all_root);
setup_matlab_toolboxes(cfg.project_root);

if exist('opt_gd_dds_mruns', 'file') ~= 2
    error('SSDI toolbox not initialized (opt_gd_dds_mruns unresolved).');
end
if exist('var_to_ss', 'file') ~= 2
    error('MVGC2 toolbox not initialized (var_to_ss unresolved).');
end
end

function verify_inputs(cfg)
files = dir(fullfile(cfg.subject_measure_root, '*', '*_subject_measures.mat'));
if isempty(files)
    error('No stage subject-measure files found in %s', cfg.subject_measure_root);
end

counts = zeros(numel(cfg.stages), 1);
for k = 1:numel(files)
    stage = string(files(k).folder);
    for s = 1:numel(cfg.stages)
        if contains(stage, filesep + string(cfg.stages{s}) + filesep) || endsWith(stage, filesep + string(cfg.stages{s}))
            counts(s) = counts(s) + 1;
        end
    end
end

fprintf('[SomnoDynamics SSDI] input files found: %d\n', numel(files));
for s = 1:numel(cfg.stages)
    fprintf('  Stage %-2s (%-4s): %d files\n', cfg.stages{s}, cfg.stage_display.(cfg.stages{s}), counts(s));
end
end
