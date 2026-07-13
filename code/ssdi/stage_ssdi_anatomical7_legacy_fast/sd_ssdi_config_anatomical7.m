function cfg = sd_ssdi_config_anatomical7()
% Configuration for SomnoDynamics stage-based SSDI (anatomical7).

cfg.project_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
results_root = sd_external_results_root(cfg.project_root);

cfg.paths.gc_subject_measures = fullfile(results_root, 'gc', 'stage_mvgc2_anatomical7', ...
    'data_intermediate', 'subject_measures');
cfg.paths.results_root = fullfile(results_root, 'ssdi', 'stage_ssdi_anatomical7');
cfg.paths.spectral_out = fullfile(cfg.paths.results_root, 'spectral');
cfg.paths.logs = fullfile(cfg.paths.results_root, 'logs');
cfg.paths.reports = fullfile(cfg.paths.results_root, 'reports');
cfg.paths.diag = fullfile(cfg.paths.results_root, 'diagnostics');

cfg.paths.ssdi_root = fullfile(cfg.project_root, 'code', 'external', 'ssdi-1');

cfg.stages = {'W_', 'N1', 'N2', 'N3', 'R_'};
cfg.stage_display = struct('W_', 'Wake', 'N1', 'N1', 'N2', 'N2', 'N3', 'N3', 'R_', 'REM');

cfg.required_regions = 7;
cfg.m_values = 2:6;

cfg.fs = 128;
cfg.fres = 128;

% Non-overlapping frequency bands requested for SomnoDynamics.
cfg.band_defs = struct( ...
    'name', {'delta','theta','alpha','sigma','beta','gamma'}, ...
    'lo',   {1,5,9,13,17,31}, ...
    'hi',   {4,8,12,16,30,44});

cfg.iseed = 0;
cfg.nrunsp = 30;
cfg.nitersp = 3000;
cfg.sig0p = 1;
cfg.gdesp = 2;
cfg.gdlsp = 2;
cfg.gdtolp = 1e-10;
cfg.histp = true;
cfg.ppp = false;

cfg.parallel = struct();
cfg.parallel.enable = true;
cfg.parallel.n_workers = [];

cfg.recompute = true;
cfg.max_files = inf;

end
