function cfg = sd_gc_config_anatomical7()
% Configuration for standalone SomnoDynamics stage-based MVGC2 analysis (DK5).

cfg.project_root = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
results_root = sd_external_results_root(cfg.project_root);

cfg.paths.source_root = fullfile(results_root, 'source_recon', 'dk5_dk68_lcmv_nai_current_epoching_2s_256hz', 'dk5');
cfg.paths.results_root = fullfile(results_root, 'gc', 'stage_mvgc2_dk5_lrt30_current_epoching_2s_256hz');

cfg.paths.data_intermediate = fullfile(cfg.paths.results_root, 'data_intermediate');
cfg.paths.subject_measures = fullfile(cfg.paths.data_intermediate, 'subject_measures');
cfg.paths.stats = fullfile(cfg.paths.results_root, 'stats');
cfg.paths.figures = fullfile(cfg.paths.results_root, 'figures');
cfg.paths.logs = fullfile(cfg.paths.results_root, 'logs');
cfg.paths.reports = fullfile(cfg.paths.results_root, 'reports');

cfg.paths.mvgc_root = fullfile(cfg.project_root, 'code', 'external', 'MVGC2');

cfg.stages = {'W_', 'N1', 'N2', 'N3', 'R_'};
cfg.stage_display = struct('W_', 'Wake', 'N1', 'N1', 'N2', 'N2', 'N3', 'N3', 'R_', 'REM');

% VAR/GC estimation settings (parity with DMT paper-style defaults).
% Sampling configuration.
cfg.fs_original = 256;
cfg.fs_target = 256;
cfg.fres = 256;
cfg.momax = 30;
cfg.moregmode = 'LWR';
cfg.mosel = 'LRT';
cfg.regmode = 'LWR';

% Multi-trial normalization strategy (Barnett 2020 style).
% 'pooled_zscore': pooled mean/variance per region across all epochs+time in a subject-stage.
% 'demean_epoch' : per-epoch demeaning only (legacy behavior).
% 'none'         : no normalization.
cfg.normalization_mode = 'pooled_zscore';
cfg.normalization_eps = 1e-12;

% Estimator-level pairwise GC significance (subject-level diagnostic).
cfg.compute_estimator_level_gc_sig = true;
cfg.alpha_estimator = 0.05;
cfg.estimator_test_stat = 'F';
cfg.estimator_mhtc = 'FDRD';

% Group-level repeated-measures settings.
cfg.alpha_group = 0.05;
cfg.fdr_primary = 'BH';
cfg.fdr_sensitivity = 'BY';

% Non-overlapping canonical sleep bands requested for SomnoDynamics DK5 GC.
% Sigma is kept separate from beta; gamma is intentionally excluded.
cfg.band_defs = {
    'delta', 1, 4;
    'theta', 5, 8;
    'alpha', 9, 12;
    'sigma', 13, 16;
    'beta', 17, 30;
};

cfg.directed_convention = 'row=to, column=from';

cfg.recompute_subject_measures = true;
cfg.recompute_group_stats = true;

cfg.parallel = struct();
cfg.parallel.enable = true;
cfg.parallel.n_workers = [];

cfg.max_files = inf; % For debugging/smoke tests only.

end
