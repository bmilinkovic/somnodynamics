function cfg = sd_ssdi_config_anatomical7()
% DMT-faithful SSDI config for standalone SomnoDynamics DK5 sleep-stage data.

cfg.project_root = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
results_root = sd_external_results_root(cfg.project_root);

cfg.results_root = fullfile(results_root, 'ssdi', 'stage_ssdi_dk5_lrt30_current_epoching_2s_256hz');
cfg.logs_dir = fullfile(cfg.results_root, 'logs');
cfg.figures_dir = fullfile(cfg.results_root, 'figures');
cfg.diag_dir = fullfile(cfg.results_root, 'diagnostics');
cfg.summary_dir = fullfile(cfg.results_root, 'summaries');
cfg.spectral_dir = fullfile(cfg.results_root, 'spectral');
cfg.subject_measure_root = fullfile(results_root, 'gc', 'stage_mvgc2_dk5_lrt30_current_epoching_2s_256hz', ...
    'data_intermediate', 'subject_measures');

cfg.stages = {'W_', 'N1', 'N2', 'N3', 'R_'};
cfg.stage_display = struct('W_', 'Wake', 'N1', 'N1', 'N2', 'N2', 'N3', 'N3', 'R_', 'REM');

cfg.required_regions = 5;
cfg.required_labels = {'Frontal','Sensorimotor','Parietal','Temporal','Occipital'};
cfg.m_values = 1:4;

cfg.fs = 256;
cfg.fres = 256;

% Sleep-adapted non-overlapping bands.
% Sigma is kept separate from beta; gamma is intentionally excluded for the
% DK5 LRT30 branch to match the corrected GC analysis.
cfg.band_defs = struct( ...
    'name', {'delta','theta','alpha','sigma','beta'}, ...
    'lo',   {1,5,9,13,17}, ...
    'hi',   {4,8,12,16,30});

cfg.iseed = 0;
cfg.nrunsp = env_num('SSDI_NRUNSP', 100);
cfg.nitersp = env_num('SSDI_NITERSP', 10000);
cfg.sig0p = 1;
cfg.gdesp = 2;
cfg.gdlsp = 2;
cfg.gdtolp = 1e-10;
cfg.histp = true;
cfg.ppp = true;

cfg.ctol = 1e-10;
cfg.niterso = env_num('SSDI_NITERSO', 10000);
cfg.gdeso = 2;
cfg.sig0o = 0.1;
cfg.gdlso = 2;
cfg.gdtolo = 1e-10;
cfg.histo = true;
cfg.ppo = true;

cfg.alpha = 0.05;
cfg.habeta_tails = 'both';
cfg.habeta_mhtc = true;

cfg.max_parpool_workers = env_num('SSDI_MAX_WORKERS', 12);
cfg.strict_n7 = true;
cfg.save_v73 = false;

cfg.parallel = struct();
cfg.parallel.enable = false;
cfg.parallel.n_workers = []; %#ok<NASGU>

cfg.recompute = true;
cfg.max_files = inf;

end

function v = env_num(name, fallback)
raw = getenv(name);
if isempty(raw)
    v = fallback;
    return;
end
x = str2double(raw);
if ~isfinite(x)
    v = fallback;
else
    v = x;
end
end
