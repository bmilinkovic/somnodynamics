function run_dual_atlas_pipelines(varargin)
%RUN_DUAL_ATLAS_PIPELINES Run legacy and Yeo7 pipelines side-by-side.
%
% Default behavior:
%   - legacy GC pipeline: healthy_sleep_mvgc.m
%   - Yeo7 GC/MI pipeline: run_yeo7_pipeline.m

p = inputParser;
addParameter(p, 'run_legacy', true);
addParameter(p, 'run_yeo7', true);
addParameter(p, 'yeo7_source_dir', '/Users/borjan/code/healthy-sleep-project/code/results/source_recon/yeo7/');
addParameter(p, 'yeo7_pwcgc_dir', '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/pwcgc/');
addParameter(p, 'yeo7_stats_dir', '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/');
parse(p, varargin{:});
cfg = p.Results;

if cfg.run_legacy
    fprintf('Running legacy pipeline via healthy_sleep_mvgc.m ...\n');
    run('healthy_sleep_mvgc.m');
end

if cfg.run_yeo7
    fprintf('Running Yeo7 pipeline ...\n');
    run_yeo7_pipeline( ...
        'source_dir', cfg.yeo7_source_dir, ...
        'pwcgc_dir', cfg.yeo7_pwcgc_dir, ...
        'stats_dir', cfg.yeo7_stats_dir);
end
end
