function run_yeo7_pipeline(varargin)
%RUN_YEO7_PIPELINE Run Yeo7 GC + MI downstream pipeline.
%
% Usage:
%   run_yeo7_pipeline
%   run_yeo7_pipeline('source_dir','/path/to/source_recon/yeo7')

p = inputParser;
addParameter(p, 'source_dir', '/Users/borjan/code/healthy-sleep-project/code/results/source_recon/yeo7/');
addParameter(p, 'pwcgc_dir', '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/pwcgc/');
addParameter(p, 'stats_dir', '/Users/borjan/code/healthy-sleep-project/code/results/yeo7/gc/');
addParameter(p, 'run_adf', true);
parse(p, varargin{:});
cfg = p.Results;

run_mvgc_atlas( ...
    'atlas', 'yeo7', ...
    'src_data', cfg.source_dir, ...
    'pwcgc_dir', cfg.pwcgc_dir, ...
    'run_adf', cfg.run_adf);

run_gc_mi_stats_atlas( ...
    'atlas', 'yeo7', ...
    'pwcgc_dir', cfg.pwcgc_dir, ...
    'output_dir', cfg.stats_dir);
end
