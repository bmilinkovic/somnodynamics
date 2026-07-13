function run_gc_sleep_stage_pipeline(varargin)
%RUN_GC_SLEEP_STAGE_PIPELINE SomnoDynamics GC pipeline entrypoint.
%
% Default behavior runs the new stage-based MVGC2 pipeline:
%   code/gc/stage_mvgc2_anatomical7/run_stage_mvgc2_anatomical7_pipeline.m
%
% Set 'mode' to 'legacy' to run the old compatibility wrapper.

p = inputParser;
addParameter(p, 'mode', 'stage_mvgc2');
addParameter(p, 'atlas', 'anatomical7');
addParameter(p, 'source_root', '');
addParameter(p, 'output_root', '');
addParameter(p, 'recompute_subject_measures', true);
addParameter(p, 'recompute_group_stats', true);
addParameter(p, 'parallel_enable', true);
addParameter(p, 'parallel_workers', []);
addParameter(p, 'max_files', inf);
addParameter(p, 'normalization_mode', 'pooled_zscore');
addParameter(p, 'normalization_eps', 1e-12);
addParameter(p, 'fs_original', 128);
addParameter(p, 'fs_target', 128);
addParameter(p, 'fres', 128);
addParameter(p, 'momax', []);
addParameter(p, 'mosel', '');
parse(p, varargin{:});
opt = p.Results;

this_file = mfilename('fullpath');
code_dir = fileparts(this_file);
repo_root = fileparts(fileparts(code_dir));

addpath(fullfile(repo_root, 'code', 'run_all'));
setup_matlab_toolboxes(repo_root);
results_root = sd_external_results_root(repo_root);

switch lower(string(opt.mode))
    case "stage_mvgc2"
        stage_dir = fullfile(code_dir, 'stage_mvgc2_anatomical7');
        addpath(stage_dir);

        if isempty(opt.source_root)
            opt.source_root = fullfile(results_root, 'source_recon', 'anatomical7');
        end
        if isempty(opt.output_root)
            opt.output_root = fullfile(results_root, 'gc', 'stage_mvgc2_anatomical7_256hz');
        end

        call_args = {
            'recompute_subject_measures', opt.recompute_subject_measures, ...
            'recompute_group_stats', opt.recompute_group_stats, ...
            'parallel_enable', opt.parallel_enable, ...
            'parallel_workers', opt.parallel_workers, ...
            'max_files', opt.max_files, ...
            'normalization_mode', opt.normalization_mode, ...
            'normalization_eps', opt.normalization_eps, ...
            'fs_original', opt.fs_original, ...
            'fs_target', opt.fs_target, ...
            'fres', opt.fres ...
        };

        if ~isempty(opt.source_root)
            call_args = [call_args, {'source_root', opt.source_root}]; %#ok<AGROW>
        end
        if ~isempty(opt.output_root)
            call_args = [call_args, {'results_root', opt.output_root}]; %#ok<AGROW>
        end
        if ~isempty(opt.momax)
            call_args = [call_args, {'momax', opt.momax}]; %#ok<AGROW>
        end
        if ~isempty(opt.mosel)
            call_args = [call_args, {'mosel', opt.mosel}]; %#ok<AGROW>
        end

        fprintf('\n[SomnoDynamics GC] mode=stage_mvgc2\n');
        fprintf('[SomnoDynamics GC] source root: %s\n', opt.source_root);
        fprintf('[SomnoDynamics GC] output root: %s\n', opt.output_root);
        run_stage_mvgc2_anatomical7_pipeline(call_args{:});

    case "legacy"
        fprintf('\n[SomnoDynamics GC] mode=legacy compatibility wrapper\n');
        run_legacy_gc_wrapper(opt, repo_root);

    otherwise
        error('Unsupported mode: %s', opt.mode);
end

end

function run_legacy_gc_wrapper(opt, repo_root)
legacy_root = fullfile(repo_root, 'code', 'legacy');
addpath(legacy_root);
addpath(fullfile(legacy_root, 'fun_conn_analysis'));

if exist('tsdata_to_varmo', 'file') ~= 2
    error(['MVGC toolbox functions are not on the MATLAB path. ' ...
        'Initialize MVGC before running this pipeline.']);
end

source_root = opt.source_root;
if isempty(source_root)
    source_root = fullfile(sd_external_results_root(repo_root), 'source_recon', 'anatomical7');
end

output_root = opt.output_root;
if isempty(output_root)
    output_root = fullfile(sd_external_results_root(repo_root), 'gc', 'legacy_wrapper');
end

pwcgc_dir = fullfile(output_root, 'pwcgc');
if exist(output_root, 'dir') ~= 7
    mkdir(output_root);
end
if exist(pwcgc_dir, 'dir') ~= 7
    mkdir(pwcgc_dir);
end

fprintf('[SomnoDynamics GC][legacy] source root: %s\n', source_root);
fprintf('[SomnoDynamics GC][legacy] pwcgc dir : %s\n', pwcgc_dir);
fprintf('[SomnoDynamics GC][legacy] output dir: %s\n', output_root);

run_mvgc_atlas('atlas', opt.atlas, 'src_data', source_root, 'pwcgc_dir', pwcgc_dir);
run_gc_mi_stats_atlas('atlas', opt.atlas, 'pwcgc_dir', pwcgc_dir, 'output_dir', output_root);

fprintf('[SomnoDynamics GC][legacy] complete.\n');
end
