function run_ssdi_sleep_stage_pipeline(varargin)
%RUN_SSDI_SLEEP_STAGE_PIPELINE SomnoDynamics SSDI pipeline entrypoint.
%
% Default behavior runs the canonical stage-based SSDI pipeline:
%   code/ssdi/stage_ssdi_anatomical7/run_stage_ssdi_anatomical7_pipeline.m
%
% Set 'mode' to 'legacy' to run the old pwcgc-driven compatibility wrapper.

p = inputParser;
addParameter(p, 'mode', 'stage_ssdi');
addParameter(p, 'input_root', '');
addParameter(p, 'output_root', '');
addParameter(p, 'recompute', true);
addParameter(p, 'parallel_enable', true);
addParameter(p, 'parallel_workers', []);
addParameter(p, 'max_files', inf);
addParameter(p, 'nrunsp', 100);
addParameter(p, 'nitersp', 10000);
addParameter(p, 'niterso', 10000);
parse(p, varargin{:});
opt = p.Results;

this_file = mfilename('fullpath');
code_dir = fileparts(this_file);
repo_root = fileparts(fileparts(code_dir));

addpath(fullfile(repo_root, 'code', 'run_all'));
setup_matlab_toolboxes(repo_root);
results_root = sd_external_results_root(repo_root);

switch lower(string(opt.mode))
    case "stage_ssdi"
        stage_dir = fullfile(code_dir, 'stage_ssdi_anatomical7');
        addpath(stage_dir);

        call_args = {
            'recompute', opt.recompute, ...
            'parallel_enable', opt.parallel_enable, ...
            'parallel_workers', opt.parallel_workers, ...
            'max_files', opt.max_files, ...
            'nrunsp', opt.nrunsp, ...
            'nitersp', opt.nitersp, ...
            'niterso', opt.niterso ...
        };
        if ~isempty(opt.input_root)
            call_args = [call_args, {'input_root', opt.input_root}]; %#ok<AGROW>
        end
        if ~isempty(opt.output_root)
            call_args = [call_args, {'results_root', opt.output_root}]; %#ok<AGROW>
        end

        fprintf('\n[SomnoDynamics SSDI] mode=stage_ssdi\n');
        run_stage_ssdi_anatomical7_pipeline(call_args{:});

    case "legacy"
        legacy_out = opt.output_root;
        if isempty(legacy_out)
            legacy_out = fullfile(results_root, 'ssdi', 'legacy_wrapper');
        end

        legacy_in = opt.input_root;
        if isempty(legacy_in)
            legacy_in = fullfile(results_root, 'gc', 'pwcgc');
        end

        fprintf('\n[SomnoDynamics SSDI] mode=legacy compatibility wrapper\n');
        addpath(code_dir);
        run_ssdi_sleep_stage_pipeline_legacy_from_pwcgc( ...
            'pwcgc_dir', legacy_in, ...
            'out_root', legacy_out, ...
            'nrunsp', opt.nrunsp, ...
            'niters', opt.nitersp);

    otherwise
        error('Unsupported mode: %s', opt.mode);
end

end
