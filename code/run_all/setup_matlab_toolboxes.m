function setup_matlab_toolboxes(repo_root)
%SETUP_MATLAB_TOOLBOXES Try to initialize MVGC2 and SSDI toolboxes.

if nargin < 1 || isempty(repo_root)
    this_file = mfilename('fullpath');
    repo_root = fileparts(fileparts(fileparts(this_file)));
end

% MVGC2 ---------------------------------------------------------------
mvgc_candidates = {
    getenv('MVGC2_DIR'), ...
    fullfile(repo_root, 'code', 'external', 'MVGC2'), ...
    fullfile(repo_root, 'MVGC2'), ...
    '/Users/borjan/code/dmt_project/MVGC2'
};

if exist('tsdata_to_varmo', 'file') ~= 2
    for i = 1:numel(mvgc_candidates)
        c = mvgc_candidates{i};
        if isempty(c) || exist(c, 'dir') ~= 7
            continue;
        end
        startup_file = fullfile(c, 'startup.m');
        if exist(startup_file, 'file') == 2
            old_pwd = pwd;
            cleanup = onCleanup(@() cd(old_pwd)); %#ok<NASGU>
            cd(c);
            run(startup_file);
            if exist('tsdata_to_varmo', 'file') == 2
                fprintf('[SomnoDynamics] MVGC2 initialized from %s\n', c);
                break;
            end
        end
    end
end

% SSDI ----------------------------------------------------------------
ssdi_candidates = {
    getenv('SSDI_TOOLBOX_DIR'), ...
    fullfile(repo_root, 'code', 'external', 'ssdi-1'), ...
    fullfile(repo_root, 'ssdi-1'), ...
    '/Users/borjan/code/dmt_project/ssdi-1', ...
    '/Users/borjan/code/matlab/ssdi-1'
};

if exist('opt_gd_dds_mruns', 'file') ~= 2
    for i = 1:numel(ssdi_candidates)
        c = ssdi_candidates{i};
        if isempty(c) || exist(c, 'dir') ~= 7
            continue;
        end
        addpath(genpath(c));
        if exist('opt_gd_dds_mruns', 'file') == 2
            fprintf('[SomnoDynamics] SSDI toolbox initialized from %s\n', c);
            break;
        end
    end
end

end
