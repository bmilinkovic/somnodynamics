function results_root = sd_external_results_root(repo_root)
%SD_EXTERNAL_RESULTS_ROOT Central SomnoDynamics external results directory.
%
% Defaults to /Volumes/ex_data/data_somnoemergence/results. Override with
% SOMNODYNAMICS_DATA_ROOT for a different external-drive mount or machine.

if nargin < 1 || isempty(repo_root)
    this_file = mfilename('fullpath');
    repo_root = fileparts(fileparts(fileparts(this_file))); %#ok<NASGU>
end

data_root = getenv('SOMNODYNAMICS_DATA_ROOT');
if isempty(data_root)
    data_root = '/Volumes/ex_data/data_somnoemergence';
end

results_root = fullfile(data_root, 'results');
end
