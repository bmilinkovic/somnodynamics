global healthy_sleep;
healthy_sleep = fileparts(mfilename('fullpath'));
addpath(healthy_sleep);

subdirs = genpath(healthy_sleep);
addpath(subdirs);
fprintf('[Emergence Pipeline startup] Added path code directory: %s\n', healthy_sleep);

global ssdiDir;
ssdiDir = '/Users/borjan/code/matlab/ssdi-1'; % I just changed this to the new ssdi-1 toolbox might not work, so please test
addpath((ssdiDir));
cd(ssdiDir);
startup;
cd(healthy_sleep);
fprintf('[Emergence Pipeline startup] Added path to State-Space Dynamical Indpendence toolbox: %s\n',ssdiDir);

global ceDir;
ceDir = '/Users/borjan/code/di_eeg/code/ReconcilingEmergences';
addpath(genpath(ceDir));
fprintf('[Emergence Pipeline startup] Added path to Strong Causal Emergence toolbox: %s\n',ceDir);

%path_rawdata = '/Volumes/dataSets/restEEGHealthySubjects/rawData';
%path_preproc = '/Volumes/dataSets/restEEGHealthySubjects/preprocessedData';
%path_source = '/Volumes/dataSets/restEEGHealthySubjects/preprocessedData/sourceReconstructions';







