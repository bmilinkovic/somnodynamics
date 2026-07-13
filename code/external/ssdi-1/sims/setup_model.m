%% Set up SS or VAR model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Define whether using a VAR model or fully-connected SS model
CON = tnet9x; % Connectivity matrix or scalar dimension
varmod = exist('CON', 'var'); 

if varmod
    % VAR model setup
    if isscalar(CON)
        n = CON;
        CON = ones(n);
    else
        n = size(CON, 1);
    end
    r = 7; % VAR model order
    w = 1; % VAR coefficients decay parameter
else
    % Fully-connected state-space model setup
    n = 9; % Microscopic dimension
    r = 3 * n; % Hidden state dimension
end

% Simulation and file parameters
moddir = 'results/nes_optimisation/';
modname = 'setup_model';

% Model parameters
rho = 0.9; % Spectral norm (< 1)
rmii = 1; % Residuals multi-information (0 for zero correlation)
fres = []; % Frequency resolution (automatic if empty)
mseed = 0; % Model random seed (0 uses current RNG state)

%% Generate random VAR or SS model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(mseed); % Seed random number generator
V0 = corr_rand(n, rmii); % Residuals covariance matrix

if varmod
    % Generate random VAR model
    ARA0 = var_rand(CON, r, rho, w);
    gc = var_to_pwcgc(ARA0, V0); % Causal graph
    [A0, C0, K0] = var_to_ss(ARA0); % Equivalent ISS model
    if isempty(fres)
        [fres, ~] = var2fres(ARA0, V0);
    end
    modcomp = r * n^2 + (n * (n + 1)) / 2; % Model complexity
    mdescript = sprintf('%d-variable VAR(%d)', n, r);
else
    % Generate random SS model
    [A0, C0, K0] = iss_rand(n, r, rho);
    gc = ss_to_pwcgc(A0, C0, K0, V0); % Causal graph
    if isempty(fres)
        [fres, ~] = ss2fres(A0, C0, K0, V0);
    end
    modcomp = (2 * n + 1) * r + (n * (n + 1)) / 2; % Model complexity
    mdescript = sprintf('%d-variable ISS(%d)', n, r);
end

% Model info
fprintf('\n---------------------------------------\n');
fprintf('Model            : %s\n', mdescript);
fprintf('---------------------------------------\n');
fprintf('Dimension        : %d\n', n);
fprintf('Complexity       : %d\n', modcomp);
fprintf('Freq. resolution : %d\n', fres);
fprintf('---------------------------------------\n\n');

% Normalize causal graph edgeweights
eweight = gc / nanmax(gc(:));

%% Save model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modfile = fullfile(moddir, [modname '.mat']);
fprintf('*** Saving model in ''%s''... ', modfile);

save_vars = {'mdescript', 'varmod', 'n', 'r', 'rho', 'rmii', 'V0', 'A0', 'C0', 'K0', 'gc', 'fres', 'modcomp', 'eweight'};
if varmod
    save_vars = [save_vars, {'ARA0'}];
end

save(modfile, save_vars{:});
fprintf('done\n\n');
