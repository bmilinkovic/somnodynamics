%% Load Model and Set Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

% Load model saved by `setup_model.m`
moddir = 'results/nes_optimisation/';
modname = 'setup_model';
modfile = fullfile(moddir, [modname '.mat']);
assert(isfile(modfile), 'Model file not found. Run `setup_model.m` first.');
load(modfile); % Load variables: A, C, K, V0, fres, varmod, etc.

% General Parameters
mdim = 2; % Macroscopic dimension
nrunsp = 100; % Pre-optimization runs
nitersp = 10000; % Pre-optimization iterations
gdsig0p = 1; % Initial step size for pre-optimization
gdlsp = 2; % Line search parameters for pre-optimization
gdtolp = 1e-8; % Convergence tolerance for pre-optimization
histp = true; % Save optimization history
ppp = false; % Parallelize pre-optimization runs

niterso = 10000; % Final optimization iterations (GD and NES)
gdsig0o = 0.1; % Initial step size for final optimization (GD)
gdlso = 2; % Line search parameters for final optimization (GD)
gdtolo = 1e-10; % Convergence tolerance for final optimization (GD)
eta_nes = 0.1; % Learning rate for NES
tol_nes = 1e-6; % Convergence tolerance for NES
histo = true; % Save optimization history
ppo = false; % Parallelize runs

poptdir = 'results/nes_optimisation/'; % Directory to save pre-optimization results
poptname = 'preopt_dd'; % Pre-optimization filename root
optdir = 'results/nes_optimisation/'; % Directory to save final optimization results
optname = 'opt_dd'; % Final optimization filename root

%% Transform Model for Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if varmod
    [ARA, V] = transform_var(ARA0, V0); % Transform VAR model
    [A, C, K] = var_to_ss(ARA); % Equivalent ISS model
    CAK = ARA;
else
    [A, C, K, V] = transform_ss(A0, C0, K0, V0); % Transform SS model
    CAK = iss2cak(A, C, K);
end

fprintf('%s: pre-optimisation for m = %d\n\n', mdescript, mdim);

% Generate initial projections
rng(0); % Seed for reproducibility
L0p = rand_orthonormal(n, mdim, nrunsp);

%% Pre-Optimization (Gradient Descent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Starting pre-optimisation...\n');
tic;
[doptp, Lp, ~] = opt_gd_ddx_mruns(CAK, L0p, nitersp, 2, gdsig0p, gdlsp, gdtolp, histp, ppp);
toc;

% Transform pre-optimized results back for original residuals
Loptp = itransform_subspace(Lp, V);

fprintf('\nPre-optimisation complete. Optimal DD values:\n');
disp(doptp);

% Save pre-optimization results
poptfile = fullfile(poptdir, [poptname '_mdim_' num2str(mdim) '_H.mat']);
save(poptfile, 'doptp', 'Lp', 'Loptp');

%% Final Optimization (Gradient Descent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nStarting final optimisation using GD...\n');
% Cluster pre-optimized hyperplanes
[uidx, ~, ~] = Lcluster(Loptp, 1e-6, doptp);
L0o = Loptp(:, :, uidx); % Use clustered projections as initialization

% Calculate transfer function
if varmod
    H = var2trfun(ARA, fres);
else
    H = ss2trfun(A, C, K, fres);
end

% Perform final optimization using Gradient Descent
tic;
[dopto, Lo, ~] = opt_gd_dds_mruns(H, L0o, niterso, 2, gdsig0o, gdlso, gdtolo, histo, ppo);
toc;

% Transform final optimized results back for original residuals
Lopto = itransform_subspace(Lo, V);

fprintf('\nFinal optimisation complete. Optimal DD values:\n');
disp(dopto);

% Save final optimization results
optfile = fullfile(optdir, [optname '_mdim_' num2str(mdim) '_H.mat']);
save(optfile, 'dopto', 'Lo', 'Lopto');

%% Natural Gradient Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nStarting optimisation using Natural Gradient...\n');

% Generate random orthonormal initial projections for NES
rng(0); % Seed for reproducibility
L0o_nes = rand_orthonormal(n, mdim, 50); % 50 initializations

% Perform optimization using NES
tic;
[dopto_nes, Lo_nes, ~] = opt_nes_dds_mruns(H, L0o_nes, niterso, eta_nes, tol_nes, histo, ppo);
toc;

% Transform final optimized results back for original residuals
Lopto_nes = itransform_subspace(Lo_nes, V);

fprintf('\nNES optimisation complete. Optimal DD values:\n');
disp(dopto_nes);

% Save NES optimization results
nesoptfile = fullfile(optdir, ['opt_nes_dd_mdim_' num2str(mdim) '.mat']);
save(nesoptfile, 'dopto_nes', 'Lo_nes', 'Lopto_nes');

%% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot pre-optimization, GD, and NES results
figure;
plot(sort(doptp), '-o', 'LineWidth', 1.5, 'DisplayName', 'Pre-Optimization');
hold on;
plot(sort(dopto), '-x', 'LineWidth', 1.5, 'DisplayName', 'Final Optimization (GD)');
plot(sort(dopto_nes), '-s', 'LineWidth', 1.5, 'DisplayName', 'Final Optimization (NES)');
xlabel('Run Index (sorted)', 'FontSize', 12);
ylabel('Dynamical Dependence', 'FontSize', 12);
title('Optimization Comparison: Pre-Optimization, GD, NES', 'FontSize', 14);
legend('show', 'Location', 'Best');
grid on;
hold off;

% Save plot
plotFile = fullfile(optdir, 'optimization_comparison_gd_vs_nes.png');
saveas(gcf, plotFile);

fprintf('\nOptimization comparison plot saved to %s\n', plotFile);
