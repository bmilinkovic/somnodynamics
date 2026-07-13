function [dopt,Lopt,conv,iopt,sopt,cput,ohist] = canonical_opt_dds_mruns( ...
    H, freqs, band_idx, L0, niters, gdes, gdsig0, gdls, gdtol, hist, pp, stage_label)
% Iterative canonical band-DD optimisation across multiple restarts.
%
% Objective per run:
%   1) compute canonical DD at every frequency bin (Eq. 25/26 form),
%   2) average over selected band indices.
%
% This keeps the two-stage SSDI structure (pre/final) while avoiding the
% non-canonical "optimise on H_band directly" path.
%
% NOTE:
%   Previous canonical_opt_single used random tangent perturbations
%   (accept/reject stochastic search). This implementation now performs
%   true gradient descent updates analogous to legacy opt_gd1_dds:
%       Ltry = orthonormalise(L - sig*(G/g))
%   with adaptive line-search-like sigma control and convergence checks.

if nargin < 12 || isempty(stage_label), stage_label = 'optimisation'; end
if nargin < 11 || isempty(pp), pp = false; end
if nargin < 10 || isempty(hist), hist = false; end
if nargin < 9 || isempty(gdtol), gdtol = 1e-10; end
if nargin < 8 || isempty(gdls), gdls = 2; end
if nargin < 7 || isempty(gdsig0), gdsig0 = 1; end
if nargin < 6 || isempty(gdes), gdes = 2; end
if nargin < 5 || isempty(niters), niters = 10000; end

nruns = size(L0,3);
dopt = inf(1,nruns);
Lopt = zeros(size(L0));
conv = zeros(1,nruns);
iopt = zeros(1,nruns);
sopt = zeros(1,nruns);
cput = zeros(1,nruns);
if hist
    ohist = cell(nruns,1);
else
    ohist = [];
end
run_logs = cell(nruns,1);

if pp
    try
        parfor k = 1:nruns
            [dopt(k), Lopt(:,:,k), conv(k), sopt(k), iopt(k), cput(k), hk, run_logs{k}] = run_single_restart( ...
                H, freqs, band_idx, L0(:,:,k), niters, gdes, gdsig0, gdls, gdtol, hist, stage_label, k, nruns);
            if hist
                ohist{k} = hk;
            end
        end
    catch ME
        warning('canonical_opt_dds_mruns:parallelFallback', ...
            ['CANON(%s) parfor failed (%s). Falling back to serial restarts for this block. ' ...
             'This is often caused by worker abort/resource instability.'], ...
            stage_label, strtrim(ME.message));
        [dopt, Lopt, conv, iopt, sopt, cput, ohist, run_logs] = run_serial_restarts( ...
            H, freqs, band_idx, L0, niters, gdes, gdsig0, gdls, gdtol, hist, stage_label);
    end
    for k = 1:nruns
        if ~isempty(run_logs{k})
            fprintf('%s\n', run_logs{k});
        end
    end
else
    [dopt, Lopt, conv, iopt, sopt, cput, ohist, run_logs] = run_serial_restarts( ...
        H, freqs, band_idx, L0, niters, gdes, gdsig0, gdls, gdtol, hist, stage_label);
    for k = 1:nruns
        fprintf('%s\n', run_logs{k});
    end
end

% Keep runs sorted by ascending objective value.
[dopt, sidx] = sort(dopt);
Lopt = Lopt(:,:,sidx);
conv = conv(sidx);
iopt = iopt(sidx);
sopt = sopt(sidx);
cput = cput(sidx);
if hist
    ohist = ohist(sidx);
end
end

function [dopt, Lopt, conv, iopt, sopt, cput, ohist, run_logs] = run_serial_restarts( ...
    H, freqs, band_idx, L0, niters, gdes, gdsig0, gdls, gdtol, hist, stage_label)

nruns = size(L0,3);
dopt = inf(1,nruns);
Lopt = zeros(size(L0));
conv = zeros(1,nruns);
iopt = zeros(1,nruns);
sopt = zeros(1,nruns);
cput = zeros(1,nruns);
if hist
    ohist = cell(nruns,1);
else
    ohist = [];
end
run_logs = cell(nruns,1);

for k = 1:nruns
    [dopt(k), Lopt(:,:,k), conv(k), sopt(k), iopt(k), cput(k), hk, run_logs{k}] = run_single_restart( ...
        H, freqs, band_idx, L0(:,:,k), niters, gdes, gdsig0, gdls, gdtol, hist, stage_label, k, nruns);
    if hist
        ohist{k} = hk;
    end
end
end

function [dout, Lout, convout, sigout, iout, cpu_secs, hout, run_log] = run_single_restart( ...
    H, freqs, band_idx, Linit, niters, gdes, gdsig0, gdls, gdtol, hist, stage_label, k, nruns)

t0 = cputime;
[dout, Lout, convout, sigout, iout, hout] = canonical_opt_single( ...
    H, freqs, band_idx, Linit, niters, gdes, gdsig0, gdls, gdtol, hist);
cpu_secs = cputime - t0;
if convout > 0
    conv_label = sprintf('converged(%d)', convout);
else
    conv_label = 'unconverged';
end
run_log = sprintf( ...
    'CANON(%s) run %4d/%4d : dd = %.4e : sig = %.4e : %s in %4d iterations : CPU secs = %6.2f', ...
    stage_label, k, nruns, dout, sigout, conv_label, iout, cpu_secs);
end

function [dd,L,converged,sig,iters,dhist] = canonical_opt_single( ...
    H, freqs, band_idx, L0, maxiters, gdes, sig0, gdls, tol, hist)

% Legacy-compatible line-search factors.
if isscalar(gdls)
    ifac = gdls;
    nfac = 1/ifac;
else
    ifac = gdls(1);
    nfac = gdls(2);
end

% Legacy-compatible convergence tolerances:
% [step-size, objective, gradient-norm].
if isscalar(tol)
    stol = tol;
    dtol = tol;
    gtol = tol;
else
    stol = tol(1);
    dtol = tol(min(2, numel(tol)));
    gtol = tol(min(3, numel(tol)));
end

[L, ~] = orthonormalise(L0);
sig = sig0;
converged = 0;

[dd, G, g, ok] = canonical_band_objective_and_grad(L, H, freqs, band_idx);
if ~ok || ~isfinite(dd)
    dd = inf;
    G = zeros(size(L));
    g = inf;
end

if hist
    % Columns: [objective, sigma, grad_norm, accepted_step_flag].
    dhist = nan(maxiters, 4);
    dhist(1,:) = [dd, sig, g, 1];
else
    dhist = [];
end

if ~(gdes == 1 || gdes == 2)
    error('unknown GD strategy');
end

for iters = 2:maxiters
    accepted = 0;

    if isfinite(g) && g > 0
        [Ltry, ~] = orthonormalise(L - sig*(G/g)); % true gradient-descent step
        [ddtry, Gtry, gtry, oktry] = canonical_band_objective_and_grad(Ltry, H, freqs, band_idx);

        if oktry && isfinite(ddtry) && ddtry < dd
            L = Ltry;
            dd = ddtry;
            G = Gtry;
            g = gtry;
            sig = ifac * sig;
            accepted = 1;
        else
            sig = nfac * sig;
        end
    else
        % Invalid/zero gradient: contract step and keep state.
        sig = nfac * sig;
    end

    if hist
        dhist(iters,:) = [dd, sig, g, accepted];
    end

    % Convergence tests mirror legacy opt_gd1_dds behavior.
    if sig < stol
        converged = 1;
        break;
    elseif dd < dtol
        converged = 2;
        break;
    elseif g < gtol
        converged = 3;
        break;
    end
end

if ~exist('iters', 'var')
    iters = 1;
end

if hist
    dhist = dhist(1:iters,:);
end
end
