function [dopt, Lopt, conv, iopt, sopt, cput, ohist] = opt_nes_dds_mruns(H, L0, niters, eta, tol, hist, pp)
% OPT_NES_DDS_MRUNS Run NES optimization across multiple runs
%
% Inputs:
%   H         - Transfer function
%   L0        - Initial subspace matrix (3D array, third dimension is number of runs)
%   niters    - Number of iterations
%   eta       - Learning rate
%   tol       - Convergence tolerance
%   hist      - If true, stores optimization history
%   pp        - If true, runs in parallel mode
%
% Outputs:
%   dopt      - Optimal dynamical dependence values (sorted)
%   Lopt      - Optimal subspace matrices (sorted)
%   conv      - Convergence status for each run
%   iopt      - Number of iterations for each run
%   sopt      - Placeholder (can store any additional result if required)
%   cput      - CPU time for each run
%   ohist     - Optimization history for each run (if hist is true)

% Borjan Milinkovic 2024

if nargin < 7 || isempty(pp), pp = false; end

nruns = size(L0, 3);

% Initialize outputs
dopt = zeros(1, nruns);
Lopt = zeros(size(L0));
conv = zeros(1, nruns);
iopt = zeros(1, nruns);
sopt = zeros(1, nruns); % Placeholder
cput = zeros(1, nruns);
if hist
    ohist = cell(nruns, 1);
else
    ohist = [];
end

% NES optimization across multiple runs
if pp
    parfor k = 1:nruns
        fprintf('NES optimization parallel run %4d of %4d : ', k, nruns);
        tcpu = cputime;
        [dopt(k), Lopt(:, :, k), conv(k), iopt(k), ohist{k}] = opt_nes_dds(H, L0(:, :, k), niters, eta, tol, hist);
        cput(k) = cputime - tcpu;
        fprintf('dd = %.4e : ', dopt(k));
        if conv(k) > 0
            fprintf('converged(%d)', conv(k));
        else
            fprintf('unconverged ');
        end
        fprintf(' in %4d iterations : CPU secs = %6.2f\n', iopt(k), cput(k));
    end
else
    for k = 1:nruns
        fprintf('NES optimization serial run %4d of %4d : ', k, nruns);
        tcpu = cputime;
        [dopt(k), Lopt(:, :, k), conv(k), iopt(k), ohist{k}] = opt_nes_dds(H, L0(:, :, k), niters, eta, tol, hist);
        cput(k) = cputime - tcpu;
        fprintf('dd = %.4e : ', dopt(k));
        if conv(k) > 0
            fprintf('converged(%d)', conv(k));
        else
            fprintf('unconverged ');
        end
        fprintf(' in %4d iterations : CPU secs = %6.2f\n', iopt(k), cput(k));
    end
end

% Sort everything by dynamical dependence
[dopt, sidx] = sort(dopt);
Lopt = Lopt(:, :, sidx);
iopt = iopt(sidx);
conv = conv(sidx);
sopt = sopt(sidx);
cput = cput(sidx);
if hist
    ohist = ohist(sidx);
end

end


