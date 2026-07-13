function [dd, L_opt, converged, iters, history] = opt_nes_dds(H, L0, maxiters, eta, tol, hist)
% OPT_NES_DDS Optimizes dynamical dependence using NES on Stiefel manifold
%
% Inputs:
%   H         - Transfer function
%   L0        - Initial orthonormal subspace (basis matrix)
%   maxiters  - Maximum number of iterations
%   eta       - Learning rate
%   tol       - Convergence tolerance
%   hist      - (Optional) If true, stores optimization history
%
% Outputs:
%   dd        - Optimal dynamical dependence value
%   L_opt     - Optimized subspace matrix
%   converged - Convergence status (1: converged, 0: max iterations reached)
%   iters     - Number of iterations performed
%   history   - (Optional) Optimization history (dynamical dependence and gradient norm)



% Borjan Milinkovic 2024

if nargin < 6
    hist = false;
end

% Initialization
L = L0;
converged = 0;

% For storing history if requested
if hist
    history.dd = zeros(maxiters, 1);
    history.grad_norm = zeros(maxiters, 1);
end

for iters = 1:maxiters
    % Calculate gradient and objective function
    [G, grad_magnitude] = trfun2ddgrad(L, H);
    dd = trfun2dd(L, H);

    % Estimate Fisher information matrix
    F = estimate_fisher(L, H);
    

    % Compute natural gradient
    grad_nat = F \ G;

    % Update L using the natural gradient
    L_new = L - eta * grad_nat;

    % Project back to Stiefel manifold using SVD
    [U, ~, V] = svd(L_new, 'econ');
    L = U * V';

    % Store history
    if hist
        history.dd(iters) = dd;
        history.grad_norm(iters) = norm(grad_nat, 'fro');
    end

    % Check for convergence
    if grad_magnitude < tol
        converged = 1;
        break;
    end
end

% Final results
L_opt = L;

% Trim history if necessary
if hist
    history.dd = history.dd(1:iters);
    history.grad_norm = history.grad_norm(1:iters);
end

end

