function F = estimate_fisher(L, H)
% ESTIMATE_FISHER Estimates the Fisher information matrix for NES
% CURRENTLY UNDER CONSTRUCTION
% 
% Inputs:
%   L - Orthonormal subspace basis (n x m)
%   H - Transfer function (n x n x h)
%
% Outputs:
%   F - Fisher information matrix ((n*m) x (n*m))

% Initialize Fisher matrix
[n, m] = size(L);
F = zeros(n * m, n * m);

% Number of frequency bands
h = size(H, 3);

% Loop over frequency bands
for k = 1:h
    Hk = H(:, :, k);      % Hk: n x n
    HLk = Hk' * L;        % HLk: n x m
    Qk = HLk' * HLk;      % Qk: m x m
    Q_inv = Qk \ eye(m);  % Q_inv: m x m

    % Compute log-derivative terms for Fisher matrix
    J = kron(Q_inv, Hk) - kron(Q_inv * (HLk' * HLk) * Q_inv, eye(n));
    F = F + J' * J;
end

% Average Fisher matrix across frequency bands
F = F / h;
end
