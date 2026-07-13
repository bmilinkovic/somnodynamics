function [d, G, mG, ok, info] = canonical_band_objective_and_grad(L, H, freqs, band_idx)
%CANONICAL_BAND_OBJECTIVE_AND_GRAD Canonical band-DD value and gradient.
%
% This computes the canonical spectral DD objective used by the canonical
% pipeline (Eq. 25/26 style per frequency, then band average), together with
% an analytic matrix-calculus gradient for gradient-descent optimisation.
%
% Inputs
%   L        : n x m orthonormal basis (real-valued in current pipeline)
%   H        : n x n x nf transfer function tensor
%   freqs    : nf-vector of frequencies corresponding to H(:,:,k)
%   band_idx : index vector of frequency bins included in band objective
%
% Outputs
%   d    : scalar band objective (Inf on invalid numerics)
%   G    : n x m Grassmann-projected gradient
%   mG   : Frobenius norm of G
%   ok   : true if objective/gradient valid
%   info : diagnostics struct (contains Euclidean gradient before projection)

if nargin < 4 || isempty(band_idx)
    d = inf;
    G = zeros(size(L));
    mG = inf;
    ok = false;
    info = struct('reason', 'empty_band_idx', 'G_euclid', G);
    return;
end

idx = band_idx(:)';
nf = size(H, 3);
if any(idx < 1) || any(idx > nf)
    d = inf;
    G = zeros(size(L));
    mG = inf;
    ok = false;
    info = struct('reason', 'band_idx_out_of_range', 'G_euclid', G);
    return;
end

fvals = nan(1, numel(idx));
Gstack = nan([size(L), numel(idx)]);
for ii = 1:numel(idx)
    [fvals(ii), Gfi, okfi] = canonical_dd_at_freq_and_grad(L, H(:,:,idx(ii)));
    if ~okfi || ~isfinite(fvals(ii)) || any(~isfinite(Gfi(:)))
        d = inf;
        G = zeros(size(L));
        mG = inf;
        ok = false;
        info = struct('reason', 'invalid_freq_term', 'freq_index', idx(ii), 'G_euclid', G);
        return;
    end
    Gstack(:,:,ii) = Gfi;
end

[d, weights, okavg] = interval_average_with_weights(fvals, freqs(idx));
if ~okavg || ~isfinite(d)
    d = inf;
    G = zeros(size(L));
    mG = inf;
    ok = false;
    info = struct('reason', 'invalid_band_average', 'G_euclid', G);
    return;
end

Ge = zeros(size(L));
for ii = 1:numel(idx)
    Ge = Ge + weights(ii) * Gstack(:,:,ii);
end

% Grassmann projection (matches legacy manifold-gradient intent).
G = Ge - L*(L'*Ge);
G = real(G);
mG = norm(G, 'fro');
ok = isfinite(mG) && all(isfinite(Ge(:)));

info = struct();
info.reason = '';
info.freq_idx = idx;
info.fvals = fvals;
info.weights = weights;
info.G_euclid = Ge;
end

function [fval, G, ok] = canonical_dd_at_freq_and_grad(L, Hk)
% Per-frequency canonical objective and analytic gradient.
%
% f(L) = log|A| - log|C|
% where
%   A = L' (H H*) L
%   B = L' H L
%   C = B B*
%
% Gradient is derived via matrix calculus on the same objective quantity.
% For real-valued optimisation variables L with complex H:
%   grad(log|A|)   = 2*Re(S*L*A^{-1}),  S = H H*
%   grad(log|C|)   = 2*Re(H*L*B^{-1} + H'*L*(B')^{-1})
% and we take their difference.

ok = false;
fval = inf;
G = zeros(size(L));

S = Hk*Hk';
A = L' * S * L;
A = (A + A')/2;
B = L' * Hk * L;
C = B * B';
C = (C + C')/2;

[fval, okf] = logdet_ratio(A, C);
if ~okf || ~isfinite(fval)
    return;
end

% A must be SPD for log|A| gradient term.
[LA, pA] = chol(A, 'lower');
if pA ~= 0
    return;
end
Ainv = LA' \ (LA \ eye(size(A), class(A)));

if rcond(B) < 1e-14
    return;
end
Binv = B \ eye(size(B), class(B));

GA = 2*real(S*L*Ainv);
GB = 2*real(Hk*L*Binv + Hk'*L*(Binv'));
G = GA - GB;

if any(~isfinite(G(:)))
    return;
end
ok = true;
end

function [avg, w, ok] = interval_average_with_weights(vals, freq_vec)
% Same averaging rule as pipeline objective, with explicit linear weights.

ok = false;
avg = nan;
w = nan(size(vals));

n = numel(vals);
if n == 0 || n ~= numel(freq_vec)
    return;
end

if n == 1
    avg = vals(1);
    w = 1;
    ok = isfinite(avg);
    return;
end

span = freq_vec(end) - freq_vec(1);
if span <= 0
    w = ones(1, n) / n;
    avg = sum(w .* vals);
    ok = isfinite(avg) && all(isfinite(w));
    return;
end

% trapz equivalent linear weights.
df = diff(freq_vec);
w = zeros(1, n);
w(1) = df(1)/2;
w(end) = df(end)/2;
for i = 2:n-1
    w(i) = (df(i-1) + df(i))/2;
end
w = w / span;
avg = sum(w .* vals);
ok = isfinite(avg) && all(isfinite(w));
end

function [v, ok] = logdet_ratio(A, C)
ok = false;
v = NaN;
A = (A + A')/2;
C = (C + C')/2;
[Lc, pc] = chol(C, 'lower');
if pc ~= 0
    return;
end
T = Lc \ A / Lc';
T = (T + T')/2;
e = eig(T);
e = real(e(:));
if any(~isfinite(e)) || any(e <= 0)
    return;
end
v = sum(log(e));
ok = isfinite(v);
end
