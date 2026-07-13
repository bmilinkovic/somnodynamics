function f = var_to_siogc(A, V, inout, fres, H)
% In/out spectral GCs per variable
%
% This function computes the Granger causality (GC) in the frequency domain,
% either **to** ('in') or **from** ('out') a specific variable,
% with respect to the rest of the system.
%
% Arguments:
%   A:      VAR coefficient matrix (n x n x p), where n is the number of variables
%           and p is the model order.
%   V:      Residual covariance matrix (n x n) from the VAR model.
%   inout:  Direction of GC ('in' for GC to the variable, 'out' for GC from the variable).
%   fres:   Frequency resolution (number of frequency points).
%   H:      (Optional) Transfer function H(f) in the frequency domain. If not provided, it will be computed.
%
% Output:
%   f:      Spectral Granger causality matrix of size (n x h), where:
%           - n is the number of variables.
%           - h is the number of frequency points (fres + 1).
%           Each entry f(x, k) or f(y, k) represents the Granger causality value
%           at frequency point k, depending on the chosen 'in' or 'out' direction:
%               - For 'in', f(x, k) is the Granger causality **to** variable x from the rest of the system at frequency k.
%               - For 'out', f(y, k) is the Granger causality **from** variable y to the rest of the system at frequency k.
%
% Notes:
%   - Spectral Granger causality is calculated based on the transfer function H(f), 
%     which can either be provided as an argument or computed from the VAR coefficients A.
%   - The Granger causality is computed as the difference between the log power spectral
%     density (SR) and the log partial spectral density (SR - HR * HR').
%   - The transfer function (H) encapsulates the system's frequency response based on the VAR model.

% Step 1: Size validation and initialization
[n, n1, p] = size(A);  % Get dimensions of the VAR coefficients matrix
assert(n1 == n, 'VAR coefficients matrix has incorrect shape');  % Ensure A is square
[n1, n2] = size(V);  % Get dimensions of the residual covariance matrix V
assert(n1 == n && n2 == n, 'Residual covariance matrix must be square and match coefficient matrix');  % Ensure V matches A's dimensions

% Step 2: Determine the direction of GC (in/out)
% 'in' calculates GC to the variable from the rest of the system.
% 'out' calculates GC from the variable to the rest of the system.
switch lower(inout)
    case 'in'
        gcin = true;
    case 'out'
        gcin = false;
    otherwise
        error('The in/out parameter must be either ''in'' or ''out''');
end

% Step 3: Calculate the number of frequency bins and initialize the output matrix
h = fres + 1;  % Number of frequency points (frequency resolution + 1)
f = nan(n, h);  % Initialize the output matrix f to store GC values for each variable and frequency

% Step 4: Compute the transfer function H(f) IF it is not provided as an argument
% The transfer function H(f) relates the system's frequency response to the
% VAR coefficients. If not provided, it will be computed using var2trfun.
if nargin < 5 || isempty(H)
    H = var2trfun(A, fres);  % H is computed over all frequencies using the VAR coefficients
end

% Step 5: Cholesky decomposition of the residual covariance matrix V
% V = LL' where L is a lower triangular matrix
VL = chol(V, 'lower');  % Cholesky decomposition of V (Lower triangular)

% Step 6: Spectral Granger Causality Calculation (In or Out)
if gcin
    % Case 1: GC **to** each variable (from the rest of the system)
    for x = 1:n
        y = 1:n; y(x) = [];  % Exclude the target variable x
        PVL = chol(parcov(V, y, x), 'lower');  % Partial covariance for the remaining system

        for k = 1:h
            % Equation: H(f) * VL gives the frequency domain representation of the system
            HVL = H(x, :, k) * VL;  % Transfer function times Cholesky decomposition at frequency k
            SR = HVL * HVL';  % Power spectral density for variable x

            % Partial Granger causality (inbound)
            HR = H(x, y, k) * PVL;  % Transfer function excluding variable x
            f(x, k) = log(SR) - log(SR - HR * HR');  % Spectral GC equation for single variable
        end
    end
else
    % Case 2: GC **from** each variable (to the rest of the system)
    for y = 1:n
        x = 1:n; x(y) = [];  % Exclude the source variable y
        PVL = chol(parcov(V, y, x), 'lower');  % Partial covariance for the remaining system

        for k = 1:h
            HVL = H(x, :, k) * VL;  % Transfer function times Cholesky decomposition at frequency k
            SR = HVL * HVL';  % Power spectral density for the remaining system

            % Partial Granger causality (outbound)
            HR = H(x, y, k) * PVL;  % Transfer function excluding variable y
            f(y, k) = logdet(SR) - logdet(SR - HR * HR');  % Spectral GC equation for outbound GC
        end
    end
end
end
