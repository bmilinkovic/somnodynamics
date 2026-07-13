function F = var_to_gwcggc(A,V,groups)
% Groupwise-conditional Global Granger Causality (GGC)
%
% This function computes the global Granger causality for each group of variables,
% conditioned on all other variables in the system.
% 
% Inputs:
% - A: VAR coefficient matrix (n x n x p), where n is the number of variables.
% - V: Residual covariance matrix (n x n) from the VAR model.
% - groups: Cell array containing groups of variable indices.
%
% Outputs:
% - F: A vector where each element corresponds to the global GC for each group.

% Step 1: Check dimensions of the VAR and residual covariance matrices
[n, n1, p] = size(A);  % Extract the size of VAR coefficient matrix
assert(n1 == n, 'VAR coefficients matrix has incorrect shape');  % Ensure the VAR matrix is square
[n1, n2] = size(V);  % Get the dimensions of the residual covariance matrix
assert(n1 == n && n2 == n, 'Residuals covariance matrix must be square and match coefficients matrix');  % Ensure V is square and matches A

% Step 2: Check the groups of variables
g = check_group(groups, n);  % Validate the group structure, return the number of groups

% Step 3: Initialize the output vector
F = nan(g, 1);  % Initialize the global GC output vector for each group

% Step 4: Extract the diagonal elements of V (residual variances)
DV = diag(V);  % Extract the diagonal of the residual covariance matrix V

% Step 5: Loop over each group to calculate Global GC
for a = 1:g
    % Get the indices of the current group
    x = groups{a};  % Indices for group 'a'
    z = 1:n; z(x) = [];  % Exclude group 'x' from the full set of variables

    nx = length(x);  % Number of variables in group 'x'
    VRx = zeros(nx, 1);  % Initialize the variance explained by the rest of the system

    % Step 6: Loop over each variable in the group
    for i = 1:nx
        y = x; y(i) = [];  % Remove variable 'x(i)' from the current group
        r = [x(i) z];  % Include the remaining group variables and other system variables

        % Step 7: Solve the Discrete Algebraic Riccati Equation (DARE) for reduced system
        [~, VR, rep] = vardare(A, V, y, r);  % Reduced system excluding variable x(i)
        
        % Step 8: Check if there are any errors in the DARE solver
        if sserror(rep, i)
            sserr = true;  % If there's an error, bail out of the loop
            break;
        end
        
        % Store the variance explained by the system for the current variable
        VRx(i) = VR(1, 1);  % Store the first element of VR for variable x(i)
    end

    % Step 9: Calculate the Global Granger Causality (GC) for group 'a'
    % Equation: Global GC for group a is the sum of the log-variances explained by the system
    % minus the sum of the log-diagonal elements of V for the group
    F(a) = sum(log(VRx)) - sum(log(DV(x)));  % Compute Global GC for group 'a'
end

