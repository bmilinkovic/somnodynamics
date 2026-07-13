function F = var_to_gwcgc1(CAK, V, groups)
    % var_to_gwcgc1 computes Group-Wise Conditional Granger-Geweke Causality
    %
    % Inputs:
    %   CAK    - VAR coefficients matrix (n x n x p)
    %   V      - Residual covariance matrix (n x n)
    %   groups - Cell array of group indices
    %
    % Output:
    %   F      - Granger causality matrix (g x g), where F(a,b) represents
    %            the causality from group a to group b conditioned on others

    % Validate input dimensions
    [n, n1, p] = size(CAK);
    assert(n1 == n, 'VAR coefficients matrix has bad shape.');

    [n1_check, n2_check] = size(V);
    assert(n1_check == n && n2_check == n, 'Residuals covariance matrix must be square and match VAR coefficients.');

    % Validate and get number of groups
    g = check_group(groups, n);

    % Initialize causality matrix as g x g
    F = nan(g, g);

    for b = 1:g
        y = groups{b};
        r = 1:n; 
        r(y) = []; % Omit target group 'b'

        % Handle case where omitting group 'b' leaves no variables
        if isempty(r)
            warning('After omitting target group %d, no variables left to condition on. Assigning NaNs to causality measures.', b);
            F(:, b) = NaN;
            continue;
        end

        % Obtain reduced residual covariance matrix
        [~, VR, rep] = vardare(CAK, V, y, r);
        if sserror(rep, b)
            warning('DARE failed for target group %d. Assigning NaNs to causality measures.', b);
            F(:, b) = NaN;
            continue;
        end

        for a = 1:g
            if a == b
                continue; % Skip self-causality
            end

            x = groups{a};

            % Check if source group 'a' is empty
            if isempty(x)
                warning('Source group %d is empty. Assigning NaN to F(%d, %d).', a, a, b);
                F(a, b) = NaN;
                continue;
            end

            xr = findin(x, r); % Find indices of group 'a' within reduced variables 'r'

            % Check if all elements of 'x' are present in 'r'
            if any(isnan(xr))
                warning('Some elements of source group %d are not in conditioning set for target group %d. Assigning NaN to F(%d, %d).', a, b, a, b);
                F(a, b) = NaN;
                continue;
            end

            % Check for empty submatrices
            if isempty(xr)
                warning('After indexing, no variables found for source group %d in target group %d. Assigning NaN to F(%d, %d).', a, b, a, b);
                F(a, b) = NaN;
                continue;
            end

            try
                % Compute log determinants
                logdet_VRxr = logdet(VR(xr, xr));
                logdet_Vxx = logdet(V(x, x));

                % Assign causality statistic
                F(a, b) = logdet_VRxr - logdet_Vxx;
            catch ME
                warning('Error computing log determinants for F(%d, %d): %s. Assigning NaN.', a, b, ME.message);
                F(a, b) = NaN;
            end
        end
    end
end
