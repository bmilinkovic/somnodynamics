% function [gnum,gsiz,idx] = check_group(group,n)
% 
% assert(iscell(group),                     'groups must be specified as a cell vector of indices');
% idx = cell2mat(group(:)')';
% assert(length(unique(idx)) == length(idx),'groups indices must be unique and non-overlapping');
% if nargin > 1 && ~isempty(n);
% 	assert(all(idx >=1 & idx <= n),       'some group indices out of range');
% else
% 	assert(all(idx >=1),                  'some group indices out of range');
% end
% gnum = length(group);
% if nargout > 1
% 	gsiz = cellfun(@length,group);
% end

function [gnum, gsiz, idx] = check_group(group, n)
    % check_group validates group assignments
    %
    % Inputs:
    %   group - Cell array of group indices
    %   n     - Total number of variables (optional)
    %
    % Outputs:
    %   gnum  - Number of groups
    %   gsiz  - Sizes of each group
    %   idx   - Concatenated group indices
    
    % Ensure 'group' is a cell array
    assert(iscell(group), 'groups must be specified as a cell vector of indices');
    
    % Convert all group cells to column vectors
    group = cellfun(@(x) x(:), group, 'UniformOutput', false);
    
    % Concatenate all group indices into a single column vector
    idx = cell2mat(group(:));
    
    % Ensure groups are non-overlapping
    assert(length(unique(idx)) == length(idx), 'groups indices must be unique and non-overlapping');
    
    % If 'n' is provided, ensure indices are within range
    if nargin > 1 && ~isempty(n)
        assert(all(idx >= 1 & idx <= n), 'some group indices are out of range');
    else
        assert(all(idx >= 1), 'some group indices are out of range');
    end
    
    % Number of groups
    gnum = length(group);
    
    % Sizes of each group
    if nargout > 1
        gsiz = cellfun(@length, group);
    end
end

