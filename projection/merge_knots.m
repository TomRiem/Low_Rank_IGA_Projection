function K = merge_knots(varargin)
% Merge knot vectors and retain only sorted breakpoints.

    allK = [];
    for q = 1:nargin
        allK = [allK, varargin{q}]; 
    end
    K = unique(sort(allK));
end