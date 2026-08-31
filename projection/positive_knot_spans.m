function spans = positive_knot_spans(U)
% Return all knot spans of positive length.

    z = unique(U);
    tol = 1e-14 * max(1, z(end)-z(1));

    spans = zeros(0,2);
    for r = 1:numel(z)-1
        if z(r+1) > z(r) + tol
            spans(end+1,:) = [z(r),z(r+1)]; 
        end
    end
end