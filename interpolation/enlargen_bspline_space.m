function [knots_det, degree_det, n_det] = enlargen_bspline_space(knots, degree, dim)
%ENLARGEN_BSPLINE_SPACE Construct the spline space containing det(grad G).
%
% Let the original univariate geometry space have degree p and let an
% interior knot have multiplicity mu. For a D-dimensional geometry, every
% term in det(grad G), considered in one parametric direction, has:
%
%   degree:                p_det  = D*p - 1,
%   interior multiplicity: mu_det = (D-1)*p + mu.
%
% The resulting open knot vector therefore contains det(grad G) exactly.
%
% Inputs:
%   knots   Open, nondecreasing knot vector.
%   degree  Original geometry degree p.
%   dim     Parametric dimension D.
%
% Outputs:
%   knots_det   Knot vector of the determinant space.
%   degree_det  Degree D*p - 1.
%   n_det       Number of basis functions.

    validateattributes(knots, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nondecreasing'});
    validateattributes(degree, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(dim, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    knots = knots(:).';

    left_endpoint = knots(1);
    right_endpoint = knots(end);

    [unique_knots, ~, knot_ids] = unique(knots);
    multiplicities = accumarray(knot_ids(:), 1).';

    degree_det = dim * degree - 1;

    % Open left endpoint: multiplicity degree_det + 1 = dim*degree.
    knots_det = repmat( ...
        left_endpoint, 1, degree_det + 1);

    % Interior knots.
    for i = 2:numel(unique_knots)-1
        multiplicity_det = ...
            (dim - 1) * degree + multiplicities(i);

        knots_det = [knots_det, ...
            repmat(unique_knots(i), ...
                   1, multiplicity_det)]; 
    end

    % Open right endpoint.
    knots_det = [knots_det, ...
        repmat(right_endpoint, 1, degree_det + 1)];

    n_det = numel(knots_det) - degree_det - 1;
end