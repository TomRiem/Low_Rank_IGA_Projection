function [geometry, nurbs, G] = rotor_gen()
%ROTORBLATT2_GEOPDES Creates a GeoPDEs-compatible NURBS volume geometry.
%
% Usage:
%   [geometry, nurbs, G] = rotorblatt2_geopdes();
%
% Then in your assembly code:
%   [stiffnessMat, massMat] = stiffness_geoPDEs(nurbs, opt);

    % -----------------------------
    % Original rotor blade geometry
    % -----------------------------
    Ax = [0.8750, 2.0000, 0.9193, -0.2500, 0.3457, -0.1615, -1.3953, -0.8000, -1.3953, -2.5406, -2.3750, -2.7075, -3.7991, -3.8799, -4.2097, -4.6283, -4.9000, -4.6681] + 4.6681;

    Ay = [5.0000, 4.5000, 4.3500, 5.3000, 4.7875, 4.3000, 5.5908, 4.8834, 4.2000, 5.8816, 4.9662, 4.1780, 6.2631, 5.2590, 4.3184, 5.5000, 5.1225, 4.5000] - 4.5;

    Az = zeros(1,18);

    Bx = Ax;       By = Ay;       Bz = ones(1,18)*6;
    Cx = Ax;       Cy = Ay;       Cz = ones(1,18)*12;
    Dx = Ax*0.8;   Dy = Ay*0.8;   Dz = ones(1,18)*20;
    Ex = Ax*0.6;   Ey = Ay*0.6;   Ez = ones(1,18)*26;
    Fx = Ax*0.4;   Fy = Ay*0.4;   Fz = ones(1,18)*30;
    Gx = Ax*0.2;   Gy = Ay*0.2;   Gz = ones(1,18)*32;

    Hx = Ax*0.5;   Hy = Ay*0.5;   Hz = ones(1,18)*(-3);
    Jx = Hx;       Jy = Hy;       Jz = ones(1,18)*(-6);

    Ax = 0.9*Ax;
    Ay = 0.9*Ay;

    X = [Jx, Hx, Ax, Bx, Cx, Dx, Ex, Fx, Gx];
    Y = [Jy, Hy, Ay, By, Cy, Dy, Ey, Fy, Gy];
    Z = [Jz, Hz, Az, Bz, Cz, Dz, Ez, Fz, Gz];

    l = numel(X)/18;

    if abs(l - round(l)) > eps
        error('The number of control points must be divisible by 18.');
    end

    l = round(l);

    % -----------------------------
    % Custom geometry structure
    % -----------------------------
    G = struct();
    G.nurbs = 0;
    G.degree = [2, 2, 2];
    G.n = [3, 6, l];
    G.dim = 3;

    G.knots = { [0, 0, 0, 1, 1, 1], [0, 0, 0, 1/4, 1/2, 3/4, 1, 1, 1], [0, 0, 0:1/(l-2):1, 1, 1] ...
    };

    G.controlPoints = reshape([X; Y; Z].', 3, 6, l, 3);
    G.weight = ones(3, 6, l);

    % Optional: use your own fitting routine if it exists
    if exist('fitNURBS', 'file') == 2
        G = fitNURBS(G, 0);
    end

    % ---------------------------------------------------------
    % Reverse one parametric direction to obtain positive det(J)
    % Same physical geometry, opposite orientation.
    % ---------------------------------------------------------
    G.controlPoints = flip(G.controlPoints, 1);
    G.weight = flip(G.weight, 1);
    G.knots{1} = 1 - fliplr(G.knots{1});

    % -----------------------------
    % Convert to GeoPDEs / NURBS Toolbox format
    % -----------------------------
    nurbs = customG_to_geopdes_nurbs(G);

    % GeoPDEs geometry object
    geometry = geo_load(nurbs);
end


function nurbs = customG_to_geopdes_nurbs(G)
%CUSTOMG_TO_GEOPDES_NURBS Converts custom G-geometry into GeoPDEs NURBS format.

    requiredFields = {'controlPoints', 'weight', 'knots', 'degree', 'n'};
    for k = 1:numel(requiredFields)
        if ~isfield(G, requiredFields{k})
            error('Missing field G.%s.', requiredFields{k});
        end
    end

    n1 = G.n(1);
    n2 = G.n(2);
    n3 = G.n(3);

    if ~isequal(size(G.controlPoints), [n1, n2, n3, 3])
        error('G.controlPoints must have size [%d %d %d 3].', n1, n2, n3);
    end

    if ~isequal(size(G.weight), [n1, n2, n3])
        error('G.weight must have size [%d %d %d].', n1, n2, n3);
    end

    order = G.degree + 1;

    for d = 1:3
        expectedKnotLength = G.n(d) + order(d);
        actualKnotLength = numel(G.knots{d});

        if actualKnotLength ~= expectedKnotLength
            error(['Invalid knot vector in direction %d. ', 'Expected length %d, got length %d.'], d, expectedKnotLength, actualKnotLength);
        end
    end

    w = G.weight;

    coefs = zeros(4, n1, n2, n3);

    % Important: homogeneous coordinates are x*w, y*w, z*w, w
    coefs(1, :, :, :) = G.controlPoints(:, :, :, 1) .* w;
    coefs(2, :, :, :) = G.controlPoints(:, :, :, 2) .* w;
    coefs(3, :, :, :) = G.controlPoints(:, :, :, 3) .* w;
    coefs(4, :, :, :) = w;

    nurbs = struct();
    nurbs.form = 'B-NURBS';
    nurbs.dim = 4;
    nurbs.number = G.n;
    nurbs.coefs = coefs;
    nurbs.order = order;
    nurbs.knots = G.knots;
end