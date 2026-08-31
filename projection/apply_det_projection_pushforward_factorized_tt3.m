function C = apply_det_projection_pushforward_factorized_tt3(C, Bmat, ops, ndof, dim)
% APPLY_DET_PROJECTION_PUSHFORWARD_FACTORIZED_TT3
%
% Applies the determinant coefficient transfer and the weighted projection
% mass integration to one physical mode of a 3D TT tensor without assembling
%
%   Lproj = Bmat * F.
%
% The selected TT mode is assumed to contain grouped raw determinant
% coefficients with size n^3 and ordering
%
%   sub2ind([n,n,n], i, j, k),
%
% where i and j are the two non-derivative geometry-spline factors and k is
% the derivative factor.  For each TT-rank slice X(i,j,k), this evaluates
%
%   Y(a,k)  = sum_{i,j} Pcoef(a,(i,j)) X(i,j,k),
%   Z(a,mu) = sum_k D(mu,k) Y(a,k),
%   g(ell)  = sum_{(a,mu)} Rcoef(ell,(a,mu)) Z(a,mu),
%   h(q)    = sum_ell Bmat(q,ell) g(ell).
%
% Hence this function is algebraically identical to applying Bmat*F, but it
% keeps the transfer in its natural factorized form.


    G = C{dim};
    [r_left, nraw, r_right] = size(G);

    n = size(ops.D, 2);
    if nraw ~= n^3
        error(['apply_det_projection_pushforward_factorized_tt3: size mismatch. ', 'TT mode has size %d, but expected n^3 = %d.'], nraw, n^3);
    end

    n2    = size(ops.Pcoef, 1);
    nd    = size(ops.D, 1);
    namp  = size(ops.AM, 1);
    nout  = size(ops.Rcoef, 1);
    m_out = size(Bmat, 1);

    if size(ops.Rcoef, 2) ~= namp
        error('apply_det_projection_pushforward_factorized_tt3: Rcoef and AM are inconsistent.');
    end
    if size(Bmat, 2) ~= nout
        error('apply_det_projection_pushforward_factorized_tt3: Bmat and Rcoef are inconsistent.');
    end

    nslices = r_left * r_right;

    % Physical index first: n^3 x (r_left*r_right)
    Gmat = reshape(permute(G, [2 1 3]), nraw, nslices);

    % Interpret the grouped physical index as (i,j,k).
    X = reshape(Gmat, [n, n, n, nslices]);

    % 1) Pair-product contraction in the two non-derivative modes.
    Y = zeros(n2, n, nslices);

    for q = 1:size(ops.IJ, 1)
        i = ops.pair_i(q);
        j = ops.pair_j(q);

        Xij = reshape(X(i,j,:,:), n, nslices);
        if i ~= j
            Xij = Xij + reshape(X(j,i,:,:), n, nslices);
        end

        pcol = ops.Pcoef(:,q);
        [a_idx, ~, pval] = find(pcol);

        if isempty(a_idx)
            continue;
        end

        Y(a_idx,:,:) = Y(a_idx,:,:) + reshape(full(pval), [], 1, 1) .* reshape(Xij, 1, n, nslices);
    end

    % 2) Derivative contraction in the k-mode.
    Ymat = reshape(permute(Y, [2 1 3]), n, n2*nslices);
    Zmat = ops.D * Ymat;                         % nd x (n2*nslices)
    Z    = permute(reshape(Zmat, [nd, n2, nslices]), [2 1 3]);

    % 3) Mixed product contraction and 4) projection mass integration.
    Zam_full = reshape(Z, n2*nd, nslices);
    am_lin   = sub2ind([n2, nd], ops.am_a, ops.am_mu);
    Zam      = Zam_full(am_lin, :);              % namp x nslices

    Gdet = ops.Rcoef * Zam;                      % nout x nslices
    Hmat = Bmat * Gdet;                          % m_out x nslices

    H = reshape(Hmat, [m_out, r_left, r_right]);
    H = permute(H, [2 1 3]);                     % r_left x m_out x r_right

    C{dim} = reshape(H, [size(H, 1), ndof, ndof, size(H, 3)]);

end