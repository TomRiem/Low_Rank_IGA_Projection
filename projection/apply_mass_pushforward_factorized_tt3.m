function T = apply_mass_pushforward_factorized_tt3(T, Bmat, ops, dim, tol)
% APPLY_MASS_PUSHFORWARD_FACTORIZED_TT3
% Apply the determinant coefficient transfer and mass integration directly
% to one TT core without assembling their combined matrix.

    if isempty(T)
        return;
    end

    C = core2cell(T);
    G = C{dim};

    [r_left, nraw, r_right] = size(G);

    n = size(ops.D,2);
    if nraw ~= n^3
        error(['apply_mass_pushforward_factorized_tt3: size mismatch. ', 'TT mode has size %d, expected %d.'], nraw, n^3);
    end

    n2       = size(ops.Pcoef,1);
    nd       = size(ops.D,1);
    namp     = size(ops.AM,1);
    nout     = size(ops.Rcoef,1);
    m_out    = size(Bmat,1);
    nslices  = r_left * r_right;

    if size(ops.Rcoef,2) ~= namp
        error('Rcoef and AM are inconsistent.');
    end
    if size(Bmat,2) ~= nout
        error('Bmat and Rcoef are inconsistent.');
    end

    % Raw grouped coefficients X(i,j,k) for all TT-rank slices.
    Gmat = reshape(permute(G,[2 1 3]), nraw, nslices);
    X = reshape(Gmat, [n,n,n,nslices]);

    % beta_i beta_j -> S_{2p}.
    Y = zeros(n2,n,nslices);

    for q = 1:size(ops.IJ,1)
        i = ops.pair_i(q);
        j = ops.pair_j(q);

        Xij = reshape(X(i,j,:,:), n, nslices);
        if i ~= j
            Xij = Xij + reshape(X(j,i,:,:), n, nslices);
        end

        [a_idx, ~, pval] = find(ops.Pcoef(:,q));
        if isempty(a_idx)
            continue;
        end

        Y(a_idx,:,:) = Y(a_idx,:,:) + reshape(full(pval),[],1,1) .* reshape(Xij,1,n,nslices);
    end

    % Derivative coefficient map in the third factor.
    Ymat = reshape(permute(Y,[2 1 3]), n, n2*nslices);
    Zmat = ops.D * Ymat;
    Z = permute(reshape(Zmat,[nd,n2,nslices]), [2 1 3]);

    % S_{2p} * S_{p-1} -> S_{3p-1}, followed by mass integration.
    Zam_full = reshape(Z, n2*nd, nslices);
    am_lin = sub2ind([n2,nd], ops.am_a, ops.am_mu);
    Zam = Zam_full(am_lin,:);

    Gdet = ops.Rcoef * Zam;
    Hmat = Bmat * Gdet;

    H = reshape(Hmat, [m_out,r_left,r_right]);
    C{dim} = permute(H,[2 1 3]);

    T = cell2core(tt_tensor, C);

    if nargin >= 5 && ~isempty(tol)
        T = round(T, tol);
    end
end