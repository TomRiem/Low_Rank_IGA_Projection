function ops = build_det_pushforward_ops(U, p)
% BUILD_DET_PUSHFORWARD_OPS
% Build the factorized coefficient transfer for
%
%   beta_i beta_j beta'_k -> S_{3p-1}(U3).
%
% The transfer is split into derivative, pair-product, and mixed-product
% operators. The full n_out-by-n^3 pushforward matrix is not formed.

    [D, Ud] = build_derivative_operator(U, p);

    U2 = build_product_space_knots(U, p, U, p);
    U3 = build_product_space_knots(U2, 2*p, Ud, p-1);

    [IJ, Pcoef] = build_product_operator(U, p, U, p, U2, 2*p, true);

    [AM, Rcoef] = build_product_operator(U2, 2*p, Ud, p-1, U3, 3*p-1, false);

    ops.U     = U;
    ops.p     = p;
    ops.Ud    = Ud;
    ops.U2    = U2;
    ops.U3    = U3;
    ops.p_new = 3*p - 1;
    ops.nout  = numel(U3) - ops.p_new - 1;

    ops.D     = D;
    ops.IJ    = IJ;
    ops.Pcoef = Pcoef;
    ops.AM    = AM;
    ops.Rcoef = Rcoef;

    % Cached indices used by the factorized TT contraction.
    ops.pair_i = IJ(:,1);
    ops.pair_j = IJ(:,2);

    ops.am_a  = AM(:,1);
    ops.am_mu = AM(:,2);
end