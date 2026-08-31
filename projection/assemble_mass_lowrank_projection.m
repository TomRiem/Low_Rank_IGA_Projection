function [M, stats] = assemble_mass_lowrank_projection(tol, geometry, space, n_quad, rounding)
% ASSEMBLE_MASS_LOWRANK_TT_ALL
% Assemble the isogeometric mass matrix directly in TT format.
%
% The determinant of the geometry map is represented by its exact spline
% coefficient tensor in TT format. In each parametric direction, the
% coefficient transfer
%
%   beta_i beta_j beta'_k  ->  S_{3p-1}
%
% and the univariate mass integration are applied directly to the TT cores.
% The large combined transfer matrices are therefore never assembled.
%
% INPUT
%   tol       TT truncation tolerance.
%   geometry  Geometry structure containing the B-spline geometry map.
%   space     Tensor-product solution space.
%   n_quad    Number of quadrature points in each parametric direction.
%   rounding  If true, TT rounding is applied after every mode transfer.
%
% OUTPUT
%   M         Mass matrix in TT-matrix format.
%   stats     Assembly time and storage statistics used in the experiments.

    time = tic;

    % Raw determinant coefficient tensor. Its three modes correspond to the
    % three parametric directions; each mode stores the grouped coefficients
    % of beta_i beta_j beta'_k.
    [M_tensor_tt, ~, ~, ~] = determinant_coefficients_tt(tol, geometry);
    Det_raw_tt = M_tensor_tt;

    % One-dimensional determinant transfers and integration operators.
    Bmat    = cell(3,1);
    det_ops = cell(3,1);

    for dim = 1:3
        U_geo = geometry.nurbs.knots{dim};
        p_geo = geometry.nurbs.order(dim) - 1;

        det_ops{dim} = build_det_pushforward_ops(U_geo, p_geo);

        Bmat{dim} = assemble_mass_Bmat_merged(space.knots{dim}, space.degree(dim), det_ops{dim}.U3, det_ops{dim}.p_new, [], n_quad(dim));
    end

    % Apply the factorized one-dimensional transfers directly to the TT cores.
    do_round = (nargin > 4) && ~isempty(rounding) && (rounding == true);

    for dim = 1:3
        if do_round
            M_tensor_tt = apply_mass_pushforward_factorized_tt3(M_tensor_tt, Bmat{dim}, det_ops{dim}, dim, tol);
        else
            M_tensor_tt = apply_mass_pushforward_factorized_tt3(M_tensor_tt, Bmat{dim}, det_ops{dim}, dim);
        end
    end

    if do_round
        M_tensor_tt = round(M_tensor_tt, tol);
    end

    % Each TT mode now contains the flattened one-dimensional matrix entries.
    M = tt_tensor_to_tt_matrix(M_tensor_tt, space.ndof_dir(:), space.ndof_dir(:));

    assembly_time = toc(time);

    stats = collect_mass_stats(assembly_time, Det_raw_tt, Bmat, det_ops);
end