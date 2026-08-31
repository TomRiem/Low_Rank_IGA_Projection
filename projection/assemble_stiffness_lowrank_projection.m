function [K_tt, stats] = assemble_stiffness_lowrank_projection(tol, geometry, space, n_quad, proj_ref, space_proj)
% ASSEMBLE_STIFFNESS_LOWRANK_TT_SEPARABLE_NEW_2
% Assemble the isogeometric stiffness matrix directly in TT format.
%
% The geometry-dependent metric is written as a numerator divided by
% det(grad G). The reciprocal determinant is obtained by an L2 projection in
% a tensor-product spline space. The metric numerators are transferred to
% separable spline spaces, and the final stiffness integrals are contracted
% directly with the numerator and denominator TT cores. Large global
% coefficient-transfer and quadrature tensors are not assembled.
%
% INPUT
%   tolTT truncation and solve tolerance.
%   geometry    Geometry structure containing the B-spline geometry map.
%   space       Tensor-product solution space.
%   n_quad      Number of quadrature points in each parametric direction.
%   proj_ref    Projection-space selector used in the numerical experiments.
%      0/empty: build the projection space from the geometry space.
%      1:       build it from space_proj.
%      2:       use space_proj directly.
%   space_proj  Optional spline space used when proj_ref is 1 or 2.
%
% OUTPUT
%   K_tt        Stiffness matrix in TT-matrix format.
%   stats       Assembly time and storage statistics used in the experiments.
%
% SHARED DEPENDENCIES
% The following helpers are also used by the mass assembly and are therefore
% intentionally not duplicated in this file:
%   determinant_coefficients_tt, build_det_pushforward_ops,
%   build_derivative_operator, build_product_space_knots,
%   build_product_operator, assemble_univariate_mass,
%   tt_tensor_to_tt_matrix, positive_knot_spans, merge_knots,
%   gauss_legendre, eval_basis_on_span, find_span, basis_funs.
% The existing build_sixfold_space_only helper is also used unchanged.

    time = tic;

    tol_coeff  = 1e-13;
    tol_A      = 1e-12;
    tol_solve  = tol;
    tol_denom  = tol / 10;
    tol_num    = tol / 10;
    tol_apply  = tol / 10;
    tol_sum    = tol / 10;
    tol_target = tol;

    % ------------------------------------------------------------------
    % Reciprocal determinant projection.
    % ------------------------------------------------------------------
    det_ops   = cell(3,1);
    proj_space = cell(3,1);
    Bmat      = cell(3,1);
    rhs       = cell(3,1);

    for dim = 1:3
        U_geo = geometry.nurbs.knots{dim};
        p_geo = geometry.nurbs.order(dim) - 1;

        det_ops{dim} = build_det_pushforward_ops(U_geo, p_geo);
        proj_space{dim} = select_projection_space(U_geo, p_geo, proj_ref, space_proj, dim);

        [Bmat{dim}, rhs{dim}] = assemble_projection_Bmat_rhs_merged(proj_space{dim}.U, proj_space{dim}.p, det_ops{dim}.U3, det_ops{dim}.p_new);
    end

    % Exact determinant coefficient tensor and the three geometry-coordinate
    % coefficient tensors used later for the metric numerator.
    [A_det_tensor_tt, C_1, C_2, C_3] = determinant_coefficients_tt(tol_coeff, geometry);

    % Keep the raw determinant tensor for the post-timing statistics.
    A_det_raw_tensor_tt = A_det_tensor_tt;

    % Assemble the weighted projection matrix
    %
    %   A_ij = int phi_i phi_j det(grad G)
    %
    % directly in TT-matrix format. The determinant coefficient transfer and
    % the one-dimensional projection integrals are applied in factorized form.
    A_det_cores = core2cell(A_det_tensor_tt);

    for dim = 1:3
        A_det_cores = apply_det_projection_pushforward_factorized_tt3(A_det_cores, Bmat{dim}, det_ops{dim}, proj_space{dim}.n, dim);
    end

    A_det_tt = cell2core(tt_matrix, A_det_cores);
    A_det_tt = 0.5 * (A_det_tt + A_det_tt');
    A_det_tt = round(A_det_tt, tol_A);

    b = tt_tensor(rhs);

    nproj = cellfun(@(s) s.n, proj_space);
    one_tt = tt_tensor({ ones(nproj(1),1), ones(nproj(2),1), ones(nproj(3),1)});

    Aone = round(A_det_tt * one_tt, tol_A);
    alpha = dot(one_tt, b) / dot(one_tt, Aone);
    x0 = alpha * one_tt;

    % L2-orthonormalize the projection basis before the AMEn solve. The
    % Cholesky route is used whenever the univariate mass matrix is
    % numerically positive definite; otherwise dependent directions are
    % filtered by an eigendecomposition.
    W = cell(3,1);

    for dim = 1:3
        Md = assemble_univariate_mass(proj_space{dim}.U, proj_space{dim}.p);
        Md = full((Md + Md.') / 2);

        [R, flag] = chol(Md);
        if flag == 0
            n = size(Md,1);
            W{dim} = R \ eye(n);
        else
            [V, lam] = eig(Md, 'vector');
            keep = lam > 1e-13 * max(lam);
            W{dim} = V(:,keep) ./ sqrt(lam(keep)).';
        end
    end

    A_w = tt_mode_congruence(A_det_tt, W);
    b_w = tt_mode_apply(b, W, true);

    [y, ~] = amen_solve2(A_w, b_w, tol_solve, 'x0', x0, 'nswp', 1000, 'kickrank', 2, 'resid_damp', 10, 'trunc_norm', 'residual', 'max_full_size', 100, 'verb', 0);

    Denom_tt = tt_mode_apply(y, W, false);
    Denom_tt = round(Denom_tt, tol_denom);

    % ------------------------------------------------------------------
    % Numerator coefficient transfers.
    % ------------------------------------------------------------------
    pattern_names = {'p12','p13','p23','m1','m2','m3'};
    kind_names    = {'BBBB','DDBB','DBBB'};

    % Cache each one-dimensional four-factor transfer once per kind and
    % parametric direction. Pattern-specific operators below only reference
    % these cached objects.
    ops_kind = struct();
    for kk = 1:numel(kind_names)
        ops_kind.(kind_names{kk}) = cell(3,1);
    end

    for dim = 1:3
        U_geo = geometry.nurbs.knots{dim};
        p_geo = geometry.nurbs.order(dim) - 1;

        for kk = 1:numel(kind_names)
   kind = kind_names{kk};
   ops_kind.(kind){dim} = build_fourfold_pushforward(U_geo, p_geo, kind);
        end
    end

    ops_num = struct();
    pattern_info = struct();

    for pp = 1:numel(pattern_names)
        pat = pattern_names{pp};
        info = stiffness_pattern_info(pat);
        pattern_info.(pat) = info;
        ops_num.(pat) = cell(3,1);

        for dim = 1:3
   ops_num.(pat){dim} = ops_kind.(info.kinds{dim}){dim};
        end
    end

    % ------------------------------------------------------------------
    % Direct stiffness assembly.
    % ------------------------------------------------------------------
    entry_names = {'Q11','Q12','Q13','Q22','Q23','Q33'};

    entry_q = struct();
    entry_q.Q11 = 1;
    entry_q.Q12 = [2 4];
    entry_q.Q13 = [3 7];
    entry_q.Q22 = 5;
    entry_q.Q23 = [6 8];
    entry_q.Q33 = 9;

    K_tt = [];

    for ee = 1:numel(entry_names)
        entry = entry_names{ee};

        % Each metric entry is expanded into the six possible derivative
        % patterns and transferred directly to its numerator spline space.
        Num_pat = build_Q_entry_num_pattern_tt(entry, C_1, C_2, C_3, ops_num, tol_num);

        for pp = 1:numel(pattern_names)
   pat = pattern_names{pp};

   if isempty(Num_pat.(pat))
       continue;
   end

   Num_tt = Num_pat.(pat);
   q_list = entry_q.(entry);

   for qq = 1:numel(q_list)
       q = q_list(qq);

       if ~ismember(q, pattern_info.(pat).needed_q)
  continue;
       end

       % Contract the denominator and numerator coefficient TTs
       % directly with the one-dimensional stiffness integrals.
       K_tensor_tt = apply_qintegrals_to_separated_coeffs_tt3(Denom_tt,Num_tt,space,proj_space,ops_num.(pat),q,n_quad,tol_apply);

       if isempty(K_tensor_tt)
  continue;
       end

       K_part_tt = tt_tensor_to_tt_matrix(K_tensor_tt,space.ndof_dir(:),space.ndof_dir(:));

       if isempty(K_tt)
  K_tt = K_part_tt;
       else
  K_tt = round(K_tt + K_part_tt, tol_sum);
       end
   end
        end
    end

    K_tt = round(K_tt, tol_target);

    assembly_time = toc(time);

    % Statistics are evaluated only after the measured assembly time.
    stats = collect_stiffness_minimal_stats(assembly_time, Bmat, det_ops, A_det_raw_tensor_tt, A_det_tensor_tt, A_det_tt, b, W, Denom_tt, {C_1, C_2, C_3}, ops_kind, kind_names, K_tt);
end