function T = apply_qintegrals_to_separated_coeffs_tt3(Denom_tt, Num_tt, space, spline_space_denom, ops_num_pat, q, nquad_dim, tol)
% APPLY_QINTEGRALS_TO_SEPARATED_COEFFS_TT3
%
% Contracts Denom_tt and Num_tt directly with the one-dimensional integrals
%
%   int (d^alpha B_i)(d^beta B_j) Phi_a^den Psi_b^num
%
% without assembling/storing the four-dimensional TT tensors Q_int{q}{dim}.
%
% The output is the same 3D TT tensor as apply_qtensors_to_separated_coeffs_tt3
% would produce from preassembled Q_int tensors, but the memory-intensive
% intermediate Q_int is avoided.

    if isempty(Denom_tt) || isempty(Num_tt)
        T = [];
        return;
    end

    if nargin < 7 || isempty(nquad_dim)
        nquad_dim = zeros(1,3);
        for dim = 1:3
            psol = space.degree(dim);
            pden = spline_space_denom{dim}.p;
            pnum = ops_num_pat{dim}.pout;
            nquad_dim(dim) = ceil((2*psol + pden + pnum + 1)/2);
        end
    end

    if nargin < 8 || isempty(tol)
        tol = 1e-10;
    end

    a_of_q = [1 2 3 1 2 3 1 2 3];
    b_of_q = [1 1 1 2 2 2 3 3 3];

    a = a_of_q(q);
    b = b_of_q(q);

    D = core2cell(Denom_tt);
    N = core2cell(Num_tt);

    if numel(D) ~= 3 || numel(N) ~= 3
        error('apply_qintegrals_to_separated_coeffs_tt3: expected three-dimensional TTs.');
    end

    C = cell(3,1);

    for dim = 1:3
        Usol = space.knots{dim};
        psol = space.degree(dim);
        nsol = space.ndof_dir(dim);

        Uden = spline_space_denom{dim}.U;
        pden = spline_space_denom{dim}.p;
        nden = spline_space_denom{dim}.n;

        Unum = ops_num_pat{dim}.Uout;
        pnum = ops_num_pat{dim}.pout;
        nnum = ops_num_pat{dim}.nout;

        spans = positive_knot_spans(merge_knots(Usol, Uden, Unum));
        nq    = nquad_dim(dim);

        C{dim} = contract_1d_integrals_with_factor_cores_direct(D{dim}, N{dim}, Usol, psol, nsol, Uden, pden, nden, Unum, pnum, nnum, spans, nq, dim, a, b);

        % Incremental compression: round the partial tensor after each new
        % dimension is appended, so the three full-rank cores never coexist.
        % The dense core for the current dim is freed before the next dim is
        % built. This is numerically identical to building all three cores
        % and rounding once at the end (round is associative over sweeps),
        % but caps peak memory at one un-rounded core rather than their kron.
        T = round(cell2core(tt_tensor, C(1:dim)), tol);
        C = core2cell(T);
    end

end


function H = contract_1d_integrals_with_factor_cores_direct(GD, GN, Usol, psol, nsol, Uden, pden, nden, Unum, pnum, nnum, spans, nq, dim, a, b)
% CONTRACT_1D_INTEGRALS_WITH_FACTOR_CORES_DIRECT
%
% Computes the final local stiffness TT core directly:
%
%   H(:,sub2ind([nsol,nsol],i,j),:) =
%       sum_{alpha,beta} Q(i,j,alpha,beta)
%       kron(GD(:,alpha,:), GN(:,beta,:)).
%
% This is evaluated in quadrature form as
%
%   sum_g w_g phi_i(x_g) phi_j(x_g)
%       kron(D_g, N_g),
%
% where
%
%   D_g = sum_alpha Phi_alpha^den(x_g) GD(:,alpha,:),
%   N_g = sum_beta  Psi_beta^num(x_g) GN(:,beta,:).

    [rd_l, nden_D, rd_r] = size(GD);
    [rn_l, nnum_N, rn_r] = size(GN);

    if nden_D ~= nden || nnum_N ~= nnum
        error(['contract_1d_integrals_with_factor_cores_direct: coefficient-core size mismatch. ', 'Denom core has %d, Num core has %d, expected %d and %d.'], nden_D, nnum_N, nden, nnum);
    end

    r_left  = rd_l * rn_l;
    r_right = rd_r * rn_r;
    H = zeros(r_left, nsol*nsol, r_right);

    for e = 1:size(spans,1)
        xa = spans(e,1);
        xb = spans(e,2);

        [xq, wq] = gauss_legendre(nq, xa, xb);
        xq = xq(:).';
        wq = wq(:);

        [Iu, Bu]  = eval_basis_on_span      (Usol, psol, xa, xb, xq);
        [~, dBu] = eval_basis_deriv_local(Usol, psol, xa, xb, xq);
        [Id, Bd]  = eval_basis_on_span      (Uden, pden, xa, xb, xq);
        [In, Bn]  = eval_basis_on_span      (Unum, pnum, xa, xb, xq);

        if a == dim
            Phi_i = dBu;
        else
            Phi_i = Bu;
        end

        if b == dim
            Phi_j = dBu;
        else
            Phi_j = Bu;
        end

        GD_loc = GD(:,Id,:);
        GN_loc = GN(:,In,:);

        for g = 1:nq
            Dg = weighted_sum_core_slices(GD_loc, Bd(g,:));
            Ng = weighted_sum_core_slices(GN_loc, Bn(g,:));
            Kg = kron(Dg, Ng);

            Wloc = wq(g) * (Phi_i(g,:).' * Phi_j(g,:));

            for ii = 1:numel(Iu)
                iglob = Iu(ii);
                for jj = 1:numel(Iu)
                    val = Wloc(ii,jj);
                    if val ~= 0
                        jglob = Iu(jj);
                        row = sub2ind([nsol,nsol], iglob, jglob);
                        H(:,row,:) = H(:,row,:) + reshape(val * Kg, [r_left, 1, r_right]);
                    end
                end
            end
        end
    end

    % H is returned at full (r_left, nsol^2, r_right). It is NOT truncated
    % here: the rank collapse happens in apply_qintegrals_...'s incremental
    % tkron/round loop, which owns the TT bond bookkeeping. Truncating the
    % right bond in place without folding the residual into the next core's
    % left bond corrupts the shared bond and breaks round().
end


function A = weighted_sum_core_slices(G, w)
% WEIGHTED_SUM_CORE_SLICES
% G has size r_left x nloc x r_right and w has length nloc.
% Returns sum_k w(k) * squeeze(G(:,k,:)).

    [r_left, nloc, r_right] = size(G);
    A = zeros(r_left, r_right);

    for k = 1:nloc
        if w(k) ~= 0
            A = A + w(k) * reshape(G(:,k,:), [r_left, r_right]);
        end
    end
end