function [pairs, Pcoef] = build_product_operator(UA, pA, UB, pB, Uout, pout, use_symmetry)
% BUILD_PRODUCT_OPERATOR
% Compute coefficients of products of two B-spline basis functions in a
% prescribed product space. If both factors belong to the same space,
% use_symmetry=true stores only overlapping pairs i<=j.

    if nargin < 7
        use_symmetry = false;
    end

    if use_symmetry
        if pA ~= pB || numel(UA) ~= numel(UB) || any(abs(UA(:)-UB(:)) > 1e-14)
            error('Symmetric product storage requires identical spline spaces.');
        end
        pairs = build_overlap_pairs_same_space(UA, pA);
    else
        pairs = build_overlap_pairs_two_spaces(UA, pA, UB, pB);
    end

    nout   = numel(Uout) - pout - 1;
    npairs = size(pairs,1);

    Mout = assemble_univariate_mass(Uout, pout);
    RHS  = zeros(nout, npairs);

    for q = 1:npairs
        i = pairs(q,1);
        j = pairs(q,2);
        RHS(:,q) = assemble_product_rhs(UA, pA, i, UB, pB, j, Uout, pout);
    end

    Pcoef = Mout \ RHS;
    Pcoef = sparsify_matrix(Pcoef, 1e-14);
end

function pairs = build_overlap_pairs_same_space(U, p)
% Overlapping pairs i<=j from one spline space.

    n = numel(U) - p - 1;
    pairs = zeros(0,2);
    tol = 1e-14 * max(1, U(end)-U(1));

    for i = 1:n
        supp_i = [U(i), U(i+p+1)];
        for j = i:n
            supp_j = [U(j), U(j+p+1)];
            if min(supp_i(2),supp_j(2)) > max(supp_i(1),supp_j(1)) + tol
                pairs(end+1,:) = [i,j]; 
            end
        end
    end
end

function pairs = build_overlap_pairs_two_spaces(UA, pA, UB, pB)
% Overlapping basis-function pairs from two spline spaces.

    nA = numel(UA) - pA - 1;
    nB = numel(UB) - pB - 1;
    pairs = zeros(0,2);

    tol = 1e-14 * max(1, max(UA(end)-UA(1), UB(end)-UB(1)));

    for i = 1:nA
        supp_i = [UA(i), UA(i+pA+1)];
        for j = 1:nB
            supp_j = [UB(j), UB(j+pB+1)];
            if min(supp_i(2),supp_j(2)) > max(supp_i(1),supp_j(1)) + tol
                pairs(end+1,:) = [i,j]; 
            end
        end
    end
end

function rhs = assemble_product_rhs(UA, pA, i, UB, pB, j, Uout, pout)
% Assemble
%   rhs_a = int gamma_a beta_i^A beta_j^B,
% where gamma_a is a basis function of the target product space.

    nout = numel(Uout) - pout - 1;
    rhs = zeros(nout,1);

    spans = positive_knot_spans(merge_knots(UA, UB, Uout));
    nq = ceil((pA + pB + pout + 1)/2);

    same_space = pA == pB && numel(UA) == numel(UB) && all(abs(UA(:)-UB(:)) < 1e-14);

    for e = 1:size(spans,1)
        xa = spans(e,1);
        xb = spans(e,2);

        [xq, wq] = gauss_legendre(nq, xa, xb);

        [IA, BA] = eval_basis_on_span(UA, pA, xa, xb, xq);

        if same_space
            IB = IA;
            BB = BA;
        else
            [IB, BB] = eval_basis_on_span(UB, pB, xa, xb, xq);
        end

        [Io, Bo] = eval_basis_on_span(Uout, pout, xa, xb, xq);

        li = find(IA == i, 1);
        lj = find(IB == j, 1);

        if isempty(li) || isempty(lj)
            continue;
        end

        product_values = BA(:,li) .* BB(:,lj);
        rhs(Io) = rhs(Io) + Bo.' * (product_values .* wq);
    end
end

function A = sparsify_matrix(A, tol)
% Remove numerical noise from coefficient-transfer matrices.

    A(abs(A) < tol) = 0;
    A = sparse(A);
end