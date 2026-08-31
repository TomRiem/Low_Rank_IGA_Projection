function M = assemble_univariate_mass(U, p)
% ASSEMBLE_UNIVARIATE_MASS
% Assemble the one-dimensional B-spline mass matrix in S_p(U).

    n = numel(U) - p - 1;
    M = sparse(n,n);

    spans = positive_knot_spans(U);
    nq = p + 1;

    for e = 1:size(spans,1)
        xa = spans(e,1);
        xb = spans(e,2);

        [xq, wq] = gauss_legendre(nq, xa, xb);
        [I, B] = eval_basis_on_span(U, p, xa, xb, xq);

        M(I,I) = M(I,I) + B.' * (B .* wq);
    end
end