function Bmat = assemble_mass_Bmat_merged(Usol, psol, Udet, pdet, tol, nq)
% ASSEMBLE_MASS_BMAT_MERGED
% Assemble the one-dimensional map from determinant coefficients to
% flattened mass-matrix entries,
%
%   Bmat((i,j),ell) = int beta_i beta_j tilde_beta_ell.

    if nargin < 5 || isempty(tol)
        tol = 1e-14;
    end

    nsol = numel(Usol) - psol - 1;
    ndet = numel(Udet) - pdet - 1;

    if nargin < 6 || isempty(nq)
        nq = ceil((2*psol + pdet + 1)/2);
    end

    spans = positive_knot_spans(merge_knots(Usol, Udet));

    chunk = 200000;
    rows = zeros(chunk,1);
    cols = zeros(chunk,1);
    vals = zeros(chunk,1);
    count = 0;

    for e = 1:size(spans,1)
        xa = spans(e,1);
        xb = spans(e,2);

        [xq, wq] = gauss_legendre(nq, xa, xb);
        wq = wq(:);

        [Iu, Bu] = eval_basis_on_span(Usol, psol, xa, xb, xq);
        [Ik, Bd] = eval_basis_on_span(Udet, pdet, xa, xb, xq);

        for k = 1:numel(Ik)
            wk = wq .* Bd(:,k);
            Mloc = Bu.' * (Bu .* wk);

            [ii, jj] = find(abs(Mloc) > tol);
            if isempty(ii)
                continue;
            end

            vv = Mloc(sub2ind(size(Mloc), ii, jj));
            m = numel(vv);

            need = count + m;
            if need > numel(rows)
                grow = max(chunk, m);
                rows = [rows; zeros(grow,1)]; 
                cols = [cols; zeros(grow,1)]; 
                vals = [vals; zeros(grow,1)]; 
            end

            idx = count + (1:m);
            rows(idx) = sub2ind([nsol,nsol], Iu(ii), Iu(jj));
            cols(idx) = Ik(k);
            vals(idx) = vv;
            count = count + m;
        end
    end

    Bmat = sparse(rows(1:count), cols(1:count), vals(1:count), nsol^2, ndet);
end