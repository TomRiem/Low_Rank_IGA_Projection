function [Bmat, rhs] = assemble_projection_Bmat_rhs_merged(Uproj, pproj, Udet, pdet, tol)
% ASSEMBLE_PROJECTION_BMAT_RHS_MERGED
%
% Assemble
%   rhs_i        = \int phi_i
%   Bmat_(ij),k  = \int phi_i phi_j psi_k
%
% where
%   {phi_i} is the projection basis S_pproj(Uproj),
%   {psi_k} is the determinant basis S_pdet(Udet).
%
% Integration is performed on the merged knot partition of Uproj and Udet.

    if nargin < 5 || isempty(tol)
        tol = 1e-14;
    end

    nproj = numel(Uproj) - pproj - 1;
    ndet  = numel(Udet)  - pdet  - 1;

    rhs = zeros(nproj, 1);

    % Exact Gauss rule for degree (2*pproj + pdet)
    nq = ceil((2*pproj + pdet + 1)/2);

    spans = positive_knot_spans(merge_knots(Uproj, Udet));

    chunk_size = 200000;
    subs = zeros(chunk_size, 3);
    vals = zeros(chunk_size, 1);
    count = 0;

    for e = 1:size(spans,1)
        xa = spans(e,1);
        xb = spans(e,2);

        [xq, wq] = gauss_legendre(nq, xa, xb);
        xq = xq(:).';
        wq = wq(:);
        
        % Local basis evaluations on the current merged span.
        % No global arrays of size nproj x nq or ndet x nq are formed.
        [Iproj, Phi] = eval_basis_on_span(Uproj, pproj, xa, xb, xq);   % Phi: nq x nproj_loc
        [Idet, Psi] = eval_basis_on_span(Udet, pdet, xa, xb, xq);   % Psi: nq x ndet_loc
        
        nloc_proj = numel(Iproj);
        nloc_det  = numel(Idet);

        % rhs contribution: int phi_i
        rhs(Iproj) = rhs(Iproj) + Phi.' * wq;

        % Bmat contribution
        for kk = 1:nloc_det
            wk = wq .* Psi(:,kk);              % nq x 1
            Mloc = Phi.' * (Phi .* wk);        % nloc_proj x nloc_proj

            [ii, jj] = find(abs(Mloc) > tol);
            if ~isempty(ii)
                lin = sub2ind(size(Mloc), ii, jj);
                vv  = Mloc(lin);
            
                m = numel(vv);
            
                subs_new = zeros(m, 3);
                subs_new(:,1) = Iproj(ii);
                subs_new(:,2) = Iproj(jj);
                subs_new(:,3) = Idet(kk);
            
                [subs, vals, count] = append_triplets(subs, vals, count, subs_new, vv, chunk_size);
            end
        end
    end

    subs = subs(1:count,:);
    vals = vals(1:count);

    rows = sub2ind([nproj, nproj], subs(:,1), subs(:,2));
    cols = subs(:,3);

    % Merge duplicates once
    Bmat = sparse(rows, cols, vals, nproj^2, ndet);
end


function [subs_buf, vals_buf, count] = append_triplets(subs_buf, vals_buf, count, subs_new, vals_new, chunk_size)

    nnew = size(subs_new,1);
    need = count + nnew;

    if need > size(subs_buf,1)
        grow = max(chunk_size, need - size(subs_buf,1));
        subs_buf = [subs_buf; zeros(grow, 3)]; 
        vals_buf = [vals_buf; zeros(grow, 1)]; 
    end

    idx = count + (1:nnew);
    subs_buf(idx,:) = subs_new;
    vals_buf(idx)   = vals_new;
    count = count + nnew;
end