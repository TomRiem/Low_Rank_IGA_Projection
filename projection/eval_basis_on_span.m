function [I, B] = eval_basis_on_span(U, p, xa, xb, xq)
% Evaluate the p+1 active B-splines on one positive knot span.

    xmid = 0.5*(xa+xb);
    span = find_span(U,p,xmid);
    I = (span-p):span;

    B = zeros(numel(xq),p+1);
    for q = 1:numel(xq)
        B(q,:) = basis_funs(span,xq(q),p,U);
    end
end