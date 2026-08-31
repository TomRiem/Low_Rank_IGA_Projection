function N = basis_funs(span, x, p, U)
% Cox-de Boor evaluation of the nonzero B-splines on a given span.

    N = zeros(1,p+1);
    N(1) = 1;

    left = zeros(1,p);
    right = zeros(1,p);

    for j = 1:p
        left(j) = x - U(span+1-j);
        right(j) = U(span+j) - x;

        saved = 0;
        for r = 1:j
            denom = right(r) + left(j-r+1);
            if abs(denom) < 1e-14
                temp = 0;
            else
                temp = N(r)/denom;
            end

            N(r) = saved + right(r)*temp;
            saved = left(j-r+1)*temp;
        end
        N(j+1) = saved;
    end
end