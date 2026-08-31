function Uout = build_product_space_knots(UA, pA, UB, pB)
% BUILD_PRODUCT_SPACE_KNOTS
% Knot vector of a spline space containing products from S_pA(UA)*S_pB(UB).

    z = unique(sort([UA(:); UB(:)])).';
    pout = pA + pB;
    Uout = [];

    for r = 1:numel(z)
        xi = z(r);

        if r == 1 || r == numel(z)
            mult_out = pout + 1;
        else
            mult_A = sum(abs(UA - xi) < 1e-14);
            mult_B = sum(abs(UB - xi) < 1e-14);

            % The product has the minimum continuity of the two factors.
            mult_out = max(pB + mult_A, pA + mult_B);
        end

        Uout = [Uout, repmat(xi, 1, mult_out)]; 
    end
end