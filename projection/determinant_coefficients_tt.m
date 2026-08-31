function [Det_raw_tt, C_1_tt, C_2_tt, C_3_tt] = determinant_coefficients_tt(tol, geometry)
% DETERMINANT_COEFFICIENTS_TT
% Construct the determinant coefficient tensor from the three coordinate coefficient
% tensors of the B-spline geometry map.

    C_1 = reshape(geometry.nurbs.coefs(1,:,:,:), geometry.nurbs.number);
    C_2 = reshape(geometry.nurbs.coefs(2,:,:,:), geometry.nurbs.number);
    C_3 = reshape(geometry.nurbs.coefs(3,:,:,:), geometry.nurbs.number);

    C_1_tt = tt_tensor(C_1, tol);
    C_2_tt = tt_tensor(C_2, tol);
    C_3_tt = tt_tensor(C_3, tol);

    % Leibniz expansion of det(grad G).
    C1_C2_C3 = round(tkron(C_1_tt, tkron(C_2_tt, C_3_tt)), tol);
    C3_C1_C2 = round(tkron(C_3_tt, tkron(C_1_tt, C_2_tt)), tol);
    C2_C3_C1 = round(tkron(C_2_tt, tkron(C_3_tt, C_1_tt)), tol);

    C3_C2_C1 = round(tkron(C_3_tt, tkron(C_2_tt, C_1_tt)), tol);
    C1_C3_C2 = round(tkron(C_1_tt, tkron(C_3_tt, C_2_tt)), tol);
    C2_C1_C3 = round(tkron(C_2_tt, tkron(C_1_tt, C_3_tt)), tol);

    Det_raw_tt =   C1_C2_C3 + C3_C1_C2 + C2_C3_C1 - C3_C2_C1 - C1_C3_C2 - C2_C1_C3;
    Det_raw_tt = round(Det_raw_tt, tol);

    % Group the three factors belonging to each parametric direction.
    Det_raw_tt = permute(Det_raw_tt, [7, 4, 1, 2, 8, 5, 3, 6, 9], tol);
    Det_raw_tt = round(Det_raw_tt, tol);

    Det_raw_tt = reshape(Det_raw_tt, [prod(Det_raw_tt.n(1:3));  prod(Det_raw_tt.n(4:6));  prod(Det_raw_tt.n(7:9))]);

    Det_raw_tt = round(Det_raw_tt, tol);
end