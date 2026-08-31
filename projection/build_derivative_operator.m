function [D, Ud] = build_derivative_operator(U, p)
% BUILD_DERIVATIVE_OPERATOR
% Coefficient map from S_p(U) to its derivative space S_{p-1}(U(2:end-1)).

    if p < 1
        error('Degree p must be at least 1.');
    end

    Ud = U(2:end-1);

    n  = length(U) - p - 1;
    nd = n - 1;

    I = zeros(2*nd,1);
    J = zeros(2*nd,1);
    V = zeros(2*nd,1);

    t = 0;
    for mu = 1:nd
        denom = U(mu+p+1) - U(mu+1);
        if abs(denom) < 1e-14
            alpha = 0.0;
        else
            alpha = p / denom;
        end

        t = t + 1;
        I(t) = mu;
        J(t) = mu;
        V(t) = -alpha;

        t = t + 1;
        I(t) = mu;
        J(t) = mu+1;
        V(t) = alpha;
    end

    D = sparse(I, J, V, nd, n);
end