function [x, w] = gauss_legendre(n, a, b)
% Cached n-point Gauss-Legendre rule on [a,b].

    persistent xhat_cache what_cache

    if isempty(xhat_cache)
        xhat_cache = cell(1,200);
        what_cache = cell(1,200);
    end

    if n > numel(xhat_cache)
        xhat_cache{n} = [];
        what_cache{n} = [];
    end

    if isempty(xhat_cache{n})
        if n == 1
            xhat = 0;
            what = 2;
        else
            k = 1:n-1;
            beta = k ./ sqrt(4*k.^2 - 1);
            T = diag(beta,1) + diag(beta,-1);

            [V,D] = eig(T);
            xhat = diag(D);
            [xhat,idx] = sort(xhat);
            V = V(:,idx);
            what = 2 * (V(1,:).^2).';
        end

        xhat_cache{n} = xhat;
        what_cache{n} = what;
    else
        xhat = xhat_cache{n};
        what = what_cache{n};
    end

    x = 0.5*((b-a)*xhat + (a+b));
    w = 0.5*(b-a)*what;
end