function T = tt_mode_apply(T, W, transp)
% Apply W' (transp = true) or W (transp = false) to every mode of a
% TT tensor.  TT ranks are unchanged.

    C = core2cell(T);

    for d = 1:numel(C)
        G  = C{d};
        r1 = size(G,1); n = size(G,2); r2 = size(G,3);

        P = reshape(permute(G, [2 1 3]), n, r1*r2);
        if transp
            P = W{d}.' * P;
        else
            P = W{d} * P;
        end

        C{d} = permute(reshape(P, size(P,1), r1, r2), [2 1 3]);
    end

    T = cell2core(tt_tensor, C);
end
