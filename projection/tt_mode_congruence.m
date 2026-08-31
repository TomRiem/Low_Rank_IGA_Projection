function A = tt_mode_congruence(A, W)
% Congruence transform  A -> (W1'xW2'xW3') A (W1xW2xW3), core by core.
% TT ranks are unchanged.

    C = core2cell(A);

    for d = 1:numel(C)
        G  = C{d};
        r1 = size(G,1); n = size(G,2); m = size(G,3); r2 = size(G,4);
        k  = size(W{d},2);

        % row index
        P = reshape(permute(G, [2 1 3 4]), n, r1*m*r2);
        P = W{d}.' * P;
        G = permute(reshape(P, k, r1, m, r2), [2 1 3 4]);   % r1 x k x m x r2

        % column index
        P = reshape(permute(G, [3 1 2 4]), m, r1*k*r2);
        P = W{d}.' * P;
        C{d} = permute(reshape(P, k, r1, k, r2), [2 3 1 4]); % r1 x k x k x r2
    end

    A = cell2core(tt_matrix, C);
end