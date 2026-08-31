function A = tt_tensor_to_tt_matrix(T, nrow, ncol)
% TT_TENSOR_TO_TT_MATRIX
% Interpret tensor modes nrow(d)*ncol(d) as paired row/column modes.

    nrow = nrow(:);
    ncol = ncol(:);

    if numel(nrow) ~= numel(ncol)
        error('nrow and ncol must have the same number of modes.');
    end

    if any(T.n(:) ~= nrow .* ncol)
        error('Mode sizes of T do not match nrow .* ncol.');
    end

    A = tt_matrix(T, nrow, ncol);
end