function stats = collect_stiffness_minimal_stats(time, Bmat, det_ops, A_det_raw_tensor_tt, A_det_tensor_tt, A_det_tt, b, W, Denom_tt, geo_coeffs, ops_kind, kind_names, K_tt)
%COLLECT_STIFFNESS_MINIMAL_STATS  Paper-level statistics after timing.
%
% The goal is to report only meaningful groups of objects, without counting
% different representations of the same mathematical object twice.
%
% Main groups:
%   1) coefficient-transfer operators,
%   2) projection tensors,
%   3) post-projection tensors.
%
% Not counted in the main totals:
%   A_det_tensor_tt : raw determinant tensor retained only as a diagnostic,
%   univariate rhs  : omitted after the TT right-hand side b has been formed.
%
% All RecursiveSize/nnz calls are made after stats.time has been fixed.

    stats = struct();
    stats.time = time;
    stats.collection = 'after_toc';
    stats.note = ['Main totals avoid double-counting temporary and equivalent ', 'representations. The reported projection operator is A_det_tt.'];

    % ============================================================
    % 1. Coefficient-transfer operators.
    % ============================================================
    % These are the constructed one-dimensional data replacing the large
    % explicit transfer matrices F^(d) and numerator transfer matrices.
    [det_coeff_mem, det_coeff_nnz, det_index_mem, det_index_nnz, det_by_dim] = det_transfer_minimal_stats(det_ops);

    num_stats = numerator_transfer_minimal_stats(ops_kind, kind_names);

    stats.transfer_det_memory = det_coeff_mem + det_index_mem;
    stats.transfer_det_nnz    = det_coeff_nnz + det_index_nnz;
    stats.transfer_num_memory = num_stats.coeff_memory + num_stats.index_memory;
    stats.transfer_num_nnz    = num_stats.coeff_nnz + num_stats.index_nnz;

    stats.transfer_total_memory = stats.transfer_det_memory + stats.transfer_num_memory;
    stats.transfer_total_nnz    = stats.transfer_det_nnz    + stats.transfer_num_nnz;

    % Optional breakdowns, useful for debugging but not necessary for a paper table.
    stats.transfer_breakdown.det_coeff_memory = det_coeff_mem;
    stats.transfer_breakdown.det_coeff_nnz    = det_coeff_nnz;
    stats.transfer_breakdown.det_index_memory = det_index_mem;
    stats.transfer_breakdown.det_index_nnz    = det_index_nnz;
    stats.transfer_breakdown.det_by_dim       = det_by_dim;

    stats.transfer_breakdown.num_coeff_memory = num_stats.coeff_memory;
    stats.transfer_breakdown.num_coeff_nnz    = num_stats.coeff_nnz;
    stats.transfer_breakdown.num_index_memory = num_stats.index_memory;
    stats.transfer_breakdown.num_index_nnz    = num_stats.index_nnz;
    stats.transfer_breakdown.num_by_kind      = num_stats.by_kind;

    % ============================================================
    % 2. Projection tensors.
    % ============================================================
    % Count the final TT-matrix projection operator, the TT right-hand side,
    % and the univariate congruence transforms used by the AMEn solve.
    stats.projection_matrix_memory = stats_bytes(A_det_tt);
    stats.projection_matrix_nnz    = stats_nnz(A_det_tt);
    stats.projection_rhs_memory    = stats_bytes(b);
    stats.projection_rhs_nnz       = stats_nnz(b);
    stats.projection_congruence_memory = stats_bytes(W);
    stats.projection_congruence_nnz = stats_nnz(W);

    stats.projection_total_memory = stats.projection_congruence_memory + stats.projection_matrix_memory + stats.projection_rhs_memory;
    stats.projection_total_nnz = stats.projection_congruence_nnz + stats.projection_matrix_nnz + stats.projection_rhs_nnz;

    % Optional diagnostic: univariate integration operators used to build the
    % projection matrix.  These are not added to projection_total_* because
    % they are operators, not assembled projection tensors.
    stats.projection_operator_memory = stats_bytes(Bmat);
    stats.projection_operator_nnz    = stats_nnz(Bmat);

    % Optional diagnostic: raw determinant tensor before projection.  This is
    % useful if you want to report the intermediate determinant construction.
    stats.det_raw_memory = stats_bytes(A_det_raw_tensor_tt);
    stats.det_raw_nnz    = stats_nnz(A_det_raw_tensor_tt);

    % Optional diagnostic retained for compatibility with existing result files.
    stats.excluded_AprojTensor_memory = stats_bytes(A_det_tensor_tt);
    stats.excluded_AprojTensor_nnz    = stats_nnz(A_det_tensor_tt);

    % ============================================================
    % 3. Post-projection tensors.
    % ============================================================
    % These are the TT objects used after solving the projection problem.
    % The geometry coefficients are included here because they are needed for
    % the numerator construction.  K_tt is the final assembled stiffness tensor.
    stats.geometry_coeff_memory = stats_bytes(geo_coeffs);
    stats.geometry_coeff_nnz    = stats_nnz(geo_coeffs);
    stats.denom_memory          = stats_bytes(Denom_tt);
    stats.denom_nnz             = stats_nnz(Denom_tt);
    stats.K_memory              = stats_bytes(K_tt);
    stats.K_nnz                 = stats_nnz(K_tt);

    stats.post_projection_total_memory = stats.geometry_coeff_memory + stats.denom_memory;
    stats.post_projection_total_nnz = stats.geometry_coeff_nnz + stats.denom_nnz;

    % ============================================================
    % Paper-level total without deliberate double-counting.
    % ============================================================
    stats.paper_total_memory = stats.transfer_total_memory + stats.projection_total_memory + stats.post_projection_total_memory;

    stats.paper_total_nnz = stats.transfer_total_nnz + stats.projection_total_nnz + stats.post_projection_total_nnz;

    % Compatibility aliases for older plotting scripts.  These are simple
    % aliases, not additional counted objects.
    stats.C_memory        = stats.geometry_coeff_memory;
    stats.C_nnz           = stats.geometry_coeff_nnz;
    stats.Denom_memory    = stats.denom_memory;
    stats.Denom_nnz       = stats.denom_nnz;
    stats.AdetMatrix_memory = stats.projection_matrix_memory;
    stats.AdetMatrix_nnz    = stats.projection_matrix_nnz;
    stats.ProjRhsTT_memory  = stats.projection_rhs_memory;
    stats.ProjRhsTT_nnz     = stats.projection_rhs_nnz;
end

function [coeff_memory, coeff_nnz, index_memory, index_nnz, by_dim] = det_transfer_minimal_stats(det_ops)
%DET_TRANSFER_MINIMAL_STATS
% Counts the determinant coefficient transfer used instead of F^(d).
%
% Coefficient maps per dimension:
%   D, Pcoef, Rcoef
% Index maps per dimension:
%   IJ, AM

    coeff_memory = 0;
    coeff_nnz = 0;
    index_memory = 0;
    index_nnz = 0;
    by_dim = struct();

    for dim = 1:numel(det_ops)
        ops = det_ops{dim};

        c_mem = stats_bytes(ops.D) + stats_bytes(ops.Pcoef) + stats_bytes(ops.Rcoef);
        c_nnz = stats_nnz(ops.D)  + stats_nnz(ops.Pcoef)  + stats_nnz(ops.Rcoef);
        i_mem = stats_bytes(ops.IJ) + stats_bytes(ops.AM);
        i_nnz = stats_nnz(ops.IJ)  + stats_nnz(ops.AM);

        coeff_memory = coeff_memory + c_mem;
        coeff_nnz    = coeff_nnz    + c_nnz;
        index_memory = index_memory + i_mem;
        index_nnz    = index_nnz    + i_nnz;

        f = sprintf('dim%d', dim);
        by_dim.(f).coeff_memory = c_mem;
        by_dim.(f).coeff_nnz    = c_nnz;
        by_dim.(f).index_memory = i_mem;
        by_dim.(f).index_nnz    = i_nnz;
        by_dim.(f).total_memory = c_mem + i_mem;
        by_dim.(f).total_nnz    = c_nnz + i_nnz;
    end
end

function out = numerator_transfer_minimal_stats(ops_kind, kind_names)
%NUMERATOR_TRANSFER_MINIMAL_STATS
% Counts all numerator coefficient transfers that are constructed in ops_kind.
%
% Counted once by kind and dimension, not through ops_num, because ops_num
% reuses these objects pattern-wise.

    out = struct();
    out.coeff_memory = 0;
    out.coeff_nnz = 0;
    out.index_memory = 0;
    out.index_nnz = 0;
    out.by_kind = struct();

    for kk = 1:numel(kind_names)
        kind = kind_names{kk};
        out.by_kind.(kind) = struct();
        out.by_kind.(kind).coeff_memory = 0;
        out.by_kind.(kind).coeff_nnz = 0;
        out.by_kind.(kind).index_memory = 0;
        out.by_kind.(kind).index_nnz = 0;
        out.by_kind.(kind).total_memory = 0;
        out.by_kind.(kind).total_nnz = 0;
        out.by_kind.(kind).by_dim = struct();

        for dim = 1:numel(ops_kind.(kind))
            ops = ops_kind.(kind){dim};
            [c_mem, c_nnz, i_mem, i_nnz] = one_numerator_transfer_minimal_stats(ops);

            out.coeff_memory = out.coeff_memory + c_mem;
            out.coeff_nnz    = out.coeff_nnz    + c_nnz;
            out.index_memory = out.index_memory + i_mem;
            out.index_nnz    = out.index_nnz    + i_nnz;

            out.by_kind.(kind).coeff_memory = out.by_kind.(kind).coeff_memory + c_mem;
            out.by_kind.(kind).coeff_nnz    = out.by_kind.(kind).coeff_nnz    + c_nnz;
            out.by_kind.(kind).index_memory = out.by_kind.(kind).index_memory + i_mem;
            out.by_kind.(kind).index_nnz    = out.by_kind.(kind).index_nnz    + i_nnz;

            f = sprintf('dim%d', dim);
            out.by_kind.(kind).by_dim.(f).coeff_memory = c_mem;
            out.by_kind.(kind).by_dim.(f).coeff_nnz    = c_nnz;
            out.by_kind.(kind).by_dim.(f).index_memory = i_mem;
            out.by_kind.(kind).by_dim.(f).index_nnz    = i_nnz;
            out.by_kind.(kind).by_dim.(f).total_memory = c_mem + i_mem;
            out.by_kind.(kind).by_dim.(f).total_nnz    = c_nnz + i_nnz;
        end

        out.by_kind.(kind).total_memory = out.by_kind.(kind).coeff_memory + out.by_kind.(kind).index_memory;
        out.by_kind.(kind).total_nnz = out.by_kind.(kind).coeff_nnz + out.by_kind.(kind).index_nnz;
    end
end

function [coeff_memory, coeff_nnz, index_memory, index_nnz] = one_numerator_transfer_minimal_stats(ops)
%ONE_NUMERATOR_TRANSFER_MINIMAL_STATS
% Counts coefficient maps and index maps for one univariate numerator transfer.

    switch ops.kind
        case 'BBBB'
            coeff_fields = {'Pbb','P22'};
            index_fields = {'IJbb','IJ22'};
        case 'DDBB'
            coeff_fields = {'D','Pdd','Pbb','Rcoef'};
            index_fields = {'IJdd','IJbb','AM'};
        case 'DBBB'
            coeff_fields = {'D','Rdb','Pbb','Rcoef'};
            index_fields = {'AMdb','IJbb','AM'};
        otherwise
            error('one_numerator_transfer_minimal_stats: unknown kind "%s".', ops.kind);
    end

    coeff_memory = 0;
    coeff_nnz = 0;
    index_memory = 0;
    index_nnz = 0;

    for k = 1:numel(coeff_fields)
        f = coeff_fields{k};
        coeff_memory = coeff_memory + stats_bytes(ops.(f));
        coeff_nnz    = coeff_nnz    + stats_nnz(ops.(f));
    end

    for k = 1:numel(index_fields)
        f = index_fields{k};
        index_memory = index_memory + stats_bytes(ops.(f));
        index_nnz    = index_nnz    + stats_nnz(ops.(f));
    end
end

function b = stats_bytes(x)
% STATS_BYTES
% Uses RecursiveSize if it exists on the MATLAB path. Otherwise, returns a
% conservative recursive byte count for ordinary MATLAB arrays/cells/structs
% and the whos byte count for opaque objects.

    if isempty(x)
        b = 0;
        return;
    end

    if exist('RecursiveSize', 'file') == 2 || exist('RecursiveSize', 'builtin') == 5
        b = RecursiveSize(x);
        return;
    end

    b = local_recursive_size(x);
end

function b = local_recursive_size(x)
    if isempty(x)
        b = 0;
    elseif isnumeric(x) || islogical(x) || ischar(x) || isa(x, 'function_handle')
        info = whos('x');
        b = info.bytes;
    elseif iscell(x)
        b = 0;
        for k = 1:numel(x)
            b = b + local_recursive_size(x{k});
        end
    elseif isstruct(x)
        b = 0;
        f = fieldnames(x);
        for j = 1:numel(x)
            for k = 1:numel(f)
                b = b + local_recursive_size(x(j).(f{k}));
            end
        end
    else
        % TT-toolbox objects are opaque here unless RecursiveSize is available.
        info = whos('x');
        b = info.bytes;
    end
end

function n = stats_nnz(x)
% STATS_NNZ
% Recursive nonzero counter. For TT-toolbox objects it first tries nnz(x). If
% this is unavailable, it counts nonzeros in the TT cores.

    if isempty(x)
        n = 0;
    elseif iscell(x)
        n = 0;
        for k = 1:numel(x)
            n = n + stats_nnz(x{k});
        end
    elseif isstruct(x)
        n = 0;
        f = fieldnames(x);
        for j = 1:numel(x)
            for k = 1:numel(f)
                n = n + stats_nnz(x(j).(f{k}));
            end
        end
    elseif isnumeric(x) || islogical(x)
        n = nnz(x);
    else
        try
            n = nnz(x);
        catch
            try
                C = core2cell(x);
                n = stats_nnz(C);
            catch
                n = NaN;
            end
        end
    end
end
