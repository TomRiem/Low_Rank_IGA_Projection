function stats = collect_mass_stats(time, Det_raw_tt, Bmat, det_ops)
% COLLECT_MASS_STATS
% Preserve the statistics fields used by the current numerical experiments.

    stats = struct();
    stats.time = time;

    stats.C_memory = RecursiveSize(Det_raw_tt);
    stats.C_nnz = nnz(Det_raw_tt);

    stats.B_memory = RecursiveSize(Bmat);
    stats.B_nnz = nnz(Bmat{1}) + nnz(Bmat{2}) + nnz(Bmat{3});

    stats.Proj_memory = 0;
    stats.Proj_nnz = 0;

    for dim = 1:3
        stats.Proj_memory = stats.Proj_memory + RecursiveSize(det_ops{dim}.D) + RecursiveSize(det_ops{dim}.Pcoef) + RecursiveSize(det_ops{dim}.Rcoef);

        stats.Proj_nnz = stats.Proj_nnz + nnz(det_ops{dim}.D) + nnz(det_ops{dim}.Pcoef) + nnz(det_ops{dim}.Rcoef);
    end
end