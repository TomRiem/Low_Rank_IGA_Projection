function [TT_K, stats] = assemble_stiffness_lowrank_interpolation(H, space, low_rank_data, method_data)
    time = tic;

    x_cor = logical(ones(space.ndof_dir(1), 1));
    y_cor = logical(ones(space.ndof_dir(2), 1));
    z_cor = logical(ones(space.ndof_dir(3), 1));

    if method_data.dirichlet == 1
        x_cor(1) = 0;
        x_cor(end) = 0;
        y_cor(1) = 0;
        y_cor(end) = 0;
        z_cor(1) = 0;
        z_cor(end) = 0;
        TT_K = tt_zeros([space.ndof_dir'-[2; 2; 2], space.ndof_dir'-[2; 2; 2]]);
    else
        TT_K = tt_zeros([space.ndof_dir', space.ndof_dir']);
    end

    if strcmp(space.space_type, 'spline')
        H = univariate_gradu_gradv_bsplines(H, space, method_data.nquad);
        for i=1:9 
            for j = 1:H.stiffness.R(H.stiffness.order(i),1)
                for k = 1:H.stiffness.R(H.stiffness.order(i),3)
                    TT_K = round(TT_K + tt_matrix({full(H.stiffness.K{1}{i}{j}(x_cor, x_cor)); ...
                        full(H.stiffness.K{2}{i}{k+(j-1)*H.stiffness.R(H.stiffness.order(i),3)}(y_cor, y_cor)); ...
                        full(H.stiffness.K{3}{i}{k}(z_cor, z_cor))}), low_rank_data.rankTol);
                end
            end
        end
    
    else
        weights = reshape(space.weights, space.ndof_dir);
        Tweights = round(tt_tensor(weights), 1e-15, 1);

        H = univariate_gradu_gradv_nurbs(H, space, Tweights);
        for i=1:9 
            for j = 1:H.stiffness.R(H.stiffness.order(i),1)
                for k = 1:H.stiffness.R(H.stiffness.order(i),3)
                    TT_K = round(TT_K + tt_matrix({full(H.stiffness.K{1}{i}{j}(x_cor, x_cor)); ...
                        full(H.stiffness.K{2}{i}{k+(j-1)*H.stiffness.R(H.stiffness.order(i),3)}(y_cor, y_cor)); ...
                        full(H.stiffness.K{3}{i}{k}(z_cor, z_cor))}), low_rank_data.rankTol);
                end
            end
        end
    
    end

    time = toc(time);

    stats = struct; 

    stats.time = time; 
    stats.cores_memory = RecursiveSize(H.stiffness.K);

    stats.cores_nnz = 0;

    for i=1:9 
        for j = 1:H.stiffness.R(H.stiffness.order(i),1)
            for k = 1:H.stiffness.R(H.stiffness.order(i),3)
                stats.cores_nnz = stats.cores_nnz + nnz(H.stiffness.K{1}{i}{j}(x_cor, x_cor)) + nnz(H.stiffness.K{2}{i}{k+(j-1)*H.stiffness.R(H.stiffness.order(i),3)}(y_cor, y_cor)) + nnz(H.stiffness.K{3}{i}{k}(z_cor, z_cor));
            end
        end
    end
end