function [TT_M, stats] = assemble_mass_lowrank_interpolation(H, space, low_rank_data, method_data, rounding)

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
        TT_M = tt_zeros([space.ndof_dir'-[2; 2; 2], space.ndof_dir'-[2; 2; 2]]);
    else
        TT_M = tt_zeros([space.ndof_dir', space.ndof_dir']);
    end
    

    if strcmp(space.space_type, 'spline')
        H = univariate_u_v_bsplines(H, space, method_data.nquad);
        
        % for i=1:H.mass.R(1)
        %     for j = 1:H.mass.R(3)
        %         TT_M = round(TT_M + tt_matrix({full(H.mass.M{1}{i}(x_cor, x_cor)); ...
        %                 full(H.mass.M{2}{i+(j-1)*H.mass.R(1)}(y_cor, y_cor)); ...
        %                 full(H.mass.M{3}{j}(z_cor, z_cor))}), low_rank_data.rankTol);
        %     end
        % end

        if nargin > 4 && rounding == true
            for i = 1:H.mass.R(1)
                for j = 1:H.mass.R(3)
                    TT_M = round(TT_M + tt_matrix({full(H.mass.M{1}{i}(x_cor, x_cor)); ...
                            full(H.mass.M{2}{j+(i-1)*H.mass.R(3)}(y_cor, y_cor)); ...
                            full(H.mass.M{3}{j}(z_cor, z_cor))}), low_rank_data.rankTol);
                end
            end        
        else
            for i = 1:H.mass.R(1)
                for j = 1:H.mass.R(3)
                    TT_M = TT_M + tt_matrix({full(H.mass.M{1}{i}(x_cor, x_cor)); ...
                            full(H.mass.M{2}{j+(i-1)*H.mass.R(3)}(y_cor, y_cor)); ...
                            full(H.mass.M{3}{j}(z_cor, z_cor))});
                end
            end      
        end
    else
        weights = reshape(space.weights, space.ndof_dir);
        Tweights = round(tt_tensor(weights), 1e-15, 1);

        H = univariate_u_v_nurbs(H, space, Tweights);
        
        % for i=1:H.mass.R(1)
        %     for j = 1:H.mass.R(3)
        %         TT_M = round(TT_M + tt_matrix({full(H.mass.M{1}{i}(x_cor, x_cor)); ...
        %                 full(H.mass.M{2}{i+(j-1)*H.mass.R(1)}(y_cor, y_cor)); ...
        %                 full(H.mass.M{3}{j}(z_cor, z_cor))}), low_rank_data.rankTol);
        %     end
        % end

        if nargin > 4 && rounding == true
            for i = 1:H.mass.R(1)
                for j = 1:H.mass.R(3)
                    TT_M = round(TT_M + tt_matrix({full(H.mass.M{1}{i}(x_cor, x_cor)); ...
                            full(H.mass.M{2}{j+(i-1)*H.mass.R(3)}(y_cor, y_cor)); ...
                            full(H.mass.M{3}{j}(z_cor, z_cor))}), low_rank_data.rankTol);
                end
            end        
        else
            for i = 1:H.mass.R(1)
                for j = 1:H.mass.R(3)
                    TT_M = TT_M + tt_matrix({full(H.mass.M{1}{i}(x_cor, x_cor)); ...
                            full(H.mass.M{2}{j+(i-1)*H.mass.R(3)}(y_cor, y_cor)); ...
                            full(H.mass.M{3}{j}(z_cor, z_cor))});
                end
            end      
        end

    end

    time = toc(time);

    stats = struct; 

    stats.time = time; 
    stats.cores_memory = RecursiveSize(H.mass.M);

    stats.cores_nnz = 0;

    for i = 1:H.mass.R(1)
        for j = 1:H.mass.R(3)
            stats.cores_nnz = stats.cores_nnz + ...
                nnz(H.mass.M{1}{i}(x_cor, x_cor)) + ...
                nnz(H.mass.M{2}{j+(i-1)*H.mass.R(3)}(y_cor, y_cor)) + ...
                nnz(H.mass.M{3}{j}(z_cor, z_cor));
        end
    end    

end