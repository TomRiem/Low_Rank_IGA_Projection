function proj = select_projection_space(U_geo, p_geo, proj_ref, space_proj, dim)
% SELECT_PROJECTION_SPACE
% Select the univariate space used for the L2 projection of 1/det(grad G).

    if ~isempty(proj_ref) && ~isempty(space_proj) && proj_ref == 1
        aux = build_sixfold_space_only(space_proj.knots{dim}, space_proj.degree(dim));
        proj.U = aux.U6;
        proj.p = aux.p_new;
        proj.n = numel(proj.U) - proj.p - 1;

    elseif ~isempty(proj_ref) && ~isempty(space_proj) && proj_ref == 2
        proj.U = space_proj.knots{dim};
        proj.p = space_proj.degree(dim);
        proj.n = space_proj.ndof_dir(dim);

    else
        aux = build_sixfold_space_only(U_geo, p_geo);
        proj.U = aux.U6;
        proj.p = aux.p_new;
        proj.n = numel(proj.U) - proj.p - 1;
    end
end
