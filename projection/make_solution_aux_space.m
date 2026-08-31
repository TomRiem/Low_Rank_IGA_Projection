function aux_space = make_solution_aux_space(space)

    dim = numel(space.degree);

    aux_space = struct;
    aux_space.knots = cell(dim, 1);
    aux_space.degree = zeros(1, dim);
    aux_space.ndof_dir = zeros(1, dim);

    for d = 1:dim
        % Exact determinant space of the original geometry.
        [U_aux, p_aux] = enlargen_bspline_space(space.knots{d}, space.degree(d), 3);


        aux_space.knots{d} = U_aux;
        aux_space.degree(d) = p_aux;
        aux_space.ndof_dir(d) = numel(U_aux) - p_aux - 1;
    end
end