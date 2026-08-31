function [space, msh] = initialize_laplace(geometry, method_data)


    if ~isfield(method_data,'degree') || isempty(method_data.degree)
        method_data.degree = geometry.nurbs.order-1;
        method_data.regularity = geometry.nurbs.order-2;
    end
    if ~isfield(method_data,'regularity') || isempty(method_data.regularity)
        method_data.regularity = method_data.degree-1;
    end
    if ~isfield(method_data,'nsub') || isempty(method_data.nsub)
        method_data.nsub = [0, 0, 0];
    end


    if isfield(method_data,'basis_functions') && strcmp(method_data.basis_functions, 'B-Splines')
        [knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub, method_data.degree, method_data.regularity);
        rule     = msh_gauss_nodes (method_data.nquad);
        [qn, qw] = msh_set_quad_nodes (zeta, rule);
        msh      = msh_cartesian (zeta, qn, qw, geometry);
        space = sp_bspline (knots, method_data.degree, msh);
    else
        nurbs = nrbdegelev(geometry.nurbs, method_data.degree - geometry.nurbs.order-1);
        [~, zeta, new_knots] = kntrefine (nurbs.knots, method_data.nsub, nurbs.order-1, nurbs.order-2);
        nurbs = nrbkntins(nurbs, new_knots);
        
        rule     = msh_gauss_nodes (method_data.nquad);
        [qn, qw] = msh_set_quad_nodes (zeta, rule);
        msh      = msh_cartesian (zeta, qn, qw, geometry);
        space = sp_nurbs (nurbs, msh);
    end
end