function H = univariate_u_v_nurbs(H, space, Tweights)
    s = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664];
    w = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]';
    H.mass.M = cell(3,1);
    for dim = 1:3
        H.mass.M{dim} = cell(H.mass.R(dim),1);
        H.mass.M{dim}(:) = {sparse(space.ndof_dir(dim), ...
            space.ndof_dir(dim))};
        for l=1:length(space.knots{dim})-1
            if space.knots{dim}(l) < space.knots{dim}(l+1)
                a = space.knots{dim}(l);
                b = space.knots{dim}(l+1);
                xx = (b-a)/2*s + (a+b)/2;
                quadValues = evalNURBS(space.knots{dim}, space.degree(dim), Tweights{dim}', xx);
                quadValues2 = evalBSpline(H.weightFun.knots{dim}, H.weightFun.degree(dim), xx);
                for i = l-space.degree(dim):l
                    for j = l-space.degree(dim):l
                        for r = 1:H.mass.R(dim)
                            H.mass.M{dim}{r}(i,j) = ...
                                H.mass.M{dim}{r}(i,j) + ...
                                ((b-a)/2)*sum(w.*quadValues(i,:)'.*quadValues(j,:)'.*quadValues2'*H.mass.SVDU{dim}(:,r));
                        end
                    end
                end
            end
        end
    end
end