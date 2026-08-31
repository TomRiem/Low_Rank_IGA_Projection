function minors = make_metric_minors()
% Cofactors:
% A = ei - fh
% B = fg - di
% C = dh - eg
% D = ch - bi
% E = ai - cg
% F = bg - ah
% G = bf - ce
% H = cd - af
% I = ae - bd

    minors.A = [make_term(+1, {'e','i'}), make_term(-1, {'f','h'})];
    minors.B = [make_term(+1, {'f','g'}), make_term(-1, {'d','i'})];
    minors.C = [make_term(+1, {'d','h'}), make_term(-1, {'e','g'})];

    minors.D = [make_term(+1, {'c','h'}), make_term(-1, {'b','i'})];
    minors.E = [make_term(+1, {'a','i'}), make_term(-1, {'c','g'})];
    minors.F = [make_term(+1, {'b','g'}), make_term(-1, {'a','h'})];

    minors.G = [make_term(+1, {'b','f'}), make_term(-1, {'c','e'})];
    minors.H = [make_term(+1, {'c','d'}), make_term(-1, {'a','f'})];
    minors.I = [make_term(+1, {'a','e'}), make_term(-1, {'b','d'})];
end

function t = make_term(coeff, vars)
    t.coeff = coeff;
    t.vars  = vars;
end