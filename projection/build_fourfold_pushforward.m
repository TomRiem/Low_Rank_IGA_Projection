function ops = build_fourfold_pushforward(U, p, kind)
% BUILD_FOURFOLD_PUSHFORWARD
% Build the one-dimensional coefficient transfer for one of the three
% four-factor numerator patterns used in the metric tensor:
%   BBBB, DDBB, or DBBB.
%
% The transfer is kept factorized into derivative and product operators.

    switch kind
        case 'BBBB'
            U2 = build_product_space_knots(U, p, U, p);
            U4 = build_product_space_knots(U2, 2*p, U2, 2*p);

            [IJbb, Pbb] = build_product_operator(U, p, U, p, U2, 2*p, true);
            [IJ22, P22] = build_product_operator(U2, 2*p, U2, 2*p, U4, 4*p, true);

            ops.kind = kind;
            ops.U = U;
            ops.p = p;
            ops.U2 = U2;
            ops.Uout = U4;
            ops.pout = 4*p;
            ops.nout = numel(U4) - ops.pout - 1;

            ops.Pbb = Pbb;
            ops.P22 = P22;
            ops.IJbb = IJbb;
            ops.IJ22 = IJ22;

            ops.bb_i = IJbb(:,1);
            ops.bb_j = IJbb(:,2);
            ops.bb_diag = (ops.bb_i == ops.bb_j);
            ops.p22_i = IJ22(:,1);
            ops.p22_j = IJ22(:,2);
            ops.p22_diag = (ops.p22_i == ops.p22_j);

            ops.nin = numel(U) - p - 1;
            ops.n2 = size(Pbb,1);

        case 'DDBB'
            [D, Ud] = build_derivative_operator(U, p);

            Udd = build_product_space_knots(Ud, p-1, Ud, p-1);
            U2  = build_product_space_knots(U, p, U, p);
            U42 = build_product_space_knots(Udd, 2*p-2, U2, 2*p);

            [IJdd, Pdd] = build_product_operator(Ud, p-1, Ud, p-1, Udd, 2*p-2, true);
            [IJbb, Pbb] = build_product_operator(U, p, U, p, U2, 2*p, true);
            [AM, Rcoef] = build_product_operator(Udd, 2*p-2, U2, 2*p, U42, 4*p-2, false);

            ops.kind = kind;
            ops.U = U;
            ops.p = p;
            ops.Ud = Ud;
            ops.Udd = Udd;
            ops.U2 = U2;
            ops.Uout = U42;
            ops.pout = 4*p - 2;
            ops.nout = numel(U42) - ops.pout - 1;

            ops.D = D;
            ops.Pdd = Pdd;
            ops.Pbb = Pbb;
            ops.Rcoef = Rcoef;
            ops.IJdd = IJdd;
            ops.IJbb = IJbb;
            ops.AM = AM;

            ops.dd_i = IJdd(:,1);
            ops.dd_j = IJdd(:,2);
            ops.dd_diag = (ops.dd_i == ops.dd_j);
            ops.bb_i = IJbb(:,1);
            ops.bb_j = IJbb(:,2);
            ops.bb_diag = (ops.bb_i == ops.bb_j);
            ops.am_dd = AM(:,1);
            ops.am_bb = AM(:,2);

            ops.nin = numel(U) - p - 1;
            ops.ndd = size(Pdd,1);
            ops.n2 = size(Pbb,1);

        case 'DBBB'
            [D, Ud] = build_derivative_operator(U, p);

            U2   = build_product_space_knots(U, p, U, p);
            Udb  = build_product_space_knots(Ud, p-1, U, p);
            Uout = build_product_space_knots(Udb, 2*p-1, U2, 2*p);

            [AMdb, Rdb] = build_product_operator(Ud, p-1, U, p, Udb, 2*p-1, false);
            [IJbb, Pbb] = build_product_operator(U, p, U, p, U2, 2*p, true);
            [AM, Rcoef] = build_product_operator(Udb, 2*p-1, U2, 2*p, Uout, 4*p-1, false);

            ops.kind = kind;
            ops.U = U;
            ops.p = p;
            ops.Ud = Ud;
            ops.U2 = U2;
            ops.Udb = Udb;
            ops.Uout = Uout;
            ops.pout = 4*p - 1;
            ops.nout = numel(Uout) - ops.pout - 1;

            ops.D = D;
            ops.Rdb = Rdb;
            ops.Pbb = Pbb;
            ops.Rcoef = Rcoef;
            ops.AMdb = AMdb;
            ops.IJbb = IJbb;
            ops.AM = AM;

            ops.am_d = AMdb(:,1);
            ops.am_b = AMdb(:,2);
            ops.bb_i = IJbb(:,1);
            ops.bb_j = IJbb(:,2);
            ops.bb_diag = (ops.bb_i == ops.bb_j);
            ops.am_db = AM(:,1);
            ops.am_bb = AM(:,2);

            ops.nin = numel(U) - p - 1;
            ops.ndb = size(Rdb,1);
            ops.n2 = size(Pbb,1);

        otherwise
            error('Unknown fourfold pushforward kind "%s".', kind);
    end
end