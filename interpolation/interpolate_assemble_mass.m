function [TT_M, stats] = interpolate_assemble_mass(geometry, low_rank_data, space, method_data, rounding)
    [H, stats_interpolate] = interpolate_system(geometry, low_rank_data);
    [TT_M, stats_assemble] = assemble_mass_lowrank_interpolation(H, space, low_rank_data, method_data, rounding);
    stats = struct; 
    stats.interpolate = stats_interpolate;
    stats.assemble = stats_assemble;
end

function [H, stats] = interpolate_system(geometry, low_rank_data)
% ADAPTIVITY_INTERPOLATION_SYSTEM
% Build low-rank (TT) interpolants for geometry-induced weights by the 
% THB-IGA low-rank assembly pipeline. Geometry may be B-splines
% or NURBS.
%
% [H, rhs, time] = ...
% ADAPTIVITY_INTERPOLATION_SYSTEM(geometry, low_rank_data, problem_data)
%
% Purpose
% -------
% Construct separated (tensor-train) representations needed for fast
% univariate quadrature on hierarchical levels:
% • H  – geometry-related weights (and auxiliary factors) in TT form,
% The routine also harmonizes degrees/regularities/subdivisions with the
% given geometry by degree elevation and knot refinement.
%
% Inputs
% ------
% geometry : GeoPDEs geometry struct (with .nurbs)
% • The routine degree-elevates and refines the knot vectors so the
%   requested interpolation spaces are available, then reloads via GEO_LOAD.
% • For NURBS, rational weights are extracted and tensorized.
%
% low_rank_data : struct (interpolation spaces & LR settings)
% • system_degree      – degree per direction for geometry weights
% • system_regularity  – continuity per direction for geometry weights
% • system_nsub        – additional dyadic subdivisions per direction
% • geometry_format    – 'B-Splines' to force the polynomial branch;
%                         otherwise NURBS is assumed
% • (any further fields are passed through to the LR helpers, e.g.,
%   rank tolerances, AMEn/TT options)
% Notes:
% • If degrees/regularities are omitted, defaults are inferred from
%   geometry.nurbs.order:  degree = order-1, regularity = degree-1.
% • *_nsub defaults to [0 0 0] if omitted.
%
%
% Outputs
% -------
% H   : struct with low-rank factors for geometry-induced weights used by
%       the univariate quadrature assembly (exact field layout matches the
%       downstream LR assembly routines, e.g. weightFun info and per-dir
%       SVD/TT factors for stiffness/mass entries).
% time: scalar, wall-clock seconds for the whole interpolation/rounding step.
%
% How it works
% ------------
% 1) Normalize interpolation choices
%    • Fill missing system_* fields from geometry.nurbs.order and
%      set *_nsub to [0 0 0] if absent.
% 2) Align the geometry to the requested spaces
%    • Degree elevation (NRBDEGELEV) and knot refinement (KNTREFINE + NRBKNTINS),
%      then reload the updated geometry (GEO_LOAD).
% 3) Prepare control points (+ weights for NURBS)
%    • B-splines: control points are taken as Cartesian; weights are 1.
%    • NURBS   : extract weights (4th homogeneous coord.), build a TT tensor
%      of weights (TT_TENSOR + rounding), and convert control points to
%      Cartesian by dividing by weights.
% 4) Interpolate geometry weights in low rank
%    • Call INTERPOLATE_WEIGHTS_BSPLINES or INTERPOLATE_WEIGHTS_NURBS to
%      obtain separable per-direction factors, then compress with LOWRANK_W.
% 5) Return H, rhs, and the elapsed time.
%
% Notes
% -----
% • Immediate TT rounding keeps intermediary ranks controlled prior to
%   assembly, which is crucial for performance.
% • The routine only prepares low-rank ingredients; the adaptive loop, error
%   estimator, marking and refinement remain those of GeoPDEs.

    time = tic;

    if ~isfield(low_rank_data,'system_nsub') || isempty(low_rank_data.system_nsub)
        low_rank_data.system_nsub = [0, 0, 0];
    end
    if ~isfield(low_rank_data,'system_degree') || isempty(low_rank_data.system_degree)
        low_rank_data.system_degree = geometry.nurbs.order-1;
        low_rank_data.system_regularity = geometry.nurbs.order-2;
    end
    if ~isfield(low_rank_data,'system_regularity') || isempty(low_rank_data.system_regularity)
        low_rank_data.system_regularity = low_rank_data.system_degree-1;
    end

    if ~isfield(low_rank_data,'rhs_nsub') || isempty(low_rank_data.rhs_nsub)
        low_rank_data.rhs_nsub = [0, 0, 0];
    end
    if ~isfield(low_rank_data,'rhs_degree') || isempty(low_rank_data.rhs_degree)
        low_rank_data.rhs_degree = geometry.nurbs.order-1;
        low_rank_data.rhs_regularity = geometry.nurbs.order-2;
    end
    if ~isfield(low_rank_data,'rhs_regularity') || isempty(low_rank_data.rhs_regularity)
        low_rank_data.rhs_regularity = low_rank_data.rhs_degree-1;
    end
    
    if isfield(low_rank_data,'geometry_format') && strcmp(low_rank_data.geometry_format, 'B-Splines')
        degelev = max (low_rank_data.system_degree - (geometry.nurbs.order-1), 0);
        nurbs = nrbdegelev (geometry.nurbs, degelev);

        [~, ~, new_knots] = kntrefine (nurbs.knots, low_rank_data.system_nsub, ...
            low_rank_data.system_degree, low_rank_data.system_regularity);

        nurbs = nrbkntins(nurbs, new_knots);

        geometry = geo_load(nurbs);


        geometry.nurbs.controlPoints = zeros(geometry.nurbs.number(1),geometry.nurbs.number(2),geometry.nurbs.number(3),3);
        geometry.nurbs.controlPoints(:,:,:,1) = reshape(geometry.nurbs.coefs(1,:,:,:), geometry.nurbs.number);
        geometry.nurbs.controlPoints(:,:,:,2) = reshape(geometry.nurbs.coefs(2,:,:,:), geometry.nurbs.number);
        geometry.nurbs.controlPoints(:,:,:,3) = reshape(geometry.nurbs.coefs(3,:,:,:), geometry.nurbs.number);
        geometry.tensor.controlPoints = reshape(geometry.nurbs.controlPoints, prod(geometry.nurbs.number),3);

        time = toc(time);

        stats = struct; 
        stats.time_prep = time; 
        stats.time = time;

        [H, ~, low_rank_data, stats_interpol] = interpolate_weights_bsplines(geometry, low_rank_data);
        
        stats.time_interpol = stats_interpol.time; 
        stats.time = stats.time + stats_interpol.time; 
        stats.A_memory = stats_interpol.A_memory; 
        stats.A_nnz = stats_interpol.A_nnz; 
        stats.b_memory = stats_interpol.b_memory; 
        stats.b_nnz = stats_interpol.b_nnz; 


        [H, time] = lowRank_w(H, low_rank_data);

        stats.time_format = time;
        stats.time = stats.time + stats.time_format; 

    else
        degelev = max (low_rank_data.system_degree - (geometry.nurbs.order-1), 0);
        nurbs = nrbdegelev (geometry.nurbs, degelev);

        [~, ~, new_knots] = kntrefine (nurbs.knots, low_rank_data.system_nsub, ...
            low_rank_data.system_degree, low_rank_data.system_regularity);

        nurbs = nrbkntins(nurbs, new_knots);
        geometry = geo_load(nurbs);

        geometry.nurbs.weight = reshape(geometry.nurbs.coefs(4,:,:,:),geometry.nurbs.number);
        geometry.nurbs.controlPoints = zeros(geometry.nurbs.number(1),geometry.nurbs.number(2),geometry.nurbs.number(3),3);
        geometry.nurbs.controlPoints(:,:,:,1) = reshape(geometry.nurbs.coefs(1,:,:,:)./geometry.nurbs.coefs(4,:,:,:), geometry.nurbs.number);
        geometry.nurbs.controlPoints(:,:,:,2) = reshape(geometry.nurbs.coefs(2,:,:,:)./geometry.nurbs.coefs(4,:,:,:), geometry.nurbs.number);
        geometry.nurbs.controlPoints(:,:,:,3) = reshape(geometry.nurbs.coefs(3,:,:,:)./geometry.nurbs.coefs(4,:,:,:), geometry.nurbs.number);
        weightR = reshape(geometry.nurbs.weight, geometry.nurbs.number(1), geometry.nurbs.number(2), geometry.nurbs.number(3));
        geometry.tensor.Tweights = round(tt_tensor(weightR), 1e-15, 1);
        geometry.tensor.weight = kron(kron(geometry.tensor.Tweights{3},geometry.tensor.Tweights{2}),geometry.tensor.Tweights{1})';
        geometry.tensor.controlPoints = reshape(geometry.nurbs.controlPoints, prod(geometry.nurbs.number),3);

        time = toc(time);

        stats = struct; 
        stats.time_prep = time; 
        stats.time = time;

        [H, ~, low_rank_data, stats_interpol] = interpolate_weights_nurbs(geometry, low_rank_data);
        
        stats.time_interpol = stats_interpol.time; 
        stats.time = stats.time + stats_interpol.time; 
        stats.A_memory = stats_interpol.A_memory; 
        stats.A_nnz = stats_interpol.A_nnz; 
        stats.b_memory = stats_interpol.b_memory; 
        stats.b_nnz = stats_interpol.b_nnz; 


        [H, time] = lowRank_w(H, low_rank_data);

        stats.time_format = time;
        stats.time = stats.time + stats.time_format; 


    end


end

function [H, rhs, opt, stats] = interpolate_weights_bsplines(G, opt)
% INTERPOLATE_WEIGHTS_BSPLINES  Low-rank interpolation of geometry weights on B-spline geometries.
%
%   [H, RHS, OPT] = INTERPOLATE_WEIGHTS_BSPLINES(G, OPT)
%
%   Purpose
%   -------
%   Builds separable (tensor) interpolants of the geometry-induced weights used in
%   univariate quadrature for mass and/or stiffness assembly when the geometry G
%   is (non-rational) B-spline. In 2D it interpolates the entries of
%       Q(s,t) = |det(JG(s,t))| * inv(JG(s,t)) * inv(JG(s,t))',
%   while in 3D it interpolates the six unique entries of the symmetric matrix
%   associated with inv(JG)*inv(JG)' scaled by |det(JG)|. Coefficients are
%   computed from values on a Greville grid in an enlarged target spline space.
%
%   Inputs
%   ------
%   G        Geometry struct with fields:
%              .rdim               spatial dimension (2 or 3)
%              .nurbs.knots{i}     open knot vectors (i = 1..rdim)
%              .nurbs.order(i)     spline orders (degree = order-1)
%              .tensor.controlPoints
%                                   tensorized control points:
%                                     - 2D: size [dimS,dimT,2] (used via matrix multiplies)
%                                     - 3D: reshaped to [prod(dim), 3] for TT ops
%
%   OPT      (optional, struct) controls interpolation/solvers. Fields:
%              .mass          (0/1) build mass weights  (default 0)
%              .stiffness     (0/1) build stiffness weights (default 0)
%              .greville      scale for Greville abscissae (1 = standard, default 1)
%              .plotW         reserved (ignored here; default 0)
%              .TT_interpolation
%                             3D only: if 1, solve in TT via AMEn; if 0, dense Kron solve.
%                             (If stiffness requested and not set, it is forced to 1.)
%              .rankTol       TT rounding/accuracy tolerance used in TT mode (required when TT_interpolation=1)
%              .splinespace2, .splinedegree2
%                             reserved (not used; kept for compatibility)
%
%   Outputs
%   -------
%   H        Struct with fields:
%              .dim           = G.rdim
%              .weightFun     description of target interpolation space:
%                               .knots{i}, .degree(i), .n(i)  (i = 1..rdim)
%              .mass.weightMat
%                             coefficients of |det(JG)| in the target space:
%                               - 2D:  matrix [n1 x n2]
%                               - 3D:  TT tensor (if TT_interpolation=1) or full [n1 x n2 x n3]
%              .stiffness.weightMat
%                             coefficients of the symmetric matrix entries:
%                               - 2D:  array [n1 x n2 x 3] storing (11,12,22); 21=12
%                               - 3D:  cell(6,1) storing (11,12,13,22,23,33) as TT or full arrays.
%                                     Entries that are numerically negligible (mean < rankTol)
%                                     may be left empty to save work.
%
%   RHS      Struct carrying the mass weight in the same representation as H.mass.weightMat.
%            In 3D:
%              * if OPT.mass==1, RHS.weightMat == H.mass.weightMat (after rounding).
%              * if OPT.mass==0, RHS.weightMat still stores the TT interpolation of |det(JG)|
%                for later reuse. In 2D, RHS is unused.
%
%   OPT      Echoed back. In 3D stiffness mode, OPT.TT_interpolation may be set to 1 internally.
%
%   Method (what the code does)
%   ---------------------------
%   1) Target space: for each parametric dir i, create an enlarged B-spline space
%      (knots/degree/size) via ENLARGEN_BSPLINE_SPACE_W(..., degree, 3). Store in H.weightFun.
%   2) Greville grid: build (scaled) Greville points in each dir and evaluate:
%         - original basis and derivatives (for JG),
%         - target basis (for interpolation system).
%   3) Evaluate Jacobian terms on the Greville tensor grid:
%         2D:   j11=jac1, j12=jac3, j21=jac2, j22=jac4; w = |j11*j22 - j12*j21|.
%               Form q11,q12,q22 for stiffness and/or w for mass, then solve
%               (Kronecker system) to get coefficient arrays.
%         3D:   assemble jac_αβ via TT-matrix times control points; compute
%               w = |det(JG)| from 3×3 minors. Interpolate w (mass) and the six
%               stiffness combinations using either:
%                 • TT AMEn (amen_block_solve with nswp=20, kickrank=2, resid_damp=1e1),
%                   then round to OPT.rankTol; or
%                 • a single dense Kron solve M\vec, reshaped to [n1 n2 n3].
%
%   Notes
%   ------------
%   * Set at least one of OPT.mass or OPT.stiffness to 1; otherwise the function errors.
%   * For large 3D problems, use TT_interpolation=1 with a sensible OPT.rankTol (e.g., 1e-10…1e-6).
%   * The target space (H.weightFun.*) is typically richer than the geometry space to
%     capture metric variations; change ENLARGEN_BSPLINE_SPACE_W if you want a different policy.
%   * In 3D stiffness, components with tiny average magnitude (relative to rankTol) are skipped.
%
%   Example (2D)
%   -----------
%     G = your_bspline_geometry();                % GeoPDEs-style struct
%     opt.mass = 1;  opt.stiffness = 1;           % both weights
%     [H, rhs, opt] = interpolate_weights_bsplines(G, opt);
%     % H.mass.weightMat is [n1 x n2]; H.stiffness.weightMat(:,:,1/2/3) = Q11/Q12/Q22
%
%   Example (3D, TT mode)
%   ---------------------
%     G = your_bspline_geometry_3d();
%     opt.mass = 1; opt.stiffness = 1;
%     opt.TT_interpolation = 1; opt.rankTol = 1e-8;
%     [H, rhs, opt] = interpolate_weights_bsplines(G, opt);
%     % H.mass.weightMat and H.stiffness.weightMat{k} are TT tensors; rhs.weightMat reused later.
%
%   See also
%   --------
%   ENLARGEN_BSPLINE_SPACE_W, GENERATEGREVILLEPOINTS, EVALBSPLINE, EVALBSPLINEDERIV,
%   TT_TENSOR, AMEN_BLOCK_SOLVE, TT_MATRIX, KRON.

    time = tic;

    if nargin < 2
        opt = struct();
    end
    if ~isfield(opt, 'plotW') || isempty(opt.plotW)
        opt.plotW = 0;
    end
    if ~isfield(opt, 'stiffness') || isempty(opt.stiffness)
        opt.stiffness = 0;
    end
    if ~isfield(opt, 'mass') || isempty(opt.mass)
        opt.mass = 0;
    end
    if ~isfield(opt,'greville')
        opt.greville = 1;
    end
    if opt.mass == 0 && opt.stiffness == 0
        error('Please specify a system matrix to be computed');
    end
    
    if ~isfield(opt,'splinespace2') || isempty(opt.splinespace2)
        opt.splinespace2 = 0;
    end
    
    if ~isfield(opt,'splinedegree2') || isempty(opt.splinedegree2)
        opt.splinedegree2 = 0;
    end
    
    H = struct;
    H.weightFun = struct;
    H.dim = G.rdim;
    
    rhs = struct;
    
    
    H.weightFun.knots = cell(G.rdim,1);
    H.weightFun.degree = zeros(G.rdim,1);
    H.weightFun.n = zeros(G.rdim,1);

    use_aux_space = isfield(opt, 'aux_space') && ~isempty(opt.aux_space);

    if use_aux_space
        required_fields = {'knots', 'degree', 'ndof_dir'};
    
        for k = 1:numel(required_fields)
            if ~isfield(opt.aux_space, required_fields{k})
                error('interpolate_weights_bsplines:InvalidAuxSpace', ...
                    'opt.aux_space must contain the field "%s".', ...
                    required_fields{k});
            end
        end
    
        if numel(opt.aux_space.knots) ~= G.rdim || ...
                numel(opt.aux_space.degree) ~= G.rdim || ...
                numel(opt.aux_space.ndof_dir) ~= G.rdim
            error('interpolate_weights_bsplines:InvalidAuxSpaceDimension', ...
                'The auxiliary space must contain one entry per parametric direction.');
        end
    
        for i = 1:G.rdim
            H.weightFun.knots{i} = opt.aux_space.knots{i};
            H.weightFun.degree(i) = opt.aux_space.degree(i);
    
            % Recompute the number of basis functions from knots and degree.
            % This avoids inconsistencies in a manually constructed aux_space.
            H.weightFun.n(i) = ...
                numel(H.weightFun.knots{i}) ...
                - H.weightFun.degree(i) - 1;
    
            if H.weightFun.n(i) ~= opt.aux_space.ndof_dir(i)
                error('interpolate_weights_bsplines:InconsistentAuxSpace', ...
                    ['In direction %d, aux_space.ndof_dir is %d, but the ', ...
                     'knot vector and degree imply %d basis functions.'], ...
                    i, opt.aux_space.ndof_dir(i), H.weightFun.n(i));
            end
        end
    
    else
        if ~isfield(opt, 'p_ref') || isempty(opt.p_ref)
            opt.p_ref = zeros(1, G.rdim);
        elseif isscalar(opt.p_ref)
            opt.p_ref = repmat(opt.p_ref, 1, G.rdim);
        elseif numel(opt.p_ref) ~= G.rdim
            error('interpolate_weights_bsplines:InvalidPRef', ...
                'opt.p_ref must be scalar or have one entry per direction.');
        end
    
        for i = 1:G.rdim
            [H.weightFun.knots{i}, ...
             H.weightFun.degree(i), ...
             H.weightFun.n(i)] = ...
                enlargen_bspline_space( ...
                    G.nurbs.knots{i}, ...
                    G.nurbs.order(i) - 1, ...
                    3);
    
            if opt.p_ref(i) > 0
                [H.weightFun.knots{i}, ...
                 H.weightFun.degree(i)] = ...
                    elevate_knot_vector( ...
                        H.weightFun.knots{i}, ...
                        H.weightFun.degree(i), ...
                        opt.p_ref(i));
    
                H.weightFun.n(i) = ...
                    numel(H.weightFun.knots{i}) ...
                    - H.weightFun.degree(i) - 1;
            end
        end
    end
    
    
    
    grevillePoints = cell(G.rdim,1);
    grevilleValues = cell(G.rdim,1);
    grevilleDerivs = cell(G.rdim,1);
    grevilleValues2 = cell(G.rdim,1);
    
    %% Generate Greville points and evaluate basis splines there
    for i = 1:G.rdim
        grevillePoints{i} = opt.greville*generateGrevillePoints(H.weightFun.knots{i}, H.weightFun.degree(i));     
        % Compute basis spline values on greville points
        grevilleValues{i} =  sparse(evalBSpline(G.nurbs.knots{i}, G.nurbs.order(i)-1, grevillePoints{i}));
        grevilleDerivs{i} =  sparse(evalBSplineDeriv(G.nurbs.knots{i}, G.nurbs.order(i)-1, grevillePoints{i}));
        grevilleValues2{i} = sparse(evalBSpline(H.weightFun.knots{i}, H.weightFun.degree(i), grevillePoints{i}));
    end
    
    %% Setup equation matrix and right hand sides
    
    if G.rdim == 2
        % Each row in grevillePoints2D contains [s t] coordinates of a
        % combination of a point in grevillePointsS and one in grevillePointsT.
        % The number of rows is dimS2*dimT2.
        M = kron(grevilleValues2{2}',grevilleValues2{1}');
        jac1 = grevilleDerivs{1}' * G.tensor.controlPoints(:,:,1) * grevilleValues{2};
        jac2 = grevilleDerivs{1}' * G.tensor.controlPoints(:,:,2) * grevilleValues{2};
        jac3 = grevilleValues{1}' * G.tensor.controlPoints(:,:,1) * grevilleDerivs{2};
        jac4 = grevilleValues{1}' * G.tensor.controlPoints(:,:,2) * grevilleDerivs{2};
        w = abs(jac1.*jac4-jac3.*jac2);
        if opt.stiffness == 1
            q1 = reshape(1./w .* (jac4.^2 + jac3.^2),prod(H.weightFun.n),1);
            q2 = reshape(-1./w .*(jac2.*jac4 + jac1.*jac3), prod(H.weightFun.n),1);
            q3 = reshape(1./w .*(jac1.^2 + jac2.^2), prod(H.weightFun.n),1);
            vecWeights = M \ [q1,q2,q3];
            H.stiffness.weightMat = reshape(vecWeights, H.weightFun.n(1), H.weightFun.n(2), 3);
        end
        if opt.mass == 1
            w1 = reshape(w, H.weightFun.n(1)*H.weightFun.n(2),1);
            H.mass.weightMat = reshape(M\w1, H.weightFun.n(1),H.weightFun.n(2));
        end
    elseif G.rdim == 3
        nswp = 20;
    
        A = tt_matrix({grevilleDerivs{1}'; grevilleValues{2}'; grevilleValues{3}'});
        jac11 = A*G.tensor.controlPoints(:,1);
        jac12 = A*G.tensor.controlPoints(:,2);
        jac13 = A*G.tensor.controlPoints(:,3);
        A = tt_matrix({grevilleValues{1}'; grevilleDerivs{2}'; grevilleValues{3}'});
        jac21 = A*G.tensor.controlPoints(:,1);
        jac22 = A*G.tensor.controlPoints(:,2);
        jac23 = A*G.tensor.controlPoints(:,3);
        A = tt_matrix({grevilleValues{1}'; grevilleValues{2}'; grevilleDerivs{3}'});
        jac31 = A*G.tensor.controlPoints(:,1);
        jac32 = A*G.tensor.controlPoints(:,2);
        jac33 = A*G.tensor.controlPoints(:,3);
        clear A;
        w = abs(jac11.*jac22.*jac33 + jac12.*jac23.*jac31 + jac21.*jac32.*jac13 - jac13.*jac22.*jac31 - jac12.*jac21.*jac33 - jac11.*jac23.*jac32);   
        
        m_1 = size(grevilleValues2{1}',1);
        m_2 = size(grevilleValues2{2}',1);
        m_3 = size(grevilleValues2{3}',1);
    
        
    
        if opt.mass == 1
            if isfield(opt, 'TT_interpolation') && opt.TT_interpolation == 1
                MM = {grevilleValues2{1}';grevilleValues2{2}';grevilleValues2{3}'};
                tt_w = tt_tensor(reshape(w, [m_1,m_2,m_3]), 1e-16);
                H.mass.weightMat= amen_block_solve({MM},{tt_w},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                clear tt_w;
                H.mass.weightMat = round(H.mass.weightMat, opt.rankTol);
            else
                M= kron(grevilleValues2{3},kron(grevilleValues2{2},grevilleValues2{1}))';
                vecWeights = M\w; 
                H.mass.weightMat = reshape(vecWeights, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                clear vecWeights;
            end
    
            
            rhs.weightMat = H.mass.weightMat;
            
        end
    
        w = 1./w;
    
        if opt.stiffness == 1
            b_memory = 0;
            b_nnz = 0;
    
            % if isfield(opt, 'TT_interpolation') && opt.TT_interpolation == 0
            %     M= kron(grevilleValues2{3},kron(grevilleValues2{2},grevilleValues2{1}))';
            %     MM = tt_matrix(M); 
            %     H.stiffness.weightMat = cell(6,1);
            % else
            %     MM = {grevilleValues2{1}';grevilleValues2{2}';grevilleValues2{3}'};
            %     H.stiffness.weightMat = cell(6,1);
            %     opt.TT_interpolation = 1;
            % end
    
            H.stiffness.weightMat = cell(6,1);
            MM = {grevilleValues2{1}';grevilleValues2{2}';grevilleValues2{3}'};
    
            m11 = jac22.*jac33 - jac32.*jac23;
            m21 = jac21.*jac33 - jac31.*jac23;
            m31 = jac21.*jac32 - jac31.*jac22;
            q11 = w.*(m11.^2 + m21.^2 + m31.^2);
            if sum(abs(q11))/numel(q11) > opt.rankTol
                if opt.TT_interpolation == 1
                    tt_11 = tt_tensor(reshape(q11, [m_1,m_2,m_3]), 1e-16);
                    H.stiffness.weightMat{1} = amen_block_solve({MM},{tt_11},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                    clear tt_11;
                else
                    M = kron(grevilleValues2{3},kron(grevilleValues2{2},grevilleValues2{1}))';
                    H.stiffness.weightMat{1} = reshape(M\q11, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                end
            end
            b_memory = b_memory + RecursiveSize(q11);
            b_nnz = b_nnz + RecursiveSize(q11);
            clear q11;
    
    
            m12 = jac12.*jac33 - jac32.*jac13;
            m22 = jac11.*jac33 - jac31.*jac13;
            m32 = jac11.*jac32 - jac31.*jac12;
            q12 = w.*(-m11.*m12 - m21.*m22 - m31.*m32);
            if sum(abs(q12))/numel(q12) > opt.rankTol
                if opt.TT_interpolation == 1
                    tt_12 = tt_tensor(reshape(q12, [m_1,m_2,m_3]), 1e-16);
                    H.stiffness.weightMat{2} = amen_block_solve({MM},{tt_12},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                    clear tt_12;
                else
                    H.stiffness.weightMat{2} = reshape(M\q12, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                end
            end
            b_memory = b_memory + RecursiveSize(q12);
            b_nnz = b_nnz + RecursiveSize(q12);
            clear q12;
            clear jac31;
            clear jac32;
            clear jac33;
    
    
            m13 = jac12.*jac23 - jac22.*jac13;
            m23 = jac11.*jac23 - jac21.*jac13;
            m33 = jac11.*jac22 - jac21.*jac12;
            q13 = w.*(m11.*m13 + m21.*m23 + m31.*m33);
            if sum(abs(q13))/numel(q13) > opt.rankTol
                if opt.TT_interpolation == 1
                    tt_13 = tt_tensor(reshape(q13, [m_1,m_2,m_3]), 1e-16);
                    H.stiffness.weightMat{3} = amen_block_solve({MM},{tt_13},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                    clear tt_13;
                else
                    H.stiffness.weightMat{3} = reshape(M\q13, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                end
            end
            b_memory = b_memory + RecursiveSize(q13);
            b_nnz = b_nnz + RecursiveSize(q13);
            clear q13; 
            clear m11;
            clear m21;
            clear m31;
            clear jac11;
            clear jac12;
            clear jac13;
            clear jac21;
            clear jac22;
            clear jac23;
    
    
            q22 = w.*(m12.^2 + m22.^2 + m32.^2);
            if sum(abs(q22))/numel(q22) > opt.rankTol
                if opt.TT_interpolation == 1
                    tt_22 = tt_tensor(reshape(q22, [m_1,m_2,m_3]), 1e-16);
                    H.stiffness.weightMat{4} = amen_block_solve({MM},{tt_22},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                    clear tt_22;
                else
                    H.stiffness.weightMat{4} = reshape(M\q22, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                end
            end
            b_memory = b_memory + RecursiveSize(q22);
            b_nnz = b_nnz + RecursiveSize(q22);
            clear q22;
    
    
            q23 = w.*(-m12.*m13 - m22.*m23 - m32.*m33);
            if sum(abs(q23))/numel(q23) > opt.rankTol
                if opt.TT_interpolation == 1
                    tt_23 = tt_tensor(reshape(q23, [m_1,m_2,m_3]), 1e-16);
                    H.stiffness.weightMat{5} = amen_block_solve({MM},{tt_23},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                    clear tt_23;
                else
                    H.stiffness.weightMat{5} = reshape(M\q23, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                end
            end     
            b_memory = b_memory + RecursiveSize(q23);
            b_nnz = b_nnz + RecursiveSize(q23);
            clear q23;
            clear m12;
            clear m22;
            clear m32;
    
    
            q33 = w.*(m13.^2 + m23.^2 + m33.^2);
            if sum(abs(q33))/numel(q33) > opt.rankTol
                if opt.TT_interpolation == 1
                    tt_33 = tt_tensor(reshape(q33, [m_1,m_2,m_3]), 1e-16);
                    H.stiffness.weightMat{6} = amen_block_solve({MM},{tt_33},opt.rankTol, 'kickrank', 2, 'resid_damp', 1e1, 'nswp', nswp, 'exitdir', -1);
                    clear tt_33;
                else
                    H.stiffness.weightMat{6} = reshape(M\q33, H.weightFun.n(1), H.weightFun.n(2), H.weightFun.n(3));
                end
            end
            b_memory = b_memory + RecursiveSize(q33);
            b_nnz = b_nnz + RecursiveSize(q33);
            clear q33;
            clear m13;
            clear m23;
            clear m33;
    
        end
    
        % clear w;
        clear M; 
        % clear MM;
        
    end

    time = toc(time);

    stats = struct; 
    stats.time = time; 

    if opt.mass == 1
        stats.b_memory = RecursiveSize(w);
        stats.b_nnz = nnz(w);
        stats.A_memory = RecursiveSize(MM);
        stats.A_nnz = nnz(grevilleValues2{1}) + nnz(grevilleValues2{2}) + nnz(grevilleValues2{3});
    end

    if opt.stiffness == 1
        stats.b_memory = b_memory;
        stats.b_nnz = b_nnz;
        stats.A_memory = RecursiveSize(MM);
        stats.A_nnz = nnz(grevilleValues2{1}) + nnz(grevilleValues2{2}) + nnz(grevilleValues2{3});
    end




end


function [knots_det, degree_det, n_det] = enlargen_bspline_space(knots, degree, dim)
%ENLARGEN_BSPLINE_SPACE Construct the spline space containing det(grad G).
%
% Let the original univariate geometry space have degree p and let an
% interior knot have multiplicity mu. For a D-dimensional geometry, every
% term in det(grad G), considered in one parametric direction, has:
%
%   degree:                p_det  = D*p - 1,
%   interior multiplicity: mu_det = (D-1)*p + mu.
%
% The resulting open knot vector therefore contains det(grad G) exactly.
%
% Inputs:
%   knots   Open, nondecreasing knot vector.
%   degree  Original geometry degree p.
%   dim     Parametric dimension D.
%
% Outputs:
%   knots_det   Knot vector of the determinant space.
%   degree_det  Degree D*p - 1.
%   n_det       Number of basis functions.

    validateattributes(knots, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nondecreasing'});
    validateattributes(degree, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(dim, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    knots = knots(:).';

    left_endpoint = knots(1);
    right_endpoint = knots(end);

    [unique_knots, ~, knot_ids] = unique(knots);
    multiplicities = accumarray(knot_ids(:), 1).';

    degree_det = dim * degree - 1;

    % Open left endpoint: multiplicity degree_det + 1 = dim*degree.
    knots_det = repmat( ...
        left_endpoint, 1, degree_det + 1);

    % Interior knots.
    for i = 2:numel(unique_knots)-1
        multiplicity_det = ...
            (dim - 1) * degree + multiplicities(i);

        knots_det = [knots_det, ...
            repmat(unique_knots(i), ...
                   1, multiplicity_det)]; 
    end

    % Open right endpoint.
    knots_det = [knots_det, ...
        repmat(right_endpoint, 1, degree_det + 1)];

    n_det = numel(knots_det) - degree_det - 1;
end

function [newKnots, newDegree] = elevate_knot_vector(knots, degree, degreeIncrease)
%ELEVATE_KNOT_VECTOR Construct the knot vector for degree elevation.
%
%   [newKnots, newDegree] = elevate_knot_vector(knots, degree, degreeIncrease)
%
% increases the spline degree from
%
%       degree
%
% to
%
%       newDegree = degree + degreeIncrease
%
% while preserving the continuity at every knot. Consequently, the
% multiplicity of every distinct knot is increased by degreeIncrease.
%
% Inputs:
%   knots          - Nondecreasing univariate knot vector.
%   degree         - Original spline degree.
%   degreeIncrease - Positive integer by which the degree is elevated.
%
% Outputs:
%   newKnots       - Knot vector associated with the elevated degree.
%   newDegree      - Elevated spline degree.
%
% Example:
%   knots = [0 0 0 0.5 1 1 1];
%   degree = 2;
%   degreeIncrease = 3;
%
%   [newKnots, newDegree] = elevate_knot_vector( ...
%       knots, degree, degreeIncrease);
%
% This gives newDegree = 5 and
%
%   newKnots =
%       [0 0 0 0 0 0 0.5 0.5 0.5 0.5 1 1 1 1 1 1]

    validateattributes(knots, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nonempty'}, ...
        mfilename, 'knots', 1);

    validateattributes(degree, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'}, ...
        mfilename, 'degree', 2);

    validateattributes(degreeIncrease, {'numeric'}, ...
        {'scalar', 'integer', 'positive'}, ...
        mfilename, 'degreeIncrease', 3);

    % Preserve whether the input was a row or column vector.
    inputIsColumn = iscolumn(knots);

    knots = knots(:).';

    if any(diff(knots) < 0)
        error('elevate_knot_vector:InvalidKnotVector', ...
            'The knot vector must be nondecreasing.');
    end

    newDegree = degree + degreeIncrease;

    % Determine the distinct knots and their multiplicities.
    runStart = [1, find(diff(knots) ~= 0) + 1];
    runEnd   = [runStart(2:end) - 1, numel(knots)];

    distinctKnots = knots(runStart);
    multiplicities = runEnd - runStart + 1;

    % Preserving C^(degree - multiplicity) regularity requires
    %
    %   newMultiplicity = multiplicity + degreeIncrease.
    newMultiplicities = multiplicities + degreeIncrease;

    knotBlocks = arrayfun( ...
        @(value, multiplicity) repmat(value, 1, multiplicity), ...
        distinctKnots, newMultiplicities, ...
        'UniformOutput', false);

    newKnots = [knotBlocks{:}];

    if inputIsColumn
        newKnots = newKnots.';
    end
end

function points = generateGrevillePoints(knots, degree)
% GENERATEGREVILLEPOINTS  Compute Greville abscissae for a univariate B-spline space.
%
%   POINTS = GENERATEGREVILLEPOINTS(KNOTS, DEGREE)
%
%   Purpose
%   -------
%   Returns the Greville abscissae (collocation/interpolation points) associated with the
%   univariate B-spline space S(KNOTS, DEGREE). These points are used to sample functions
%   for interpolation in your low-rank setup.
%
%   Inputs
%   ------
%   KNOTS    Nondecreasing knot vector (row or column). Typically open/clamped on [0,1].
%   DEGREE   Polynomial degree p (must satisfy p >= 1).
%
%   Output
%   -------
%   POINTS   Row vector of length dim = numel(KNOTS) - DEGREE - 1 containing the Greville
%            abscissae:
%               POINTS(i) = (KNOTS(i+1) + KNOTS(i+2) + ... + KNOTS(i+DEGREE)) / DEGREE,
%            for i = 1, …, dim.
%
%   Details
%   -------
%   * dim equals the number of B-spline basis functions in S(KNOTS, DEGREE).
%   * For open/clamped KNOTS, POINTS lie in [KNOTS(1), KNOTS(end)] and increase
%     monotonically (with repeats at multiple knots).
%   * Repeated interior knots cause clusters of repeated Greville points, reflecting
%     reduced continuity there.
%
%   Notes
%   -----
%   * DEGREE must be at least 1 (the formula divides by DEGREE).
%   * The function returns a row vector. Use POINTS(:) if you prefer a column vector.
%   * Works for any nondecreasing KNOTS; open/clamped is recommended for standard IgA setups.
%
%   Example
%   -------
%     % Quadratic (p=2) open knot vector with simple interior knots:
%     K = [0 0 0  0.25  0.5  0.75  1 1 1];
%     p = 2;
%     pts = generateGrevillePoints(K, p)
%     % pts(i) = (K(i+1) + K(i+2)) / 2,  i = 1..numel(K)-p-1
%
%   See also
%   --------
%   EVALBSPLINE, EVALNURBS, ENLARGEN_BSPLINE_SPACE, KNTREFINE, NRBDEGELEV.

    dim = numel(knots) - degree - 1;
    
    points = zeros(1, dim);

    for i = 1:dim
        points(i) = sum(knots(i+1:i+degree)) / (degree);
    end
end


function splines = evalBSpline(knots, degree, eta)
% EVALBSPLINE  Evaluate all univariate B-spline basis functions at given points.
%
%   SPL = EVALBSPLINE(KNOTS, DEGREE, ETA)
%
%   Purpose
%   -------
%   Computes the values of every B-spline basis function N_{i,p}(η) from the space
%   S(KNOTS, p=DEGREE) at the query points ETA. The result is returned as a matrix
%   with one row per basis function and one column per evaluation point.
%
%   Inputs
%   ------
%   KNOTS    Nondecreasing knot vector (row or column). Typically open/clamped.
%   DEGREE   Polynomial degree p (p >= 0).
%   ETA      Evaluation points (row or column vector). Values are typically within
%            [KNOTS(1), KNOTS(end)].
%
%   Output
%   -------
%   SPL      Matrix of size [m x n], where
%              m = numel(KNOTS) - DEGREE - 1   (number of basis functions),
%              n = numel(ETA).
%            Entry SPL(i,j) = N_{i,p}(ETA(j)).
%
%   Details
%   -------
%   * Uses the Cox–de Boor recursion:
%       For p = 0:
%         N_{i,0}(η) = 1  if  KNOTS(i) ≤ η < KNOTS(i+1),  else 0.
%         To ensure left-continuity at the last nonempty span, the code sets
%         SPL(i_end, η == KNOTS(end)) = 1.
%       For p > 0:
%         N_{i,p}(η) = (η - KNOTS(i)) / (KNOTS(i+p)   - KNOTS(i  )) * N_{i,  p-1}(η)
%                    + (KNOTS(i+p+1) - η) / (KNOTS(i+p+1) - KNOTS(i+1)) * N_{i+1,p-1}(η),
%         with each fraction skipped when its denominator is zero (knot multiplicity).
%
%   Conventions & boundary handling
%   -------------------------------
%   * Degree-0 basis functions are 1 on half-open intervals [KNOTS(i), KNOTS(i+1)).
%   * At the right endpoint η = KNOTS(end), the last active basis is set to 1 to
%     preserve partition of unity for open/clamped knot vectors.
%
%   Performance notes
%   -----------------
%   * Complexity is O(m * p * n). The routine is vectorized over ETA inside each
%     degree loop. For very large problems, consider making SPL sparse:
%       % SPL = sparse(SPL);   % optional, basis is locally supported
%
%   Example
%   -------
%     % Quadratic (p=2) open knot vector with simple interior knots:
%     K = [0 0 0  0.25  0.5  0.75  1 1 1];
%     p = 2;
%     xi = linspace(0,1,11);
%     N = evalBSpline(K, p, xi);   % size(N) == [numel(K)-p-1, numel(xi)]
%     % Check partition of unity:
%     max(abs(sum(N,1) - 1))   % ~ 0 up to roundoff (including xi(end)=1)
%
%   See also
%   --------
%   EVALNURBS, EVALNURBSDERIV, GENERATEGREVILLEPOINTS.

    m = numel(knots) - degree - 1;
    
    splines = zeros(m, length(eta));
    if degree == 0
        for i = 1:m
            splines(i, :) = (knots(i) <= eta) & (eta < knots(i+1));
        end
        if knots(1) < knots(end)
            i = find(knots < knots(end), 1, 'last');
            splines(i, eta == knots(end)) = 1;
        end
    else
        prevSplines = evalBSpline(knots, degree-1, eta);
        for i = 1:m
            if knots(i+degree) > knots(i)
                splines(i, :) = (eta - knots(i)) / (knots(i+degree) - knots(i)) .* prevSplines(i, :);
            end
            if knots(i+degree+1) > knots(i+1)
                splines(i, :) = splines(i, :) + (knots(i+degree+1) - eta) / (knots(i+degree+1) - knots(i+1)) .* prevSplines(i+1, :);
            end
        end
        splines(m, eta == knots(end)) = 1;
    end
end
   
function derivatives = evalBSplineDeriv(knots, degree, eta)
% EVALBSPLINEDERIV  Evaluate first derivatives of univariate B-spline basis functions.
%
%   dN = EVALBSPLINEDERIV(KNOTS, DEGREE, ETA)
%
%   Purpose
%   -------
%   Computes d/dη N_{i,p}(η) for all basis functions in the B-spline space
%   S(KNOTS, p=DEGREE) at the query points ETA. Returns one row per basis
%   function and one column per evaluation point.
%
%   Inputs
%   ------
%   KNOTS    Nondecreasing knot vector (row or column). Typically open/clamped.
%   DEGREE   Polynomial degree p (integer, p >= 0).
%   ETA      Evaluation points (row or column vector).
%
%   Output
%   -------
%   dN       Matrix of size [m x n], where
%              m = numel(KNOTS) - DEGREE - 1   (number of basis functions),
%              n = numel(ETA).
%            Entry dN(i,j) = d/dη N_{i,p}(ETA(j)).
%
%   Details
%   -------
%   * Uses the Cox–de Boor derivative formula (with zero-denominator guards):
%       For p = 0:   dN_{i,0}(η) = 0.
%       For p > 0:   dN_{i,p}(η) =
%           p/(KNOTS(i+p)   - KNOTS(i  )) * N_{i,  p-1}(η)  ...
%         - p/(KNOTS(i+p+1) - KNOTS(i+1)) * N_{i+1,p-1}(η),
%     with each fraction omitted when its denominator is zero (due to knot multiplicity).
%   * Internally evaluates the degree-(p-1) basis via EVALBSPLINE and combines terms.
%
%   Properties / sanity checks
%   --------------------------
%   * Partition of unity derivative:  sum_i dN_{i,p}(η) = 0  for all η.
%   * On open spans, dN is continuous up to C^{p-2}; at multiple knots, continuity
%     drops accordingly.
%
%   Example
%   -------
%     % Quadratic (p=2) open knot vector with simple interior knots:
%     K = [0 0 0  0.25  0.5  0.75  1 1 1];  p = 2;
%     xi = linspace(0,1,11);
%     dN = evalBSplineDeriv(K, p, xi);          % size(dN) = [(numel(K)-p-1) x 11]
%     max(abs(sum(dN,1)))                        % ~ 0 (partition-of-unity derivative)
%
%   See also
%   --------
%   EVALBSPLINE, EVALNURBS, EVALNURBSDERIV, GENERATEGREVILLEPOINTS.


    m = numel(knots) - degree - 1;

    derivatives = zeros(m, length(eta));
    
    if degree == 0
        % all derivatives zero, do nothing
    else
        prevSplines = evalBSpline(knots, degree-1, eta);
        for i = 1:m
            if knots(i+degree) > knots(i)
                derivatives(i, :) = degree / (knots(i+degree) - knots(i)) .* prevSplines(i, :);
            end
            if knots(i+degree+1) > knots(i+1)
                derivatives(i, :) = derivatives(i, :) - degree / (knots(i+degree+1) - knots(i+1)) .* prevSplines(i+1, :);
            end
        end
    end
end
    
function [H, time] = lowRank_w(H, opt)
% LOWRANK_W  Convert/interleave weight tensors into 1D factors for fast univariate integration.
%
%   [H, OPT] = LOWRANK_W(H, OPT)
%
%   Purpose
%   -------
%   Takes the weight data produced by INTERPOLATE_WEIGHTS_* (mass |det(JF)| and/or stiffness
%   entries) and converts them into separable 1D factors that are convenient for univariate
%   quadrature. In 2D this is done via SVD of each coefficient matrix; in 3D it reshapes TT cores
%   into per-direction factor blocks. Optional truncation (rank selection) is performed using a
%   tolerance, and (to save memory) the original full/TT objects can be discarded.
%
%   Inputs
%   ------
%   H    Struct returned from INTERPOLATE_WEIGHTS_BSPLINES or *_NURBS with fields:
%          .dim                = 2 or 3
%          .weightFun.n        vector of target-space sizes per direction
%          .mass.weightMat     (2D: full [n1 x n2]; 3D: TT tensor)  — if OPT.mass==1
%          .stiffness.weightMat
%                               (2D: full [n1 x n2 x 3] for 11/12/22;
%                                3D: cell(6,1) of TT tensors for 11,12,13,22,23,33) — if OPT.stiffness==1
%
%   OPT  Options (fields are optional unless noted):
%          .mass, .stiffness   (0/1) indicate which parts are present (must match H)
%          .rankTol            singular-value / TT-round tolerance (default:
%                               1e-5 if lowRank==1, otherwise Inf)
%          .discardFull        (0/1) remove original full/TT objects after factorization (default 1)
%
%   Outputs
%   -------
%   H  (augmented) factors for univariate integration:
%     If H.dim == 2
%       H.stiffness:
%         .Rmax                 = min(H.weightFun.n)
%         .R(3x1)               selected ranks for (11,12,22)
%         .SVDU (n1 x max(R) x 3), .SVDV (n2 x max(R) x 3)
%         .SVDWeights (max(R) x 3)  with sqrt of singular values
%         (Optional) .weightMat removed when OPT.discardFull==1
%       H.mass:
%         .Rmax                 = min(H.weightFun.n)
%         .R (scalar)           selected rank
%         .SVDU (n1 x R), .SVDV (n2 x R), .SVDWeights (R x 1) with sqrt singular values
%         (Optional) .weightMat removed
%
%     If H.dim == 3
%       H.stiffness:
%         .SVDU{l}{d}          for entry l ∈ {1..6} and direction d ∈ {1,2,3}:
%                                size = [n_d  ,  r_d * r_{d+1}]
%                                (constructed by stacking TT core slices)
%         .R(6 x 3)            directional ranks r_d * r_{d+1} for each l
%         .order               3x3 index map from (i,j) to l:
%                                [1 2 3; 2 4 5; 3 5 6]
%         (Optional) .weightMat removed
%       H.mass:
%         .SVDU{d}             per direction d, size [n_d , r_d * r_{d+1}]
%         .R(1 x 3)            directional ranks r_d * r_{d+1}
%         (Optional) .weightMat removed
%
%   How it works
%   ------------
%   * 2D (matrix SVD):
%       For each matrix A (mass or a stiffness component), compute [U,S,V] = svd(A).
%       Store:
%           SVDU = U(:,1:R),  SVDV = V(:,1:R),  SVDWeights = sqrt(diag(S(1:R,1:R))).
%       This yields A ≈ Σ_{r=1}^R (SVDWeights(r)^2) * (SVDU(:,r) * SVDV(:,r)').
%       (Taking sqrt lets you attach one factor to each 1D side of the quadrature.)
%
%   * 3D (TT to per-direction blocks):
%       Round each TT tensor to OPT.rankTol, then for each direction d and TT ranks r_d, r_{d+1}
%       build a 2D block by stacking slices of the TT core:
%           block_d(:, l + (k-1)*r_{d+1}) = core_d(k, :, l)'   for k=1..r_d, l=1..r_{d+1}.
%       The resulting blocks (SVDU) are exactly the matrices you contract with univariate
%       basis/quad vectors; the directional "ranks" are R_d = r_d * r_{d+1}.
%       For stiffness, six entries are handled separately; H.stiffness.order encodes the
%       (i,j) -> l mapping for reconstructing the symmetric 3x3 matrix.
%
%   Why this format is useful
%   -------------------------
%   Univariate integration with separable weights reduces to 1D contractions.
%   With the factors above, a bilinear form can be applied/assembled as a sum over ranks:
%       Σ_r  (B_1' * SVDU{1}(:,r)) ⊗ (B_2' * SVDU{2}(:,r))   [and ⊗ (B_3' * SVDU{3}(:,r)) in 3D]
%   which avoids ever forming dense multidimensional tensors.
%
%   Notes
%   -----
%   * OPT.mass / OPT.stiffness must reflect which fields exist in H.
%   * If a 3D stiffness component is numerically empty, its factors are set to 0 and its
%     directional ranks to 0.
%   * With OPT.discardFull==1 the original .weightMat fields (full or TT) are removed to
%     save memory; set it to 0 if you still need them.
%
%   Example (2D)
%   ------------
%     % After interpolation:
%     [H, ~, optW] = interpolate_weights_bsplines(G, struct('mass',1,'stiffness',1));
%     opt.rankTol = 1e-8; opt.mass = 1; opt.stiffness = 1;
%     [H, opt] = lowRank_w(H, opt);
%     % Now use H.mass.SVDU/SVDV/SVDWeights and H.stiffness.* per component in 1D quadrature.
%
%   Example (3D)
%   ------------
%     % H.mass.weightMat, H.stiffness.weightMat{k} are TT tensors from interpolation:
%     opt.mass = 1; opt.stiffness = 1; opt.rankTol = 1e-8; opt.discardFull = 1;
%     [H, opt] = lowRank_w(H, opt);
%     % Contract each H.*.SVDU{d} with your univariate basis/quad vectors per direction.
%
%   See also
%   --------
%   INTERPOLATE_WEIGHTS_BSPLINES, INTERPOLATE_WEIGHTS_NURBS,
%   LOWRANKSVD_THB (if used earlier), TT_TENSOR, SVD.

    time = tic;

    if nargin < 2
        opt = struct();
    end

    if ~isfield(opt, 'rankTol') || isempty(opt.rankTol)
        opt.rankTol = 1e-5;
    end
    if ~isfield(opt, 'discardFull') || isempty(opt.discardFull)
        opt.discardFull = 1;
    end
    
    
    if H.dim == 2
        if opt.stiffness == 1
            H.stiffness.Rmax = min(H.weightFun.n);
            H.stiffness.SVDU = zeros(H.weightFun.n(1), H.stiffness.Rmax, 3);
            H.stiffness.SVDV = zeros(H.weightFun.n(2), H.stiffness.Rmax, 3);
            H.stiffness.SVDWeights = zeros(H.stiffness.Rmax, 3);
            H.stiffness.R = ones(3,1)*H.stiffness.Rmax;
            for i = 1:3
                % svd for each entry
                [U,weights1D,V] = svd(H.stiffness.weightMat(:,:,i));
                weights1D = diag(weights1D);

                    % choose low rank according to truncation tolerance
                for j = 1:numel(weights1D)
                    if weights1D(j) < opt.rankTol
                        H.stiffness.R(i) = j-1;
                        break;
                    end
                end

    
                % truncation of the svd matrices
                H.stiffness.SVDU(:,1:H.stiffness.R(i),i) = U(:,1:H.stiffness.R(i));
                H.stiffness.SVDV(:,1:H.stiffness.R(i),i) = V(:,1:H.stiffness.R(i));
                H.stiffness.SVDWeights(1:H.stiffness.R(i),i) = sqrt(weights1D(1:H.stiffness.R(i)));
            end
            H.stiffness.SVDU = H.stiffness.SVDU(:,1:max(H.stiffness.R),:);
            H.stiffness.SVDV = H.stiffness.SVDV(:,1:max(H.stiffness.R),:);
            H.stiffness.SVDWeights = H.stiffness.SVDWeights(1:max(H.stiffness.R),:);
    %         fprintf('Low ranks are %i, %i, %i for tolerance %d\n', H.stiffness.R(1), H.stiffness.R(2), H.stiffness.R(3), opt.rankTol);
            if opt.discardFull == 1
                H.stiffness = rmfield(H.stiffness, 'weightMat');
            end
        end
        if opt.mass == 1
            H.mass.Rmax = min(H.weightFun.n);
            H.mass.R = 0;
            [U, weights1D, V] = svd(H.mass.weightMat);
            weights1D = diag(weights1D);

            H.mass.R = find(weights1D<opt.rankTol,1)-1;

            H.mass.SVDU = U(:,1:H.mass.R);
            H.mass.SVDV = V(:,1:H.mass.R);
            H.mass.SVDWeights = sqrt(weights1D(1:H.mass.R));
            if opt.discardFull == 1
                H.mass = rmfield(H.mass, 'weightMat');
            end
        end
    elseif H.dim == 3
        if opt.stiffness == 1
            H.stiffness.SVDU = cell(6,1);
            H.stiffness.R = zeros(6,3);
            for i = 1:6
                if ~isempty(H.stiffness.weightMat{i})
                    H.stiffness.weightMat{i} = round(H.stiffness.weightMat{i}, opt.rankTol);
                    H.stiffness.SVDU{i} = cell(1,3);
                    for j = 1:3
                        H.stiffness.SVDU{i}{j} = zeros(H.stiffness.weightMat{i}.n(j),H.stiffness.weightMat{i}.r(j)*H.stiffness.weightMat{i}.r(j+1));
                        H.stiffness.R(i,j) = H.stiffness.weightMat{i}.r(j)*H.stiffness.weightMat{i}.r(j+1);
                        for k = 1:H.stiffness.weightMat{i}.r(j)
                            for l = 1:H.stiffness.weightMat{i}.r(j+1)
                                H.stiffness.SVDU{i}{j}(:,l+H.stiffness.weightMat{i}.r(j+1)*(k-1)) = H.stiffness.weightMat{i}{j}(k,:,l)';
                            end
                        end
                    end
                else
                    H.stiffness.SVDU{i} = 0;
                    H.stiffness.R(i,:) = 0;
                end
            end                 
        
            H.stiffness.order = [1,2,3,2,4,5,3,5,6];
            if opt.discardFull == 1
                H.stiffness = rmfield(H.stiffness, 'weightMat');
            end
        end
        if opt.mass == 1
            H.mass.R = zeros(1,3);
            if ~isempty(H.mass.weightMat)
                H.mass.weightMat = round(H.mass.weightMat, opt.rankTol);
                H.mass.SVDU = cell(1,3);
                for j = 1:3
                    H.mass.SVDU{j} = zeros(H.mass.weightMat.n(j),H.mass.weightMat.r(j)*H.mass.weightMat.r(j+1));
                    H.mass.R(j) = H.mass.weightMat.r(j)*H.mass.weightMat.r(j+1);
                    for k = 1:H.mass.weightMat.r(j)
                        for l = 1:H.mass.weightMat.r(j+1)
                            H.mass.SVDU{j}(:,l+H.mass.weightMat.r(j+1)*(k-1)) = H.mass.weightMat{j}(k,:,l)';
                        end
                    end
                end
            else
                H.mass.SVDU = 0;
                H.mass.R(:) = 0;
            end      
            if opt.discardFull == 1
                H.mass = rmfield(H.mass, 'weightMat');
            end
        end
    end
    
    time = toc(time);

end

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

function H = univariate_u_v_bsplines(H, space, n)
    % s = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664];
    % w = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]';
    
    H.mass.M = cell(3,1);
    for dim = 1:3
        H.mass.M{dim} = cell(H.mass.R(dim),1);
        H.mass.M{dim}(:) = {sparse(space.ndof_dir(dim), ...
            space.ndof_dir(dim))};
        for l=1:length(space.knots{dim})-1
            if space.knots{dim}(l) < space.knots{dim}(l+1)
                a = space.knots{dim}(l);
                b = space.knots{dim}(l+1);
                
                % xx = (b-a)/2*s + (a+b)/2;
                
                [xx, w] = gauss_legendre(n(dim), a, b);

                quadValues = evalBSpline(space.knots{dim}, space.degree(dim), xx');
                quadValues2 = evalBSpline(H.weightFun.knots{dim}, H.weightFun.degree(dim), xx');
    
                for i = l-space.degree(dim):l
                    for j = l-space.degree(dim):l
                        for r = 1:H.mass.R(dim)
                            H.mass.M{dim}{r}(i,j) = ...
                                H.mass.M{dim}{r}(i,j) + ...
                                sum(w.*quadValues(i,:)'.*quadValues(j,:)'.*quadValues2'*H.mass.SVDU{dim}(:,r));
                        end
                    end
                end
            end
        end
    end
end

function [x, w] = gauss_legendre(n, a, b)
% GAUSS_LEGENDRE
%
% Cached n-point Gauss-Legendre rule on [a,b].
% The expensive reference rule on [-1,1] is computed only once per n.

    persistent xhat_cache what_cache

    if isempty(xhat_cache)
        xhat_cache = cell(1, 200);
        what_cache = cell(1, 200);
    end

    if n > numel(xhat_cache)
        xhat_cache{n} = [];
        what_cache{n} = [];
    end

    if isempty(xhat_cache{n})
        if n == 1
            xhat = 0;
            what = 2;
        else
            k = 1:n-1;
            beta = k ./ sqrt(4*k.^2 - 1);
            T = diag(beta,1) + diag(beta,-1);

            [V, D] = eig(T);
            xhat = diag(D);
            [xhat, idx] = sort(xhat);
            V = V(:,idx);
            what = 2 * (V(1,:).^2).';
        end

        xhat_cache{n} = xhat;
        what_cache{n} = what;
    else
        xhat = xhat_cache{n};
        what = what_cache{n};
    end

    x = 0.5*((b-a)*xhat + (a+b));
    w = 0.5*(b-a)*what;
end