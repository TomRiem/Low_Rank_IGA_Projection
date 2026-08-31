clear;
clc;
rng("default");

% ========================================================================
% Numerical experiments for the stiffness tensor assembly
%
% This script reproduces the stiffness-tensor experiments from the
% manuscript. For each geometry, solution-space degree, quadrature rule,
% refinement level, and auxiliary-space policy, we compare:
%
%   1. GeoPDEs:       full reference assembly,
%   2. Interpolation: interpolation-based low-rank assembly,
%   3. Projection:    proposed projection-based low-rank assembly.
%
% The low-rank methods are tested for the tolerances
%       1e-5, 1e-7.
%
% Two auxiliary-space policies are considered:
%
%   default:
%       geometry-based space of Eq. (57),
%
%   refined:
%       solution-aligned space of Eq. (58).
%
% In each comparison, interpolation and projection use the same auxiliary
% spline space. Two quadrature choices are used:
%
%   - [5,5,5] as a practical fixed rule,
%   - the sufficient Gaussian rule from Eq. (55).
%
% The script stores assembly time, operator storage, nonzero counts, and
% errors with respect to the GeoPDEs reference. Peak-memory measurements
% reported in the manuscript are obtained separately in isolated processes
% and are therefore not part of this script.
% ========================================================================


%% Common experiment parameters

settings = stiffness_experiment_settings();


%% Thick flag

cfg = flag_configuration(settings);
run_stiffness_experiment(cfg, settings);


%% Rotor blade

cfg = rotor_configuration(settings);
run_stiffness_experiment(cfg, settings);


%% Almost singular geometry

cfg = singularity_configuration(settings);
run_stiffness_experiment(cfg, settings);


fprintf('\nAll stiffness-tensor experiments completed.\n');


function settings = stiffness_experiment_settings()
    settings.tol = [1e-5, 1e-7];

    % Solution-space degrees used in the manuscript.
    settings.degree = [3, 3, 3; 5, 5, 5];

    settings.practical_nquad = [5, 5, 5];

    % Constant diffusion coefficient in the GeoPDEs stiffness form.
    settings.c_diff = @(x, y, z) ones(size(x));

    % Assembly routines used in the comparison.
    settings.interpolate_system = @interpolate_system;
    settings.assemble_interpolation = @assemble_stiffness_lowrank_interpolation;
    settings.assemble_projection = @assemble_stiffness_lowrank_projection;
end



function cfg = flag_configuration(~)
    data = load('3Dflag_GeoPDEs.mat', 'g');

    cfg.name = 'thick flag';
    cfg.geometry = geo_load(data.g.nurbs);
    cfg.result_file = 'flag_stiffness.mat';
    cfg.compare_interior = false;

    % The exact-quadrature experiments use one fewer refinement step because
    % of their larger computational cost.
    p3_practical = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12; 15, 15, 15; 18, 18, 18];

    p5_practical = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12; 15, 15, 15];

    p3_exact = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12; 15, 15, 15];

    p5_exact = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12];

    cfg.n_sub.default = { p3_practical, p3_exact; p5_practical, p5_exact};

    % Both auxiliary-space policies use the same solution refinements.
    cfg.n_sub.refined = cfg.n_sub.default;
end


function cfg = rotor_configuration(~)
    [geometry, ~, ~] = rotor_gen();

    cfg.name = 'rotor blade';
    cfg.geometry = geometry;
    cfg.result_file = 'rotor_stiffness.mat';
    cfg.compare_interior = false;

    p3 = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; 4, 4, 4; 5, 5, 5; 6, 6, 6; 7, 7, 7];

    p5 = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; 4, 4, 4; 5, 5, 5];

    cfg.n_sub.default = {p3, p3; p5, p5};
    cfg.n_sub.refined = cfg.n_sub.default;
end


function cfg = singularity_configuration(~)
    % Trilinear B-spline geometry from the manuscript. The deformation
    % parameter c > -1 keeps the map orientation-preserving. For
    % c = -1 + 1e-5, the Jacobian determinant ranges from 1e-5 to 1.
    c = -1 + 1e-5;

    knots = { [0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]};

    coefs = zeros(4, 2, 2, 2);

    for i = 1:2
        xi = i - 1;

        for j = 1:2
            eta = j - 1;

            for k = 1:2
                zeta = k - 1;

                x = (1 + c * eta * zeta) * xi;
                y = eta;
                z = zeta;

                coefs(:,i,j,k) = [x; y; z; 1];
            end
        end
    end

    cfg.name = 'almost singular geometry';
    cfg.geometry = geo_load(nrbmak(coefs, knots));
    cfg.result_file = 'singularity_stiffness.mat';
    cfg.compare_interior = true;
    cfg.c = c;

    % For this geometry the default and solution-aligned auxiliary spaces
    % are tested over different refinement ranges, as in the manuscript.
    p3_default = [4, 4, 4;  8, 8, 8; 12, 12, 12; 16, 16, 16; 20, 20, 20; 24, 24, 24; 28, 28, 28];

    p5_default = [4, 4, 4;  8, 8, 8; 12, 12, 12; 16, 16, 16; 20, 20, 20];

    p3_refined = [4, 4, 4;  8, 8, 8; 12, 12, 12; 16, 16, 16; 20, 20, 20];

    p5_refined = [4, 4, 4;  8, 8, 8; 12, 12, 12; 16, 16, 16];

    cfg.n_sub.default = { p3_default, p3_default; p5_default, p5_default};

    cfg.n_sub.refined = { p3_refined, p3_refined; p5_refined, p5_refined};
end


function run_stiffness_experiment(cfg, settings)
    fprintf('\n============================================================\n');
    fprintf('Stiffness tensor: %s\n', cfg.name);
    fprintf('============================================================\n');

    tol_n = numel(settings.tol);
    degree_n = size(settings.degree, 1);

    geometry = cfg.geometry;
    p_geo = geometry.nurbs.order(:).' - 1;

    % The default auxiliary space of Eq. (57) depends only on the geometry.
    default_aux_space = make_paper_default_aux_space(geometry);

    experiments = build_auxiliary_experiments(cfg, settings, default_aux_space, p_geo);

    [nquad_all, nsub_all] = collect_unique_configs(experiments, degree_n);

    n_quad_all_max = max(cellfun(@(x) size(x,1), nquad_all));

    % --------------------------------------------------------------------
    % Result structure
    % --------------------------------------------------------------------
    results = struct;

    results.config.geometry = cfg.name;
    results.config.degree = settings.degree;
    results.config.tol = settings.tol;
    results.config.experiments = experiments;
    results.config.nquad_all = nquad_all;
    results.config.nsub_all = nsub_all;
    results.config.p_geo = p_geo;
    results.config.exact_quadrature_formula = 'ceil((2*p_sol + p_aux + 4*p_geo + 1)/2)';

    if isfield(cfg, 'c')
        results.config.c = cfg.c;
    end

    results.gp = initialize_gp_results(degree_n, n_quad_all_max, cfg.compare_interior);
    results.gp.nquad = nquad_all;
    results.gp.nsub = nsub_all;

    for i_exp = 1:numel(experiments)

        tag = experiments(i_exp).tag;
        n_quad_max = max(cellfun(@(x) size(x,1), experiments(i_exp).nquad));

        results.(tag).interpolation = initialize_lr_results(degree_n, n_quad_max, tol_n, cfg.compare_interior);

        results.(tag).projection = initialize_lr_results(degree_n, n_quad_max, tol_n, cfg.compare_interior);

        results.(tag).nquad = experiments(i_exp).nquad;
        results.(tag).nsub = experiments(i_exp).nsub;
        results.(tag).aux_policy = experiments(i_exp).aux_policy;
        results.(tag).aux_degree = cell(degree_n, n_quad_max);
        results.(tag).aux_ndof_dir = cell(degree_n, n_quad_max);
        results.(tag).gp_quad_id = cell(degree_n, 1);

        for i_deg = 1:degree_n
            [is_present, gp_ids] = ismember(experiments(i_exp).nquad{i_deg}, nquad_all{i_deg}, 'rows');

            if ~all(is_present)
                error(['Could not map all quadrature rules of ', 'experiment "%s" to the GeoPDEs rules.'], tag);
            end

            results.(tag).gp_quad_id{i_deg} = gp_ids;
        end
    end

    % --------------------------------------------------------------------
    % Assemble all distinct GeoPDEs references and reuse them for every
    % matching low-rank configuration.
    % --------------------------------------------------------------------
    for i_deg = 1:degree_n

        degree_i = settings.degree(i_deg,:);
        q_all_i = nquad_all{i_deg};

        fprintf('\nSolution degree [%d %d %d]\n', degree_i);

        for i_quad_gp = 1:size(q_all_i,1)

            nquad_i = q_all_i(i_quad_gp,:);
            nsub_all_i = nsub_all{i_deg,i_quad_gp};

            for i_nsub_gp = 1:size(nsub_all_i,1)

                nsub_i = nsub_all_i(i_nsub_gp,:);

                fprintf(['  quadrature [%d %d %d] | ', 'inserted knots [%d %d %d]\n'], nquad_i, nsub_i);

                method_data = build_method_data(degree_i, nsub_i, nquad_i);

                [space, msh] = initialize_laplace(geometry, method_data);

                % --------------------------------------------------------
                % GeoPDEs reference matrix
                % --------------------------------------------------------
                time_gp = tic;
                K_ref = op_gradu_gradv_tp(space, space, msh, settings.c_diff);
                time_gp = toc(time_gp);

                boundary_dofs = [];
                interior_dofs = [];

                if cfg.compare_interior
                    [boundary_dofs, interior_dofs] = determine_boundary_dofs(space, msh);
                end

                [results.gp, ref_data] = store_gp_result(results.gp, i_deg, i_quad_gp, i_nsub_gp, K_ref, time_gp, space, boundary_dofs, interior_dofs);

                % The solution-aligned auxiliary space depends on the
                % current solution discretization, but not on the TT
                % tolerance or on the quadrature rule.
                refined_aux_space = [];

                % --------------------------------------------------------
                % Low-rank methods
                % --------------------------------------------------------
                for i_exp = 1:numel(experiments)

                    tag = experiments(i_exp).tag;
                    q_list = experiments(i_exp).nquad{i_deg};

                    local_quad_ids = find(ismember(q_list, nquad_i, 'rows'));

                    for i_quad = local_quad_ids(:).'

                        nsub_list = experiments(i_exp).nsub{i_deg,i_quad};

                        local_nsub_ids = find(ismember(nsub_list, nsub_i, 'rows'));

                        if isempty(local_nsub_ids)
                            continue;
                        end

                        switch experiments(i_exp).aux_policy

                            case 'paper_default'
                                aux_space = default_aux_space;

                            case 'solution_aligned'
                                if isempty(refined_aux_space)
                                    refined_aux_space = make_solution_aux_space(space);
                                end
                                aux_space = refined_aux_space;

                            otherwise
                                error('Unknown auxiliary-space policy "%s".', experiments(i_exp).aux_policy);
                        end

                        validate_aux_space(aux_space, geometry.rdim, tag);

                        for i_nsub = local_nsub_ids(:).'

                            results.(tag).aux_degree{i_deg,i_quad}(i_nsub,:) = aux_space.degree(:).';

                            results.(tag).aux_ndof_dir{i_deg,i_quad}(i_nsub,:) = aux_space.ndof_dir(:).';

                            for i_tol = 1:tol_n

                                tol_i = settings.tol(i_tol);
                                fprintf('    %s | tolerance %.1e\n', tag, tol_i);

                                % ----------------------------------------
                                % Interpolation-based method
                                % ----------------------------------------
                                [K_interp, time_interp, stats_interp] = run_interpolation_stiffness(geometry, space, method_data, tol_i, aux_space, settings);

                                results.(tag).interpolation = store_lr_result(results.(tag).interpolation, i_deg, i_quad, i_tol, i_nsub, K_interp, time_interp, stats_interp, K_ref, ref_data);

                                clear K_interp

                                % ----------------------------------------
                                % Projection-based method
                                % ----------------------------------------
                                [K_proj, time_proj, stats_proj] = run_projection_stiffness(geometry, space, nquad_i, tol_i, aux_space, settings);

                                results.(tag).projection = store_lr_result(results.(tag).projection, i_deg, i_quad, i_tol, i_nsub, K_proj, time_proj, stats_proj, K_ref, ref_data);

                                clear K_proj

                                save(cfg.result_file, 'results', '-v7.3');
                            end
                        end
                    end
                end

                clear K_ref space msh method_data refined_aux_space
            end
        end
    end

    fprintf('\nCompleted stiffness experiment for %s.\n', cfg.name);
end


function experiments = build_auxiliary_experiments(cfg, settings, default_aux_space, p_geo)
    degree_n = size(settings.degree,1);

    experiments = struct([]);

    experiments(1).tag = 'default';
    experiments(1).aux_policy = 'paper_default';
    experiments(1).nsub = cfg.n_sub.default;

    experiments(2).tag = 'refined';
    experiments(2).aux_policy = 'solution_aligned';
    experiments(2).nsub = cfg.n_sub.refined;

    for i_deg = 1:degree_n

        p_sol = settings.degree(i_deg,:);

        % Eq. (57): fixed geometry-based auxiliary space.
        p_aux_default = default_aux_space.degree(:).';
        nquad_default = sufficient_stiffness_nquad(p_sol, p_aux_default, p_geo);

        % Eq. (58): auxiliary space generated from the solution space.
        p_aux_refined = 3 * p_sol - 1;
        nquad_refined = sufficient_stiffness_nquad(p_sol, p_aux_refined, p_geo);

        experiments(1).nquad{i_deg,1} = [  settings.practical_nquad;     nquad_default];

        experiments(2).nquad{i_deg,1} = [  settings.practical_nquad;     nquad_refined];
    end
end


function aux_space = make_paper_default_aux_space(geometry)
%MAKE_PAPER_DEFAULT_AUX_SPACE  Construct the space from Eq. (57).
%
% In direction d, for geometry degree p and an interior knot of
% multiplicity mu, the space has
%
%   degree:                6*p - 2,
%   endpoint multiplicity: 6*p - 1,
%   interior multiplicity: 5*p + mu - 1.

    dim = geometry.rdim;

    aux_space = struct;
    aux_space.knots = cell(dim,1);
    aux_space.degree = zeros(1,dim);
    aux_space.ndof_dir = zeros(1,dim);

    for d = 1:dim

        U = geometry.nurbs.knots{d}(:).';
        p = geometry.nurbs.order(d) - 1;

        [distinct_knots, ~, knot_ids] = unique(U, 'stable');
        multiplicities = accumarray(knot_ids(:), 1).';

        p_aux = 6 * p - 2;
        U_aux = [];

        for r = 1:numel(distinct_knots)

            if r == 1 || r == numel(distinct_knots)
                m_aux = 6 * p - 1;
            else
                m_aux = 5 * p + multiplicities(r) - 1;
            end

            U_aux = [U_aux, repmat(distinct_knots(r), 1, m_aux)]; %#ok<AGROW>
        end

        aux_space.knots{d} = U_aux;
        aux_space.degree(d) = p_aux;
        aux_space.ndof_dir(d) = numel(U_aux) - p_aux - 1;
    end
end


function nquad = sufficient_stiffness_nquad(p_sol, p_aux, p_geo)
%SUFFICIENT_STIFFNESS_NQUAD  Gaussian rule from Eq. (55).
%
% The maximal polynomial degree in direction d is bounded by
%
%       2*p_sol + p_aux + 4*p_geo,
%
% so that the number of Gauss points is chosen as
%
%       ceil((2*p_sol + p_aux + 4*p_geo + 1)/2).

    p_sol = p_sol(:).';
    p_aux = p_aux(:).';
    p_geo = p_geo(:).';

    if numel(p_sol) ~= numel(p_aux) ||     numel(p_sol) ~= numel(p_geo)
        error('Degree vectors must have the same length.');
    end

    max_degree = 2 * p_sol + p_aux + 4 * p_geo;
    nquad = ceil((max_degree + 1) / 2);
end


function [nquad_all, nsub_all] = collect_unique_configs(experiments, degree_n)
%COLLECT_UNIQUE_CONFIGS
% Collect all distinct (quadrature, refinement) combinations required by
% either auxiliary-space experiment. This avoids repeated GeoPDEs assembly.

    nquad_all = cell(degree_n,1);
    nsub_all = cell(degree_n,1);

    for i_deg = 1:degree_n

        q_all = [];

        for i_exp = 1:numel(experiments)
            q_all = [q_all;         experiments(i_exp).nquad{i_deg}]; %#ok<AGROW>
        end

        nquad_all{i_deg} = unique(q_all, 'rows', 'stable');

        for i_quad_gp = 1:size(nquad_all{i_deg},1)

            q_i = nquad_all{i_deg}(i_quad_gp,:);
            nsub_rows = [];

            for i_exp = 1:numel(experiments)

                q_list = experiments(i_exp).nquad{i_deg};
                local_quad_ids = find(ismember(q_list, q_i, 'rows'));

                for i_quad = local_quad_ids(:).'
                    nsub_rows = [nsub_rows;                 experiments(i_exp).nsub{i_deg,i_quad}]; %#ok<AGROW>
                end
            end

            nsub_all{i_deg,i_quad_gp} = unique(nsub_rows, 'rows', 'stable');
        end
    end
end


function validate_aux_space(aux_space, dim, tag)
    required_fields = {'knots', 'degree', 'ndof_dir'};

    for i_field = 1:numel(required_fields)
        if ~isfield(aux_space, required_fields{i_field})
            error(['Auxiliary space "%s" is missing the ', 'field "%s".'], tag, required_fields{i_field});
        end
    end

    if numel(aux_space.knots) ~= dim ||     numel(aux_space.degree) ~= dim ||     numel(aux_space.ndof_dir) ~= dim
        error(['Auxiliary space "%s" must contain one entry ', 'per parametric direction.'], tag);
    end

    for d = 1:dim

        n_implied = numel(aux_space.knots{d})     - aux_space.degree(d) - 1;

        if n_implied ~= aux_space.ndof_dir(d)
            error(['Auxiliary space "%s" is inconsistent in ', 'direction %d.'], tag, d);
        end
    end
end


function method_data = build_method_data(degree, n_sub, n_quad)
    method_data = struct;
    method_data.degree = degree;
    method_data.regularity = degree - 1;
    method_data.nsub = n_sub;
    method_data.nquad = n_quad;
    method_data.basis_functions = 'B-Splines';
    method_data.dirichlet = 0;
end


function [boundary_dofs, interior_dofs] = determine_boundary_dofs(space, msh)
%DETERMINE_BOUNDARY_DOFS
% Determine all boundary degrees of freedom. This is used only for the
% almost singular geometry, for which the manuscript reports the error on
% the interior degrees of freedom.

    drchlt_sides = [1, 2, 3, 4, 5, 6];
    h = @(x, y, z, ind) 0 .* x; %#ok<INUSD>

    [~, boundary_dofs] = sp_drchlt_l2_proj(space, msh, h, drchlt_sides);

    boundary_dofs = unique(boundary_dofs(:)).';
    interior_dofs = setdiff(1:space.ndof, boundary_dofs);
end


function results = initialize_gp_results(degree_n, n_quad_max, compare_interior)
    results = struct;
    results.memory = cell(degree_n, n_quad_max);
    results.nnz = cell(degree_n, n_quad_max);
    results.time = cell(degree_n, n_quad_max);
    results.norm = cell(degree_n, n_quad_max);
    results.norm_est = cell(degree_n, n_quad_max);
    results.ndof = cell(degree_n, n_quad_max);

    if compare_interior
        % Descriptive fields used for the almost singular geometry.
        results.norm_with_boundary = cell(degree_n, n_quad_max);
        results.norm_est_with_boundary = cell(degree_n, n_quad_max);
        results.norm_without_boundary = cell(degree_n, n_quad_max);
        results.norm_est_without_boundary = cell(degree_n, n_quad_max);

        % Backward-compatible aliases.
        results.norm_inner = cell(degree_n, n_quad_max);
        results.norm_est_inner = cell(degree_n, n_quad_max);

        results.n_boundary_dofs = cell(degree_n, n_quad_max);
        results.n_interior_dofs = cell(degree_n, n_quad_max);
    end
end


function results = initialize_lr_results(degree_n, n_quad_max, tol_n, compare_interior)
    sz = [degree_n, n_quad_max, tol_n];

    results = struct;
    results.memory = cell(sz);
    results.time = cell(sz);
    results.abs_err = cell(sz);
    results.rel_err = cell(sz);
    results.abs_err_est = cell(sz);
    results.rel_err_est = cell(sz);

    if compare_interior
        % Descriptive fields for errors with and without boundary DoFs.
        results.abs_err_with_boundary = cell(sz);
        results.rel_err_with_boundary = cell(sz);
        results.abs_err_est_with_boundary = cell(sz);
        results.rel_err_est_with_boundary = cell(sz);

        results.abs_err_without_boundary = cell(sz);
        results.rel_err_without_boundary = cell(sz);
        results.abs_err_est_without_boundary = cell(sz);
        results.rel_err_est_without_boundary = cell(sz);

        % Backward-compatible aliases for the interior error.
        results.abs_err_inner = cell(sz);
        results.rel_err_inner = cell(sz);
        results.abs_err_est_inner = cell(sz);
        results.rel_err_est_inner = cell(sz);
    end
end


function [results, ref] = store_gp_result(results, i_deg, i_quad, i_nsub, K_ref, assembly_time, space, boundary_dofs, interior_dofs)
    results.time{i_deg,i_quad}(i_nsub) = assembly_time;
    results.memory{i_deg,i_quad}(i_nsub) = RecursiveSize(K_ref);
    results.nnz{i_deg,i_quad}(i_nsub) = nnz(K_ref);
    results.norm{i_deg,i_quad}(i_nsub) = norm(K_ref, 'fro');
    results.norm_est{i_deg,i_quad}(i_nsub) = normest(K_ref, 1e-10);
    results.ndof{i_deg,i_quad}(i_nsub) = space.ndof;

    ref = struct;
    ref.norm = results.norm{i_deg,i_quad}(i_nsub);
    ref.norm_est = results.norm_est{i_deg,i_quad}(i_nsub);
    ref.interior_dofs = interior_dofs;

    if ~isempty(interior_dofs)

        K_inner = K_ref(interior_dofs, interior_dofs);

        norm_inner = norm(K_inner, 'fro');
        norm_est_inner = normest(K_inner, 1e-10);

        results.norm_with_boundary{i_deg,i_quad}(i_nsub) = ref.norm;
        results.norm_est_with_boundary{i_deg,i_quad}(i_nsub) = ref.norm_est;
        results.norm_without_boundary{i_deg,i_quad}(i_nsub) = norm_inner;
        results.norm_est_without_boundary{i_deg,i_quad}(i_nsub) = norm_est_inner;

        results.norm_inner{i_deg,i_quad}(i_nsub) = norm_inner;
        results.norm_est_inner{i_deg,i_quad}(i_nsub) = norm_est_inner;

        results.n_boundary_dofs{i_deg,i_quad}(i_nsub) = numel(boundary_dofs);
        results.n_interior_dofs{i_deg,i_quad}(i_nsub) = numel(interior_dofs);

        ref.norm_inner = norm_inner;
        ref.norm_est_inner = norm_est_inner;
    end
end


function results = store_lr_result(results, i_deg, i_quad, i_tol, i_nsub, K_tt, time_total, memory_stats, K_ref, ref)
    idx = {i_deg, i_quad, i_tol};

    results.memory{idx{:}}(i_nsub) = RecursiveSize(K_tt);
    results.time{idx{:}}(i_nsub) = time_total;

    % Store all method-specific memory statistics without hard-coding the
    % individual fields in the experiment driver.
    stat_names = fieldnames(memory_stats);

    for i_stat = 1:numel(stat_names)

        stat_name = stat_names{i_stat};

        if ~isfield(results, stat_name)
            results.(stat_name) = cell(size(results.time));
        end

        results.(stat_name){idx{:}}(i_nsub) = memory_stats.(stat_name);
    end

    % Conversion to a sparse matrix is used only for the error evaluation.
    K_full = sparse(full(K_tt));
    D = K_full - K_ref;

    abs_err = norm(D, 'fro');
    results.abs_err{idx{:}}(i_nsub) = abs_err;
    results.rel_err{idx{:}}(i_nsub) = abs_err / ref.norm;

    abs_err_est = normest(D, 1e-10);
    results.abs_err_est{idx{:}}(i_nsub) = abs_err_est;
    results.rel_err_est{idx{:}}(i_nsub) = abs_err_est / ref.norm_est;

    % For the almost singular geometry, the manuscript additionally
    % compares the matrices after removing all boundary rows and columns.
    if ~isempty(ref.interior_dofs)

        K_inner = K_full(ref.interior_dofs, ref.interior_dofs);
        K_ref_inner = K_ref(ref.interior_dofs, ref.interior_dofs);
        D_inner = K_inner - K_ref_inner;

        abs_err_inner = norm(D_inner, 'fro');
        rel_err_inner = abs_err_inner / ref.norm_inner;

        abs_err_est_inner = normest(D_inner, 1e-10);
        rel_err_est_inner = abs_err_est_inner / ref.norm_est_inner;

        % Explicit with-boundary fields.
        results.abs_err_with_boundary{idx{:}}(i_nsub) = results.abs_err{idx{:}}(i_nsub);
        results.rel_err_with_boundary{idx{:}}(i_nsub) = results.rel_err{idx{:}}(i_nsub);
        results.abs_err_est_with_boundary{idx{:}}(i_nsub) = results.abs_err_est{idx{:}}(i_nsub);
        results.rel_err_est_with_boundary{idx{:}}(i_nsub) = results.rel_err_est{idx{:}}(i_nsub);

        % Explicit without-boundary fields.
        results.abs_err_without_boundary{idx{:}}(i_nsub) = abs_err_inner;
        results.rel_err_without_boundary{idx{:}}(i_nsub) = rel_err_inner;
        results.abs_err_est_without_boundary{idx{:}}(i_nsub) = abs_err_est_inner;
        results.rel_err_est_without_boundary{idx{:}}(i_nsub) = rel_err_est_inner;

        % Backward-compatible interior aliases.
        results.abs_err_inner{idx{:}}(i_nsub) = abs_err_inner;
        results.rel_err_inner{idx{:}}(i_nsub) = rel_err_inner;
        results.abs_err_est_inner{idx{:}}(i_nsub) = abs_err_est_inner;
        results.rel_err_est_inner{idx{:}}(i_nsub) = rel_err_est_inner;
    end
end


function [K_tt, time_total, memory_stats] = run_interpolation_stiffness(geometry, space, method_data, tol, aux_space, settings)
    time_total = tic;

    low_rank_data = struct;
    low_rank_data.mass = 0;
    low_rank_data.stiffness = 1;
    low_rank_data.TT_interpolation = 1;
    low_rank_data.boundary_conditions = 'Dirichlet';
    low_rank_data.geometry_format = 'B-Splines';
    low_rank_data.rankTol = tol;

    % For a fair comparison, interpolation and projection use exactly the
    % same auxiliary spline space.
    low_rank_data.aux_space = aux_space;

    [H, stats_interp] = settings.interpolate_system(geometry, low_rank_data);

    [K_tt, stats_assembly] = settings.assemble_interpolation(H, space, low_rank_data, method_data);

    time_total = toc(time_total);

    memory_stats = struct;

    memory_stats.interpolation_memory = stats_interp.A_memory + stats_interp.b_memory;
    memory_stats.interpolation_nnz = stats_interp.A_nnz + stats_interp.b_nnz;

    memory_stats.cores_memory = stats_assembly.cores_memory;
    memory_stats.cores_nnz = stats_assembly.cores_nnz;

    memory_stats.assembly_memory = memory_stats.interpolation_memory + memory_stats.cores_memory;
    memory_stats.assembly_nnz = memory_stats.interpolation_nnz + memory_stats.cores_nnz;

    memory_stats.K_memory = RecursiveSize(K_tt);
    memory_stats.K_nnz = count_tt_nnz(K_tt);

    memory_stats.total_memory = memory_stats.assembly_memory + memory_stats.K_memory;
    memory_stats.total_nnz = memory_stats.assembly_nnz + memory_stats.K_nnz;
end


function [K_tt, time_total, memory_stats] = run_projection_stiffness(geometry, space, n_quad, tol, aux_space, settings)
    time_total = tic;

    % proj_ref = 2 selects the explicitly supplied auxiliary space.
    proj_ref = 2;

    [K_tt, stats] = settings.assemble_projection(tol, geometry, space, n_quad, proj_ref, aux_space);

    time_total = toc(time_total);

    memory_stats = struct;

    memory_stats.transfer_memory = stats.transfer_total_memory;
    memory_stats.transfer_nnz = stats.transfer_total_nnz;

    memory_stats.projection_memory = stats.projection_total_memory;
    memory_stats.projection_nnz = stats.projection_total_nnz;

    % The geometry coefficient tensors and the projected reciprocal are
    % stored after the projection and are therefore reported separately.
    if isfield(stats, 'post_projection_total_memory')
        memory_stats.post_projection_memory = stats.post_projection_total_memory;
        memory_stats.post_projection_nnz = stats.post_projection_total_nnz;
    else
        memory_stats.post_projection_memory = stats.geometry_coeff_memory + stats.denom_memory;
        memory_stats.post_projection_nnz = stats.geometry_coeff_nnz + stats.denom_nnz;
    end

    memory_stats.assembly_memory = memory_stats.transfer_memory + memory_stats.projection_memory + memory_stats.post_projection_memory;

    memory_stats.assembly_nnz = memory_stats.transfer_nnz + memory_stats.projection_nnz + memory_stats.post_projection_nnz;

    memory_stats.K_memory = stats.K_memory;
    memory_stats.K_nnz = stats.K_nnz;

    memory_stats.total_memory = memory_stats.assembly_memory + memory_stats.K_memory;
    memory_stats.total_nnz = memory_stats.assembly_nnz + memory_stats.K_nnz;
end


function n = count_tt_nnz(T)
    try
        n = nnz(T);
    catch
        cores = core2cell(T);
        n = sum(cellfun(@nnz, cores));
    end
end