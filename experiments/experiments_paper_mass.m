clear;
clc;
rng("default");

% ========================================================================
% Numerical experiments for the mass tensor assembly
%
% This script reproduces the mass-tensor experiments from the manuscript.
% For each geometry, solution-space degree, quadrature rule, and refinement
% level, we compare:
%
%   1. GeoPDEs:       full reference assembly,
%   2. Interpolation: interpolation-based low-rank assembly,
%   3. Projection:    proposed projection-based low-rank assembly.
%
% The low-rank methods are tested for the tolerances
%       1e-3, 1e-5, 1e-7.
%
% Two quadrature choices are used in the final mass assembly:
%   - [5,5,5] as a practical fixed rule,
%   - the exact Gaussian rule from Eq. (54) of the manuscript.
%
% The script stores assembly time, operator storage, nonzero counts, and
% errors with respect to the GeoPDEs reference. Peak-memory measurements
% reported in the manuscript are obtained separately in isolated processes
% and are therefore not part of this script.
% ========================================================================


%% Common experiment parameters

settings = mass_experiment_settings();


%% Twisted beam

cfg = twisted_beam_configuration(settings);
run_mass_experiment(cfg, settings);


%% Thick flag

cfg = flag_configuration(settings);
run_mass_experiment(cfg, settings);


%% Rotor blade

cfg = rotor_configuration(settings);
run_mass_experiment(cfg, settings);


fprintf('\nAll mass-tensor experiments completed.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMON EXPERIMENT SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function settings = mass_experiment_settings()
    settings.tol = [1e-3, 1e-5, 1e-7];

    % Solution-space degrees used in the manuscript.
    settings.degree = [3, 3, 3; 5, 5, 5];

    % Final quadrature rules.
    %
    % For geometry degree p = 2, Eq. (54) gives
    %   ceil((2*phat + 3*p)/2)
    % quadrature points per direction. Hence:
    %
    %   phat = 3  -> 6 points,
    %   phat = 5  -> 8 points.
    %
    % The first row for each degree is the fixed practical rule [5,5,5].
    settings.n_quadrature = cell(2,1);
    settings.n_quadrature{1} = [5, 5, 5; 6, 6, 6];

    settings.n_quadrature{2} = [5, 5, 5; 8, 8, 8];

    % Constant coefficient in the GeoPDEs mass form.
    settings.c_diff = @(x, y, z) ones(size(x));

    % Assembly routines used in the comparison.
    settings.interpolate_system = @interpolate_system;

    % Resolve the assembly routine names once. 
    settings.assemble_interpolation = @assemble_mass_lowrank_interpolation;

    settings.assemble_projection = @assemble_mass_lowrank_projection;
end


function cfg = twisted_beam_configuration(settings)
    data = load('geo_twisted_pipe.mat', 'geo');

    cfg.name = 'twisted beam';
    cfg.geometry = geo_load(data.geo);
    cfg.file_prefix = 'beam_mass';

    % Number of uniformly inserted knots in each parametric direction.
    % The same refinement sequence is used for both quadrature rules.
    ref_p3 = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; 4, 4, 4; 5, 5, 5; 6, 6, 6; 7, 7, 7];

    ref_p5 = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; 4, 4, 4; 5, 5, 5];

    cfg.n_sub = repeat_refinements_for_quadrature({ref_p3, ref_p5}, settings.n_quadrature);
end


function cfg = flag_configuration(settings)
    data = load('3Dflag_GeoPDEs.mat', 'g');

    cfg.name = 'thick flag';
    cfg.geometry = geo_load(data.g.nurbs);
    cfg.file_prefix = 'flag_mass';

    % The exact quadrature experiments use one fewer refinement step in
    % the manuscript because of their larger computational cost.
    cfg.n_sub = cell(2,2);

    cfg.n_sub{1,1} = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12; 15, 15, 15; 18, 18, 18];

    cfg.n_sub{2,1} = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12; 15, 15, 15];

    cfg.n_sub{1,2} = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12; 15, 15, 15];

    cfg.n_sub{2,2} = [0, 0, 0;  3, 3, 3;  6, 6, 6;  9, 9, 9; 12, 12, 12];
end


function cfg = rotor_configuration(settings)
    [geometry, ~, ~] = rotor_gen();

    cfg.name = 'rotor blade';
    cfg.geometry = geometry;
    cfg.file_prefix = 'rotor_mass';

    ref_p3 = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; 4, 4, 4; 5, 5, 5; 6, 6, 6; 7, 7, 7];

    ref_p5 = [0, 0, 0; 1, 1, 1; 2, 2, 2; 3, 3, 3; 4, 4, 4; 5, 5, 5];

    cfg.n_sub = repeat_refinements_for_quadrature({ref_p3, ref_p5}, settings.n_quadrature);
end


function n_sub = repeat_refinements_for_quadrature(refinement_by_degree, n_quadrature)
    degree_n = numel(refinement_by_degree);
    n_quad_max = max(cellfun(@(q) size(q,1), n_quadrature));

    n_sub = cell(degree_n, n_quad_max);

    for i_deg = 1:degree_n
        n_quad = size(n_quadrature{i_deg}, 1);

        for i_quad = 1:n_quad
            n_sub{i_deg, i_quad} = refinement_by_degree{i_deg};
        end
    end
end


function run_mass_experiment(cfg, settings)
    fprintf('\n============================================================\n');
    fprintf('Mass tensor: %s\n', cfg.name);
    fprintf('============================================================\n');

    tol_n = numel(settings.tol);
    degree_n = size(settings.degree, 1);
    n_quad_max = max(cellfun(@(q) size(q,1), settings.n_quadrature));

    results_gp = initialize_gp_results(degree_n, n_quad_max);
    results_interpolation = initialize_interpolation_results(degree_n, n_quad_max, tol_n);
    results_projection = initialize_projection_results(degree_n, n_quad_max, tol_n);

    for i_deg = 1:degree_n

        degree = settings.degree(i_deg,:);
        quadrature_rules = settings.n_quadrature{i_deg};

        for i_quad = 1:size(quadrature_rules,1)

            n_quad = quadrature_rules(i_quad,:);
            refinements = cfg.n_sub{i_deg, i_quad};

            for i_ref = 1:size(refinements,1)

                n_sub = refinements(i_ref,:);

                fprintf(['\n%s | degree [%d %d %d] | quadrature [%d %d %d] ', '| inserted knots [%d %d %d]\n'], cfg.name, degree, n_quad, n_sub);

                method_data = build_method_data(degree, n_sub, n_quad);
                [space, msh] = initialize_laplace(cfg.geometry, method_data);

                % --------------------------------------------------------
                % GeoPDEs reference matrix
                % --------------------------------------------------------
                time_gp = tic;
                M_ref = op_u_v_tp(space, space, msh, settings.c_diff);
                time_gp = toc(time_gp);

                results_gp = store_gp_result(results_gp, i_deg, i_quad, M_ref, space.ndof, time_gp);

                save_results(cfg.file_prefix, 'gp', results_gp);

                % The reference norms are reused for all low-rank
                % tolerances at the current discretization.
                ref_norm = results_gp.norm{i_deg,i_quad}(end);
                ref_norm_est = results_gp.norm_est{i_deg,i_quad}(end);

                % --------------------------------------------------------
                % Low-rank methods
                % --------------------------------------------------------
                for i_tol = 1:tol_n

                    tol = settings.tol(i_tol);

                    fprintf('  tolerance %.1e\n', tol);

                    % ----------------------------------------------------
                    % Interpolation-based method
                    % ----------------------------------------------------
                    low_rank_data = build_low_rank_data(tol);

                    [H, stats_interp] = settings.interpolate_system(cfg.geometry, low_rank_data);

                    [M_interp, stats_interp_assembly] = settings.assemble_interpolation(H, space, low_rank_data, method_data);

                    results_interpolation = store_interpolation_result(results_interpolation, i_deg, i_quad, i_tol, M_interp, M_ref, ref_norm, ref_norm_est, stats_interp, stats_interp_assembly);

                    save_results(cfg.file_prefix, 'interpolation', results_interpolation);

                    % ----------------------------------------------------
                    % Projection-based method proposed in the manuscript
                    % ----------------------------------------------------
                    [M_proj, stats_proj] = settings.assemble_projection(tol, cfg.geometry, space, n_quad);

                    results_projection = store_projection_result(results_projection, i_deg, i_quad, i_tol, M_proj, M_ref, ref_norm, ref_norm_est, stats_proj);

                    save_results(cfg.file_prefix, 'projection', results_projection);
                end
            end
        end
    end

    fprintf('\nCompleted mass experiment for %s.\n', cfg.name);
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


function low_rank_data = build_low_rank_data(tol)
    low_rank_data = struct;
    low_rank_data.mass = 1;
    low_rank_data.stiffness = 0;
    low_rank_data.TT_interpolation = 1;
    low_rank_data.boundary_conditions = 'Dirichlet';
    low_rank_data.geometry_format = 'B-Splines';
    low_rank_data.rankTol = tol;
end


function results = initialize_gp_results(degree_n, n_quad_max)
    results = struct;
    results.M_memory = cell(degree_n, n_quad_max);
    results.M_nnz = cell(degree_n, n_quad_max);
    results.time = cell(degree_n, n_quad_max);
    results.norm = cell(degree_n, n_quad_max);
    results.norm_est = cell(degree_n, n_quad_max);
    results.ndof = cell(degree_n, n_quad_max);
end


function results = initialize_interpolation_results(degree_n, n_quad_max, tol_n)
    sz = [degree_n, n_quad_max, tol_n];

    results = struct;
    results.M_memory = cell(sz);
    results.time = cell(sz);
    results.rel_err = cell(sz);
    results.abs_err = cell(sz);
    results.rel_err_est = cell(sz);
    results.abs_err_est = cell(sz);
    results.M_nnz = cell(sz);

    results.time_interpol = cell(sz);
    results.A_memory = cell(sz);
    results.A_nnz = cell(sz);
    results.b_memory = cell(sz);
    results.b_nnz = cell(sz);
    results.cores_memory = cell(sz);
    results.cores_nnz = cell(sz);
    results.interpolation_memory = cell(sz);
    results.interpolation_nnz = cell(sz);
end


function results = initialize_projection_results(degree_n, n_quad_max, tol_n)
    sz = [degree_n, n_quad_max, tol_n];

    results = struct;
    results.M_memory = cell(sz);
    results.time = cell(sz);
    results.rel_err = cell(sz);
    results.abs_err = cell(sz);
    results.rel_err_est = cell(sz);
    results.abs_err_est = cell(sz);
    results.M_nnz = cell(sz);

    % Stored components of the projection-based assembly.
    results.C_memory = cell(sz);
    results.C_nnz = cell(sz);
    results.B_memory = cell(sz);
    results.B_nnz = cell(sz);
    results.Proj_memory = cell(sz);
    results.Proj_nnz = cell(sz);
    results.assembly_memory = cell(sz);
    results.assembly_nnz = cell(sz);
end


function results = store_gp_result(results, i_deg, i_quad, M_ref, ndof, assembly_time)
    results.time{i_deg,i_quad}(end+1) = assembly_time;
    results.M_memory{i_deg,i_quad}(end+1) = RecursiveSize(M_ref);
    results.M_nnz{i_deg,i_quad}(end+1) = nnz(M_ref);
    results.norm{i_deg,i_quad}(end+1) = norm(M_ref, 'fro');
    results.norm_est{i_deg,i_quad}(end+1) = normest(M_ref, 1e-10);
    results.ndof{i_deg,i_quad}(end+1) = ndof;
end


function results = store_interpolation_result(results, i_deg, i_quad, i_tol, M_tt, M_ref, ref_norm, ref_norm_est, stats_interp, stats_assembly)
    idx = {i_deg, i_quad, i_tol};

    results.M_memory{idx{:}}(end+1) = RecursiveSize(M_tt);
    results.M_nnz{idx{:}}(end+1) = nnz(M_tt);

    % The interpolation time contains both the construction of the
    % interpolation coefficients and the subsequent low-rank assembly.
    results.time{idx{:}}(end+1) = stats_interp.time + stats_assembly.time;

    results.time_interpol{idx{:}}(end+1) = stats_interp.time_interpol;

    results.A_memory{idx{:}}(end+1) = stats_interp.A_memory;
    results.A_nnz{idx{:}}(end+1) = stats_interp.A_nnz;
    results.b_memory{idx{:}}(end+1) = stats_interp.b_memory;
    results.b_nnz{idx{:}}(end+1) = stats_interp.b_nnz;

    results.cores_memory{idx{:}}(end+1) = stats_assembly.cores_memory;
    results.cores_nnz{idx{:}}(end+1) = stats_assembly.cores_nnz;

    results.interpolation_memory{idx{:}}(end+1) = stats_interp.A_memory + stats_interp.b_memory;

    results.interpolation_nnz{idx{:}}(end+1) = stats_interp.A_nnz + stats_interp.b_nnz;

    [abs_err, rel_err, abs_err_est, rel_err_est] = mass_errors(M_tt, M_ref, ref_norm, ref_norm_est);

    results.abs_err{idx{:}}(end+1) = abs_err;
    results.rel_err{idx{:}}(end+1) = rel_err;
    results.abs_err_est{idx{:}}(end+1) = abs_err_est;
    results.rel_err_est{idx{:}}(end+1) = rel_err_est;
end


function results = store_projection_result(results, i_deg, i_quad, i_tol, M_tt, M_ref, ref_norm, ref_norm_est, stats)
    idx = {i_deg, i_quad, i_tol};

    results.M_memory{idx{:}}(end+1) = RecursiveSize(M_tt);
    results.M_nnz{idx{:}}(end+1) = nnz(M_tt);
    results.time{idx{:}}(end+1) = stats.time;

    results.C_memory{idx{:}}(end+1) = stats.C_memory;
    results.C_nnz{idx{:}}(end+1) = stats.C_nnz;

    results.B_memory{idx{:}}(end+1) = stats.B_memory;
    results.B_nnz{idx{:}}(end+1) = stats.B_nnz;

    results.Proj_memory{idx{:}}(end+1) = stats.Proj_memory;
    results.Proj_nnz{idx{:}}(end+1) = stats.Proj_nnz;

    results.assembly_memory{idx{:}}(end+1) = stats.C_memory + stats.B_memory + stats.Proj_memory;

    results.assembly_nnz{idx{:}}(end+1) = stats.C_nnz + stats.B_nnz + stats.Proj_nnz;

    [abs_err, rel_err, abs_err_est, rel_err_est] = mass_errors(M_tt, M_ref, ref_norm, ref_norm_est);

    results.abs_err{idx{:}}(end+1) = abs_err;
    results.rel_err{idx{:}}(end+1) = rel_err;
    results.abs_err_est{idx{:}}(end+1) = abs_err_est;
    results.rel_err_est{idx{:}}(end+1) = rel_err_est;
end


function [abs_err, rel_err, abs_err_est, rel_err_est] = mass_errors(M_tt, M_ref, ref_norm, ref_norm_est)
    % Conversion to a sparse matrix is used only for the error evaluation,
    % as described in the numerical-experiment section of the manuscript.
    M_full = sparse(full(M_tt));
    D = M_full - M_ref;

    abs_err = norm(D, 'fro');
    rel_err = abs_err / ref_norm;

    % Retained as an auxiliary diagnostic for compatibility with the
    % existing result files. The manuscript reports the Frobenius error.
    abs_err_est = normest(D, 1e-10);
    rel_err_est = abs_err_est / ref_norm_est;
end


function save_results(prefix, method, value)
    % Preserve both the original file names and the variable names stored
    % in the MAT files so that existing plotting scripts remain usable.
    switch method
        case 'gp'
            file_name = [prefix, '_gp.mat'];
            variable_name = 'results_gp';

        case 'interpolation'
            file_name = [prefix, '_interpolation.mat'];
            variable_name = 'results_interpolation';

        case 'projection'
            file_name = [prefix, '_projection.mat'];
            variable_name = 'results_projection';

        otherwise
            error('Unknown method "%s".', method);
    end

    S = struct;
    S.(variable_name) = value;
    save(file_name, '-struct', 'S', variable_name);
end

