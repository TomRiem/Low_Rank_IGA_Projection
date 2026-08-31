function varargout = measure_stiffness_peak_memory(action, varargin)
%MEASURE_STIFFNESS_PEAK_MEMORY  Peak-memory workflow for the stiffness tensor.
%
% This file is the stiffness counterpart of measure_mass_peak_memory.m and
% the memory counterpart of experiments_paper_stiffness.m. That script
% measures assembly time, operator storage and accuracy; this one measures
% the transient memory required during the assembly itself.
%
% For each geometry, auxiliary-space policy, quadrature rule, solution-space
% degree, refinement level and rank tolerance the peak memory of
%
%   1. GeoPDEs:       full reference assembly,
%   2. Interpolation: interpolation-based low-rank assembly,
%   3. Projection:    proposed projection-based low-rank assembly,
%
% is measured in a separate MATLAB process. Peak memory cannot be measured
% reliably inside one long session, because the resident set of a process
% does not shrink when MATLAB frees an array; every measurement would then
% depend on the history of the session. Each case therefore runs in a fresh
% process whose resident-set high-water mark is read from
% /proc/self/status. The measurement requires Linux.
%
% Per-case results are stored as
%
%   peak_out/<experiment>_stiff/peak_<method>_<policy>_<quad>_<ideg>_<insub>_<itol>.mat
%
% with method in {gp, interpolation, projection}, policy in {default,
% refined} and quad in {practical, exact}.
%
% -------------------------------------------------------------------------
% AUXILIARY-SPACE POLICIES
% -------------------------------------------------------------------------
%
%   'default'  geometry-based space of degree 6*p_geo - 2, built once per
%              geometry by make_paper_default_aux_space,
%   'refined'  solution-aligned space of degree 3*p_sol - 1, built per
%              solution space by make_solution_aux_space.
%
% Interpolation and projection always use the same auxiliary space, so that
% the comparison isolates the assembly strategy.
%
% -------------------------------------------------------------------------
% QUADRATURE
% -------------------------------------------------------------------------
%
%   'practical'  the fixed rule [5 5 5],
%   'exact'      the sufficient Gaussian rule
%                   ceil((2*p_sol + p_aux + 4*p_geo + 1)/2).
%
% Both rules are kept by default, unlike the mass experiments: the exact
% rule is large and drives the stiffness peak. As in
% experiments_paper_stiffness.m, the GeoPDEs reference is assembled with
% the SAME rule as the low-rank methods it is compared against, so the
% exact rule depends on the auxiliary-space policy and GeoPDEs is measured
% once per policy.
%
% -------------------------------------------------------------------------
% WORKING WITH EXISTING MEASUREMENTS
% -------------------------------------------------------------------------
%
% Collect one slice, that is one policy and one quadrature kind:
%
%   results = measure_stiffness_peak_memory('collect', 'beam', 'refined', 'exact');
%
% The returned structure contains
%
%   results.gp
%   results.interpolation
%   results.projection
%   results.config
%
% Collect every slice for which measurements exist:
%
%   results = measure_stiffness_peak_memory('collect_all');
%
% Plot an already measured slice:
%
%   measure_stiffness_peak_memory('plot', 'beam', 'refined', 'exact');
%
% Plot the stored operator memory instead of the peak memory:
%
%   measure_stiffness_peak_memory('plot', 'beam', 'refined', 'exact', 'metric', 'K_memory');
%
% -------------------------------------------------------------------------
% RUNNING THE SWEEP
% -------------------------------------------------------------------------
%
% Prepare the isolated jobs from an interactive session that has GeoPDEs,
% the TT toolbox and this repository on the path:
%
%   measure_stiffness_peak_memory('prepare', 'beam');
%
% This saves a snapshot of the search path in session_path.mat, so that
% every child process reproduces the present environment, and writes
% run_peak_jobs_beam_stiff.sh. Run that script in a terminal. The
% interactive session may be closed first, which frees the licence seat for
% the children. Completed cases are skipped, so the script may be
% interrupted and restarted.
%
% Subsets are selected with
%
%   measure_stiffness_peak_memory('prepare', 'beam', 'policies', {'refined'});
%   measure_stiffness_peak_memory('prepare', 'beam', 'quads', {'exact'});
%   measure_stiffness_peak_memory('prepare', 'beam', 'methods', {'gp','projection'});
%
% -------------------------------------------------------------------------
% WHAT IS MEASURED
% -------------------------------------------------------------------------
%
% For every isolated case the sequence is
%
%   1. build geometry, mesh, solution space and auxiliary space,
%   2. record VmRSS as the baseline,
%   3. reset VmHWM through /proc/self/clear_refs and verify the reset,
%   4. execute exactly one assembly,
%   5. immediately read VmHWM,
%   6. only afterwards evaluate RecursiveSize and nnz,
%
% so that
%
%   peak_memory = VmHWM - VmRSS_baseline.
%
% The discretisation and the auxiliary space of step 1 are common to the
% methods being compared and are deliberately excluded. Nothing that
% allocates is inserted between steps 3 and 5.
%
% Supported actions
% -----------------
%   'collect'      collect one policy and quadrature slice of one geometry
%   'collect_all'  collect every slice with existing case files
%   'plot'         collect and plot one slice
%   'prepare'      generate the isolated shell jobs
%   'config'       return the experiment configuration
%   'worker'       internal mode used by the generated child jobs
%   'help'         print a short usage summary
%
% Supported experiments
% ---------------------
%   'beam', 'flag', 'rotor', 'singularity'.
%
% The manuscript reports the flag, the rotor and the almost singular
% geometry. The twisted beam is available as an additional case; its
% interpolation system is singular, so interpolation is omitted for it.
%
% The refinement grids depend on the policy and on the quadrature kind and
% are taken from experiments_paper_stiffness.m.

    if nargin == 0
        action = 'help';
    end

    action = lower(char(action));

    switch action

        case 'help'
            print_usage();

        case 'config'
            require_nargs(varargin, 1, 'config');

            varargout{1} = stiffness_peak_config(varargin{1});

        case 'collect'
            require_nargs(varargin, 3, 'collect');

            varargout{1} = collect_experiment(varargin{1}, varargin{2}, varargin{3});

        case 'collect_all'
            varargout{1} = collect_all_experiments();

        case 'plot'
            require_nargs(varargin, 3, 'plot');

            experiment = varargin{1};
            policy = varargin{2};
            quad_kind = varargin{3};
            plot_options = varargin(4:end);

            cfg = stiffness_peak_config(experiment);

            [results_gp, results_interpolation, results_projection] = collect_peaks(experiment, policy, quad_kind);

            out = plot_peak_scaling(results_gp, results_interpolation, results_projection, cfg, ...
                                    policy, quad_kind, plot_options{:});

            if nargout > 0
                varargout{1} = out;
            end

        case 'prepare'
            require_nargs(varargin, 1, 'prepare');

            experiment = varargin{1};
            prepare_options = varargin(2:end);

            script_file = write_peak_jobs(experiment, prepare_options{:});

            if nargout > 0
                varargout{1} = script_file;
            end

        case 'worker'
            % Internal entry point used only by the generated shell script.
            %
            % Arguments:
            %   experiment, policy, i_deg, i_nsub, quad_kind, i_tol, method, outfile
            require_nargs(varargin, 8, 'worker');

            run_peak_case(varargin{1}, varargin{2}, varargin{3}, varargin{4}, ...
                          varargin{5}, varargin{6}, varargin{7}, varargin{8});

        otherwise
            error('measure_stiffness_peak_memory:action', 'Unknown action "%s".', action);
    end
end


function results = collect_experiment(experiment, policy, quad_kind)
    [results_gp, results_interpolation, results_projection] = collect_peaks(experiment, policy, quad_kind);

    results = struct;
    results.config = stiffness_peak_config(experiment);
    results.policy = lower(char(policy));
    results.quad_kind = lower(char(quad_kind));
    results.gp = results_gp;
    results.interpolation = results_interpolation;
    results.projection = results_projection;
end


function results = collect_all_experiments()
    % Every slice is collected for which measurements exist. An experiment
    % that has never been run is skipped instead of producing empty entries.
    experiments = {'beam', 'flag', 'rotor', 'singularity'};
    policies = {'default', 'refined'};
    quad_kinds = {'practical', 'exact'};

    results = struct;

    for i_exp = 1:numel(experiments)

        experiment = experiments{i_exp};
        cfg = stiffness_peak_config(experiment);

        if ~exist(cfg.outdir, 'dir')
            fprintf('No measurements for %s, skipped.\n', experiment);
            continue;
        end

        for i_policy = 1:numel(policies)

            for i_quad = 1:numel(quad_kinds)
                results.(experiment).(policies{i_policy}).(quad_kinds{i_quad}) = ...
                    collect_experiment(experiment, policies{i_policy}, quad_kinds{i_quad});
            end
        end
    end
end


function script_file = write_peak_jobs(experiment, varargin)
%WRITE_PEAK_JOBS  Generate the terminal script for the isolated sweep.
%
% One MATLAB call per case. Existing result files are skipped by the shell
% script, so an interrupted sweep can simply be restarted.

    p = inputParser;
    p.addParameter('policies', {'default', 'refined'});
    p.addParameter('quads', {'practical', 'exact'});
    p.addParameter('methods', {'gp', 'interpolation', 'projection'});
    p.parse(varargin{:});

    policies = normalize_cellstr(p.Results.policies);
    quad_kinds = normalize_cellstr(p.Results.quads);
    methods = normalize_cellstr(p.Results.methods);

    validate_members(policies, {'default', 'refined'}, 'policy');
    validate_members(quad_kinds, {'practical', 'exact'}, 'quadrature kind');
    validate_members(methods, {'gp', 'interpolation', 'projection'}, 'method');

    cfg = stiffness_peak_config(experiment);
    here = fileparts(mfilename('fullpath'));

    if ~exist(cfg.outdir, 'dir')
        mkdir(cfg.outdir);
    end

    % Snapshot the current MATLAB search path. Each child restores it before
    % constructing the geometry and the solution space.
    saved_path = path;
    save(fullfile(here, 'session_path.mat'), 'saved_path');

    matlab_executable = fullfile(matlabroot, 'bin', 'matlab');

    jobs = build_jobs(cfg, policies, quad_kinds, methods);

    script_file = fullfile(here, sprintf('run_peak_jobs_%s_stiff.sh', cfg.experiment));

    fid = fopen(script_file, 'w');

    if fid < 0
        error('measure_stiffness_peak_memory:writeJobs', 'Could not create "%s".', script_file);
    end

    fprintf(fid, '#!/usr/bin/env bash\n');
    fprintf(fid, 'set -u\n');
    fprintf(fid, 'cd "%s" || { echo "cannot cd to code dir"; exit 1; }\n', here);
    fprintf(fid, 'ML="%s"\n', matlab_executable);
    fprintf(fid, 'if [ ! -x "$ML" ]; then echo "MATLAB not found at $ML"; exit 1; fi\n');
    fprintf(fid, 'mkdir -p "%s"\n\n', cfg.outdir);
    fprintf(fid, 'echo "Stiffness %s: %d jobs"\n\n', cfg.experiment, numel(jobs));

    for k = 1:numel(jobs)

        job = jobs(k);
        tag = job_tag(job);

        outfile = fullfile(cfg.outdir, sprintf('peak_%s.mat', tag));
        logfile = fullfile(cfg.outdir, sprintf('peak_%s.log', tag));

        % The child enters the worker through this file.
        statement = sprintf('measure_stiffness_peak_memory(''worker'',''%s'',''%s'',%d,%d,''%s'',%d,''%s'',''%s'');', ...
                            cfg.experiment, job.policy, job.i_deg, job.i_nsub, job.quad_kind, job.i_tol, ...
                            job.method, outfile);

        fprintf(fid, 'if [ ! -f "%s" ]; then\n', outfile);
        fprintf(fid, '  echo "[%d/%d] %s"\n', k, numel(jobs), tag);
        fprintf(fid, '  "$ML" -nodisplay -batch "%s" > "%s" 2>&1\n', statement, logfile);
        fprintf(fid, '  [ -f "%s" ] || echo "   -> NO OUTPUT, see %s"\n', outfile, logfile);
        fprintf(fid, 'fi\n');
    end

    fclose(fid);

    try
        fileattrib(script_file, '+x');
    catch
        % Not critical, the script is started with "bash".
    end

    fprintf('\nPrepared the stiffness peak-memory sweep.\n');
    fprintf('Experiment : %s\n', cfg.experiment);
    fprintf('Jobs       : %d\n', numel(jobs));
    fprintf('Shell file : %s\n', script_file);
    fprintf('Path file  : %s\n', fullfile(here, 'session_path.mat'));
    fprintf('\nRun in a terminal:\n');
    fprintf('  bash "%s"\n\n', script_file);
end


function jobs = build_jobs(cfg, policies, quad_kinds, methods)
%BUILD_JOBS  Job grid of the isolated sweep.
%
% The GeoPDEs reference does not depend on the rank tolerance, but it does
% depend on the policy: the exact rule is built from the auxiliary degree of
% that policy, and the refinement grid may differ between policies. It is
% therefore emitted once per policy, quadrature kind, degree and refinement
% level.

    requested = @(method) any(strcmp(method, methods));

    jobs = repmat(make_job('', '', '', 1, 1, 1), 0, 1);

    for i_policy = 1:numel(policies)

        policy = policies{i_policy};

        for i_quad_kind = 1:numel(quad_kinds)

            quad_kind = quad_kinds{i_quad_kind};
            i_quad = quadrature_index(quad_kind);

            for i_deg = 1:size(cfg.degree, 1)

                refinements = cfg.n_sub.(policy){i_deg, i_quad};

                for i_nsub = 1:size(refinements, 1)

                    if requested('gp')
                        jobs(end+1,1) = make_job('gp', policy, quad_kind, i_deg, i_nsub, 1); %#ok<AGROW>
                    end

                    for i_tol = 1:numel(cfg.tol)

                        if requested('interpolation') && cfg.interpolation_supported
                            jobs(end+1,1) = make_job('interpolation', policy, quad_kind, i_deg, i_nsub, i_tol); %#ok<AGROW>
                        end

                        if requested('projection')
                            jobs(end+1,1) = make_job('projection', policy, quad_kind, i_deg, i_nsub, i_tol); %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    if requested('interpolation') && ~cfg.interpolation_supported
        fprintf('Interpolation jobs are omitted for %s: the interpolation system is singular.\n', cfg.name);
    end
end


function job = make_job(method, policy, quad_kind, i_deg, i_nsub, i_tol)
    job = struct;
    job.method = method;
    job.policy = policy;
    job.quad_kind = quad_kind;
    job.i_deg = i_deg;
    job.i_nsub = i_nsub;
    job.i_tol = i_tol;
end


function tag = job_tag(job)
    tag = sprintf('%s_%s_%s_%d_%d_%d', job.method, job.policy, job.quad_kind, ...
                  job.i_deg, job.i_nsub, job.i_tol);
end


function run_peak_case(experiment, policy, i_deg, i_nsub, quad_kind, i_tol, method, outfile)
%RUN_PEAK_CASE  Measure one stiffness assembly in a fresh MATLAB process.
%
% Do not insert any operation between reset_peak_rss() and the final
% vmhwm_bytes() call. Everything that allocates in between enters the
% reported peak.

    restore_session_path();

    policy = lower(char(policy));
    quad_kind = lower(char(quad_kind));
    method = lower(char(method));

    results = struct;
    results.ok = false;
    results.experiment = experiment;
    results.policy = policy;
    results.quad_kind = quad_kind;
    results.method = method;
    results.i_deg = i_deg;
    results.i_nsub = i_nsub;
    results.i_tol = i_tol;

    try
        % The scratch space of the multithreaded BLAS grows with the number
        % of computational threads and would otherwise enter the peak.
        maxNumCompThreads(1);

        cfg = stiffness_peak_config(experiment);
        geometry = cfg.make_geometry();

        p_geo = geometry.nurbs.order(:).' - 1;
        p_sol = cfg.degree(i_deg, :);

        i_quad = quadrature_index(quad_kind);
        n_sub = cfg.n_sub.(policy){i_deg, i_quad}(i_nsub, :);
        tol = cfg.tol(i_tol);

        is_low_rank = ~strcmp(method, 'gp');

        % -----------------------------------------------------------------
        % Auxiliary degree and quadrature rule
        %
        % The auxiliary degree enters the exact rule of both the low-rank
        % methods and the GeoPDEs reference, because the reference is
        % assembled with the same rule as the methods it is compared
        % against, exactly as in experiments_paper_stiffness.m.
        % -----------------------------------------------------------------
        default_aux_space = [];

        switch policy

            case 'default'
                % Eq. (57), a fixed geometry-based space.
                default_aux_space = make_paper_default_aux_space(geometry);
                p_aux = default_aux_space.degree(:).';

            case 'refined'
                % Eq. (58). The degree is known in advance; the space itself
                % needs the solution space and is built below.
                p_aux = 3 * p_sol - 1;

            otherwise
                error('measure_stiffness_peak_memory:policy', ...
                      'Unknown auxiliary-space policy "%s".', policy);
        end

        switch quad_kind

            case 'practical'
                n_quad = practical_nquad(numel(p_geo));

            case 'exact'
                n_quad = sufficient_stiffness_nquad(p_sol, p_aux, p_geo);

            otherwise
                error('measure_stiffness_peak_memory:quadrature', ...
                      'Unknown quadrature kind "%s".', quad_kind);
        end

        method_data = build_method_data(p_sol, n_sub, n_quad);

        [space, msh] = initialize_laplace(geometry, method_data);

        % Only the low-rank methods use an auxiliary space. The
        % solution-aligned one is available only now.
        aux_space = [];

        if is_low_rank

            switch policy

                case 'default'
                    aux_space = default_aux_space;

                case 'refined'
                    aux_space = make_solution_aux_space(space);
            end
        end

        % A fixed seed makes the random initial guesses of the AMEn solver,
        % and therefore the resulting TT ranks and the peak, reproducible.
        rng('default');

        % -----------------------------------------------------------------
        % MEMORY MEASUREMENT
        % -----------------------------------------------------------------

        base = vmrss_bytes();

        if ~isfinite(base)
            error('measure_stiffness_peak_memory:vmrss', 'Could not read VmRSS from /proc/self/status.');
        end

        hwm_after_reset = reset_peak_rss(base);

        t_assembly = tic;

        switch method

            case 'gp'
                K = op_gradu_gradv_tp(space, space, msh, cfg.c_diff);

            case 'interpolation'
                low_rank_data = build_low_rank_data(tol, aux_space);

                H = cfg.interpolate_system(geometry, low_rank_data);

                K = cfg.assemble_interpolation(H, space, low_rank_data, method_data);

            case 'projection'
                K = cfg.assemble_projection(tol, geometry, space, n_quad, cfg.projection_refinement, aux_space);

            otherwise
                error('measure_stiffness_peak_memory:method', 'Unknown method "%s".', method);
        end

        hwm = vmhwm_bytes();

        % -----------------------------------------------------------------
        % POST-PROCESSING
        %
        % Only operations below this point may allocate additional memory;
        % they cannot influence the recorded VmHWM.
        % -----------------------------------------------------------------

        assembly_time = toc(t_assembly);

        results.peak_memory = hwm - base;
        results.peak_abs = hwm;
        results.base_rss = base;
        results.hwm_after_reset = hwm_after_reset;
        results.K_memory = RecursiveSize(K);
        results.K_nnz = count_tt_nnz(K);
        results.ndof = space.ndof;
        results.assembly_time = assembly_time;
        results.degree = p_sol;
        results.n_sub = n_sub;
        results.nquad = n_quad;
        results.tol = tol;

        if isempty(aux_space)
            % GeoPDEs uses no auxiliary space.
            results.aux_degree = nan(1, numel(p_geo));
            results.aux_ndof_dir = nan(1, numel(p_geo));
        else
            results.aux_degree = aux_space.degree(:).';
            results.aux_ndof_dir = aux_space.ndof_dir(:).';
        end

        results.ok = true;

        fprintf('OK  %-11s %-7s %-13s q=%-9s deg=%d nsub=%d tol=%d  nquad=[%s]  peak=%.1f MB  ndof=%d\n', ...
                experiment, policy, method, quad_kind, i_deg, i_nsub, i_tol, ...
                strtrim(sprintf('%d ', n_quad)), results.peak_memory/1e6, results.ndof);

    catch ME
        results.ok = false;
        results.err = getReport(ME, 'extended', 'hyperlinks', 'off');

        fprintf(2, 'FAIL %s %s %s %s: %s\n', experiment, policy, quad_kind, method, results.err);
    end

    save(outfile, '-struct', 'results');
end


function restore_session_path()
    % The snapshot written by 'prepare' lies next to this file, so that the
    % child finds GeoPDEs, the TT toolbox, the assembly routines and the
    % geometry data.
    snapshot = fullfile(fileparts(mfilename('fullpath')), 'session_path.mat');

    if ~exist(snapshot, 'file')
        return;
    end

    stored = load(snapshot);

    if ~isfield(stored, 'saved_path') || isempty(stored.saved_path)
        return;
    end

    directories = strsplit(stored.saved_path, pathsep);
    directories = directories(~cellfun(@isempty, directories));

    addpath(directories{:});
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


function low_rank_data = build_low_rank_data(tol, aux_space)
    low_rank_data = struct;
    low_rank_data.mass = 0;
    low_rank_data.stiffness = 1;
    low_rank_data.TT_interpolation = 1;
    low_rank_data.boundary_conditions = 'Dirichlet';
    low_rank_data.geometry_format = 'B-Splines';
    low_rank_data.rankTol = tol;
    low_rank_data.aux_space = aux_space;
end


function n = count_tt_nnz(K)
    % The TT toolbox does not overload nnz for every operator class, so the
    % core-wise count is used as a fallback.
    try
        n = nnz(K);
    catch
        try
            cores = core2cell(K);
            n = sum(cellfun(@nnz, cores));
        catch
            n = NaN;
        end
    end
end


function bytes = vmrss_bytes()
    % Resident set size at this moment.
    bytes = status_field('VmRSS:');
end


function bytes = vmhwm_bytes()
    % Largest resident set size since the last reset.
    bytes = status_field('VmHWM:');
end


function hwm_after_reset = reset_peak_rss(base)
%RESET_PEAK_RSS  Reset the high-water resident set size of this process.
%
% Writing 5 to /proc/self/clear_refs sets VmHWM to the current VmRSS. This
% is available from Linux 4.0 on and requires CONFIG_PROC_PAGE_MONITOR. On
% a kernel that ignores the request, VmHWM would still reflect the
% construction of the discretisation and the reported peak would be far too
% large, so the reset is verified rather than assumed.

    fid = fopen('/proc/self/clear_refs', 'w');

    if fid > 0
        fprintf(fid, '5');
        fclose(fid);
    end

    hwm_after_reset = vmhwm_bytes();

    if ~isfinite(hwm_after_reset)
        error('measure_stiffness_peak_memory:vmhwm', 'Could not read VmHWM after the reset.');
    end

    % The reset brings VmHWM down to the current VmRSS. A small gap is
    % normal, a large one means the kernel ignored the request.
    tolerance = max(5 * 1024^2, 0.05 * base);

    if hwm_after_reset - base > tolerance
        error('measure_stiffness_peak_memory:peakReset', ...
              ['VmHWM was not reset to the current resident set size. ', ...
               'Baseline %.1f MB, VmHWM after the reset %.1f MB.'], base/1e6, hwm_after_reset/1e6);
    end
end


function bytes = status_field(key)
%STATUS_FIELD  Read a kB-valued field from /proc/self/status.

    contents = read_proc('/proc/self/status');

    if isempty(contents)
        bytes = NaN;
        return;
    end

    token = regexp(contents, [key, '\s*(\d+)\s*kB'], 'tokens', 'once');

    if isempty(token)
        bytes = NaN;
    else
        bytes = str2double(token{1}) * 1024;
    end
end


function contents = read_proc(file_name)
    fid = fopen(file_name, 'r');

    if fid < 0
        contents = '';
        return;
    end

    contents = fread(fid, inf, '*char').';
    fclose(fid);
end


function cfg = stiffness_peak_config(experiment)
%STIFFNESS_PEAK_CONFIG  Stiffness peak-memory experiment configuration.
%
% Single source of truth for the isolated worker, the job generator, the
% collector and the plotter. Degrees, tolerances and refinement grids agree
% with experiments_paper_stiffness.m, so that peak memory, assembly time and
% accuracy refer to the same cases.
%
% The refinement grids are indexed as
%
%   cfg.n_sub.<policy>{i_deg, i_quad}
%
% with i_quad = 1 for the practical rule and i_quad = 2 for the exact rule.

    cfg.experiment = lower(char(experiment));
    cfg.tol = [1e-5, 1e-7];
    cfg.degree = [3, 3, 3; 5, 5, 5];

    % Constant coefficient in the GeoPDEs stiffness form.
    cfg.c_diff = @(x, y, z) ones(size(x));

    % Assembly routines used in the comparison. These are the only names
    % that have to be adapted if the assembly functions are renamed.
    cfg.interpolate_system = @interpolate_system;
    cfg.assemble_interpolation = @assemble_stiffness_lowrank_interpolation;
    cfg.assemble_projection = @assemble_stiffness_lowrank_projection;

    % The projection method is given an explicit auxiliary space rather than
    % constructing one itself.
    cfg.projection_refinement = 2;

    % A collected slice pins one quadrature kind, so it has a single
    % quadrature slot per degree. This field lets the plotter treat mass and
    % stiffness results alike.
    cfg.n_quadrature = repmat({[0, 0, 0]}, size(cfg.degree, 1), 1);

    % The output folder is anchored to this file, so that the collector
    % finds the case files independently of the working directory.
    cfg.outdir = fullfile(fileparts(mfilename('fullpath')), 'peak_out', [cfg.experiment, '_stiff']);

    switch cfg.experiment

        case 'flag'
            cfg.name = 'thick flag';
            cfg.make_geometry = @make_flag_geometry;
            cfg.interpolation_supported = true;

            % The exact-quadrature experiments use one fewer refinement step
            % because of their larger computational cost.
            p3_practical = step_grid(0, 18, 3);
            p5_practical = step_grid(0, 15, 3);
            p3_exact = step_grid(0, 15, 3);
            p5_exact = step_grid(0, 12, 3);

            cfg.n_sub.default = {p3_practical, p3_exact; p5_practical, p5_exact};

            % Both auxiliary-space policies use the same solution refinements.
            cfg.n_sub.refined = cfg.n_sub.default;

        case 'rotor'
            cfg.name = 'rotor blade';
            cfg.make_geometry = @make_rotor_geometry;
            cfg.interpolation_supported = true;

            p3 = step_grid(0, 7, 1);
            p5 = step_grid(0, 5, 1);

            cfg.n_sub.default = {p3, p3; p5, p5};
            cfg.n_sub.refined = cfg.n_sub.default;

        case 'singularity'
            cfg.name = 'almost singular geometry';
            cfg.make_geometry = @make_singularity_geometry;
            cfg.interpolation_supported = true;

            % For this geometry the default and the solution-aligned
            % auxiliary spaces are tested over different refinement ranges.
            p3_default = step_grid(4, 28, 4);
            p5_default = step_grid(4, 20, 4);
            p3_refined = step_grid(4, 20, 4);
            p5_refined = step_grid(4, 16, 4);

            cfg.n_sub.default = {p3_default, p3_default; p5_default, p5_default};
            cfg.n_sub.refined = {p3_refined, p3_refined; p5_refined, p5_refined};

        case 'beam'
            cfg.name = 'twisted beam';
            cfg.make_geometry = @make_beam_geometry;

            % The interpolation system of this geometry is singular.
            cfg.interpolation_supported = false;

            p3 = step_grid(0, 7, 1);
            p5 = step_grid(0, 5, 1);

            cfg.n_sub.default = {p3, p3; p5, p5};
            cfg.n_sub.refined = cfg.n_sub.default;

        otherwise
            error('measure_stiffness_peak_memory:experiment', 'Unknown experiment "%s".', experiment);
    end
end


function refinements = step_grid(first, last, step)
    % Rows of isotropically inserted knots, [n n n].
    levels = (first:step:last).';
    refinements = [levels, levels, levels];
end


function i_quad = quadrature_index(quad_kind)
    % Column index of the refinement grids: the practical rule first, the
    % exact rule second.
    switch lower(char(quad_kind))

        case 'practical'
            i_quad = 1;

        case 'exact'
            i_quad = 2;

        otherwise
            error('measure_stiffness_peak_memory:quadrature', ...
                  'Unknown quadrature kind "%s".', quad_kind);
    end
end


function geometry = make_beam_geometry()
    data = load('geo_twisted_pipe.mat', 'geo');
    geometry = geo_load(data.geo);
end


function geometry = make_flag_geometry()
    data = load('3Dflag_GeoPDEs.mat', 'g');
    geometry = geo_load(data.g.nurbs);
end


function geometry = make_rotor_geometry()
    [geometry, ~, ~] = rotor_gen();
end


function geometry = make_singularity_geometry()
    % Trilinear map whose Jacobian determinant degenerates for c -> -1.
    % The value below keeps the determinant positive.
    c = -1.0 + 1e-5;

    knots = {[0 0 1 1], [0 0 1 1], [0 0 1 1]};
    coefs = zeros(4,2,2,2);

    for i = 1:2
        xi = i - 1;

        for j = 1:2
            eta = j - 1;

            for k = 1:2
                zeta = k - 1;

                x = (1 + c*eta*zeta) * xi;
                y = eta;
                z = zeta;
                w = 1;

                coefs(:,i,j,k) = [x*w; y*w; z*w; w];
            end
        end
    end

    geometry = geo_load(nrbmak(coefs, knots));
end


function aux_space = make_paper_default_aux_space(geometry)
%MAKE_PAPER_DEFAULT_AUX_SPACE  Geometry-based auxiliary space.
%
% For a geometry of degree p with interior knot multiplicity mu, the
% auxiliary space has
%
%   degree                 6*p - 2,
%   endpoint multiplicity  6*p - 1,
%   interior multiplicity  5*p + mu - 1.

    dim = geometry.rdim;

    aux_space = struct;
    aux_space.knots = cell(dim, 1);
    aux_space.degree = zeros(1, dim);
    aux_space.ndof_dir = zeros(1, dim);

    for d = 1:dim

        U = geometry.nurbs.knots{d}(:).';
        p = geometry.nurbs.order(d) - 1;

        [distinct_knots, ~, knot_ids] = unique(U, 'stable');
        multiplicities = accumarray(knot_ids(:), 1).';

        p_aux = 6*p - 2;
        U_aux = [];

        for r = 1:numel(distinct_knots)

            if r == 1 || r == numel(distinct_knots)
                m_aux = 6*p - 1;
            else
                m_aux = 5*p + multiplicities(r) - 1;
            end

            U_aux = [U_aux, repmat(distinct_knots(r), 1, m_aux)]; %#ok<AGROW>
        end

        aux_space.knots{d} = U_aux;
        aux_space.degree(d) = p_aux;
        aux_space.ndof_dir(d) = numel(U_aux) - p_aux - 1;
    end
end


function n_quad = sufficient_stiffness_nquad(p_sol, p_aux, p_geo)
%SUFFICIENT_STIFFNESS_NQUAD  Sufficient Gaussian rule per direction.
%
% The stiffness integrand has degree at most 2*p_sol + p_aux + 4*p_geo, and
% a Gaussian rule with n points is exact up to degree 2*n - 1.

    p_sol = p_sol(:).';
    p_aux = p_aux(:).';
    p_geo = p_geo(:).';

    if numel(p_sol) ~= numel(p_aux) || numel(p_sol) ~= numel(p_geo)
        error('measure_stiffness_peak_memory:degrees', 'The degree vectors must have the same length.');
    end

    max_integrand_degree = 2*p_sol + p_aux + 4*p_geo;
    n_quad = ceil((max_integrand_degree + 1) / 2);
end


function n_quad = practical_nquad(dim)
    n_quad = 5 * ones(1, dim);
end


function [results_gp, results_interpolation, results_projection] = collect_peaks(experiment, policy, quad_kind)
%COLLECT_PEAKS  Collect one policy and quadrature slice of one experiment.
%
% The output follows the layout of the mass results:
%
%   results_gp.<field>{i_deg, 1}
%   results_interpolation.<field>{i_deg, 1, i_tol}
%   results_projection.<field>{i_deg, 1, i_tol}
%
% Each entry is a vector over the refinement levels of that slice. The
% quadrature index is fixed to one, because the slice already pins one
% quadrature kind. Cases that were not run, or that ran out of memory, stay
% NaN, so that a missing measurement never shifts the refinement indexing.
%
% The GeoPDEs reference belongs to the slice as well: its exact rule is
% built from the auxiliary degree of the policy, so its curve is the one
% that the low-rank curves of this slice must be compared against. Fields
% that hold one value per direction, such as nquad and aux_degree, are
% reduced to their first component.

    cfg = stiffness_peak_config(experiment);

    policy = lower(char(policy));
    quad_kind = lower(char(quad_kind));

    degree_n = size(cfg.degree, 1);
    tol_n = numel(cfg.tol);
    i_quad = quadrature_index(quad_kind);

    fields = {'peak_memory', 'peak_abs', 'K_memory', 'K_nnz', 'ndof', 'assembly_time', 'nquad', 'aux_degree'};

    results_gp = allocate_results(fields, {degree_n, 1});
    results_interpolation = allocate_results(fields, {degree_n, 1, tol_n});
    results_projection = allocate_results(fields, {degree_n, 1, tol_n});

    for i_deg = 1:degree_n

        refinement_n = size(cfg.n_sub.(policy){i_deg, i_quad}, 1);

        for i_field = 1:numel(fields)

            field = fields{i_field};

            results_gp.(field){i_deg,1} = nan(1, refinement_n);

            for i_tol = 1:tol_n
                results_interpolation.(field){i_deg,1,i_tol} = nan(1, refinement_n);
                results_projection.(field){i_deg,1,i_tol} = nan(1, refinement_n);
            end
        end
    end

    files = dir(fullfile(cfg.outdir, 'peak_*.mat'));

    expression = ['^peak_(gp|interpolation|projection)_(default|refined)_', ...
                  '(practical|exact)_(\d+)_(\d+)_(\d+)\.mat$'];

    n_ok = 0;
    n_fail = 0;

    for k = 1:numel(files)

        tokens = regexp(files(k).name, expression, 'tokens', 'once');

        if isempty(tokens)
            continue;
        end

        method = tokens{1};
        file_policy = tokens{2};
        file_quad_kind = tokens{3};
        i_deg = str2double(tokens{4});
        i_nsub = str2double(tokens{5});
        i_tol = str2double(tokens{6});

        if ~strcmp(file_quad_kind, quad_kind) || ~strcmp(file_policy, policy)
            continue;
        end

        case_results = load(fullfile(files(k).folder, files(k).name));

        if ~isfield(case_results, 'ok') || ~case_results.ok
            n_fail = n_fail + 1;
            continue;
        end

        n_ok = n_ok + 1;

        for i_field = 1:numel(fields)

            field = fields{i_field};
            value = first_component(case_results, field);

            switch method

                case 'gp'
                    results_gp = put_result(results_gp, field, {i_deg,1}, i_nsub, value);

                case 'interpolation'
                    results_interpolation = put_result(results_interpolation, field, {i_deg,1,i_tol}, i_nsub, value);

                case 'projection'
                    results_projection = put_result(results_projection, field, {i_deg,1,i_tol}, i_nsub, value);
            end
        end
    end

    fprintf('collect_peaks(%s, %s, %s): %d ok, %d failed, %d files in %s\n', ...
            cfg.experiment, policy, quad_kind, n_ok, n_fail, numel(files), cfg.outdir);
end


function value = first_component(case_results, field)
    % Per-direction quantities are stored as rows. Only the first direction
    % is kept, so that every collected field is a vector over refinements.
    if ~isfield(case_results, field)
        value = NaN;
        return;
    end

    stored = case_results.(field);

    if isempty(stored)
        value = NaN;
    elseif isscalar(stored)
        value = stored;
    else
        value = stored(1);
    end
end


function results = allocate_results(fields, dimensions)
    results = struct;

    for i = 1:numel(fields)
        results.(fields{i}) = cell(dimensions{:});
    end
end


function results = put_result(results, field, idx, i_nsub, value)
    try
        entry = results.(field){idx{:}};

        if i_nsub <= numel(entry)
            entry(i_nsub) = value;
            results.(field){idx{:}} = entry;
        end

    catch
        % A case file outside the configured grid, for instance from an
        % earlier version of the sweep. It is ignored.
    end
end


function out = plot_peak_scaling(results_gp, results_interpolation, results_projection, cfg, policy, quad_kind, varargin)
%PLOT_PEAK_SCALING  Peak memory against the number of degrees of freedom.
%
% One figure per degree, on doubly logarithmic axes, with an empirical
% slope fitted to each curve. The refinement range in which the GeoPDEs
% reference no longer completes is shaded.

    p = inputParser;
    p.addParameter('degrees', []);
    p.addParameter('tols', []);
    p.addParameter('metric', 'peak_memory');
    p.addParameter('unit', 'MB');
    p.addParameter('savedir', '');
    p.parse(varargin{:});

    opt = p.Results;

    [scale, unit_label] = memory_unit(opt.unit);

    if isempty(opt.degrees)
        opt.degrees = 1:size(cfg.degree, 1);
    end

    if isempty(opt.tols)
        opt.tols = 1:numel(cfg.tol);
    end

    color_gp = [0.00, 0.00, 0.00];
    color_interpolation = [0.16, 0.44, 0.78];
    color_projection = [0.85, 0.33, 0.10];

    markers = {'o', 's', '^', 'd', 'v', '>'};

    out = struct('degree', {}, 'slope_gp', {}, 'slope_interpolation', {}, 'slope_projection', {}, ...
                 'ndof_max_gp', {}, 'ndof_max_interpolation', {}, 'ndof_max_projection', {});

    for i_deg = opt.degrees

        degree_value = cfg.degree(i_deg, 1);

        fig = figure('Color', 'w', 'Name', sprintf('stiffness peak scaling, degree %d', degree_value));

        ax = axes(fig);
        hold(ax, 'on');
        grid(ax, 'on');
        box(ax, 'on');

        set(ax, 'XScale', 'log', 'YScale', 'log');

        % The fields are created in the order of the output structure, so
        % that the record can be appended to it.
        rec = struct;
        rec.degree = degree_value;
        rec.slope_gp = NaN;
        rec.slope_interpolation = nan(1, numel(opt.tols));
        rec.slope_projection = nan(1, numel(opt.tols));
        rec.ndof_max_gp = NaN;
        rec.ndof_max_interpolation = NaN;
        rec.ndof_max_projection = NaN;

        % --- GeoPDEs reference -------------------------------------------
        [x_gp, y_gp] = memory_curve(results_gp, {i_deg,1}, opt.metric, scale);

        rec.slope_gp = loglog_slope(x_gp, y_gp);
        rec.ndof_max_gp = safe_max(x_gp);

        if ~isempty(x_gp)
            plot(ax, x_gp, y_gp, '-o', 'Color', color_gp, 'LineWidth', 1.8, ...
                 'MarkerFaceColor', color_gp, 'MarkerSize', 6, ...
                 'DisplayName', sprintf('GeoPDEs  (slope %.2f)', rec.slope_gp));
        end

        % --- Low-rank methods, one curve per tolerance -------------------
        for i = 1:numel(opt.tols)

            i_tol = opt.tols(i);
            marker = markers{min(i, numel(markers))};
            tol = cfg.tol(i_tol);

            [x_int, y_int] = memory_curve(results_interpolation, {i_deg,1,i_tol}, opt.metric, scale);

            slope_interpolation = loglog_slope(x_int, y_int);

            rec.slope_interpolation(i) = slope_interpolation;
            rec.ndof_max_interpolation = max([rec.ndof_max_interpolation, safe_max(x_int)]);

            if ~isempty(x_int)
                plot(ax, x_int, y_int, ['-', marker], 'Color', color_interpolation, ...
                     'LineWidth', 1.4, 'MarkerSize', 5, ...
                     'DisplayName', sprintf('interpolation, tol=%g  (slope %.2f)', tol, slope_interpolation));
            end

            [x_proj, y_proj] = memory_curve(results_projection, {i_deg,1,i_tol}, opt.metric, scale);

            slope_projection = loglog_slope(x_proj, y_proj);

            rec.slope_projection(i) = slope_projection;
            rec.ndof_max_projection = max([rec.ndof_max_projection, safe_max(x_proj)]);

            if ~isempty(x_proj)
                plot(ax, x_proj, y_proj, ['-', marker], 'Color', color_projection, ...
                     'LineWidth', 1.4, 'MarkerSize', 5, ...
                     'DisplayName', sprintf('projection, tol=%g  (slope %.2f)', tol, slope_projection));
            end
        end

        % --- Range beyond the largest completed GeoPDEs case -------------
        all_dofs = all_ndof(results_gp, results_interpolation, results_projection, i_deg, opt.tols);

        shade_reference_limit(ax, rec.ndof_max_gp, safe_max(all_dofs));

        xlabel(ax, 'solution dofs  (space.ndof)');
        ylabel(ax, sprintf('%s  [%s]', pretty_metric(opt.metric), unit_label));
        title(ax, sprintf('%s, degree %d, %s auxiliary space, %s quadrature', ...
                          cfg.name, degree_value, policy, quad_kind));
        legend(ax, 'Location', 'northwest', 'Box', 'off');

        if ~isempty(opt.savedir)

            if ~exist(opt.savedir, 'dir')
                mkdir(opt.savedir);
            end

            file_base = fullfile(opt.savedir, sprintf('%s_stiff_%s_%s_%s_deg%02d', ...
                                 cfg.experiment, policy, quad_kind, opt.metric, degree_value));

            savefig(fig, [file_base, '.fig']);
            exportgraphics(ax, [file_base, '.png'], 'Resolution', 200);
        end

        out(end+1) = rec; %#ok<AGROW>
    end

    print_plot_summary(out, cfg, policy, quad_kind);
end


function shade_reference_limit(ax, ndof_max_gp, ndof_max_all)
    if ~isfinite(ndof_max_gp) || ~isfinite(ndof_max_all) || ndof_max_all <= ndof_max_gp
        return;
    end

    y_limits = ylim(ax);

    patch_handle = patch(ax, ...
        [ndof_max_gp, ndof_max_all, ndof_max_all, ndof_max_gp], ...
        [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
        [0.90, 0.90, 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');

    uistack(patch_handle, 'bottom');

    text(ax, sqrt(ndof_max_gp*ndof_max_all), y_limits(2), sprintf('GeoPDEs cannot complete\n(out of memory)'), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9, 'Color', [0.35, 0.35, 0.35]);

    ylim(ax, y_limits);
end


function [x, y] = memory_curve(results, idx, memory_field, scale)
    x = cell_at(results, 'ndof', idx);
    y = cell_at(results, memory_field, idx);

    x = x(:).';
    y = y(:).';

    n = min(numel(x), numel(y));

    x = x(1:n);
    y = y(1:n) * scale;

    valid = isfinite(x) & isfinite(y) & x > 0 & y > 0;

    x = x(valid);
    y = y(valid);

    [x, order] = sort(x);
    y = y(order);
end


function value = cell_at(results, field, idx)
    value = [];

    if ~isfield(results, field)
        return;
    end

    try
        value = results.(field){idx{:}};
    catch
        value = [];
    end
end


function slope = loglog_slope(x, y)
    valid = x > 0 & y > 0 & isfinite(x) & isfinite(y);

    if nnz(valid) < 2
        slope = NaN;
        return;
    end

    coefficients = polyfit(log10(x(valid)), log10(y(valid)), 1);
    slope = coefficients(1);
end


function ndof = all_ndof(results_gp, results_interpolation, results_projection, i_deg, tol_ids)
    ndof = finite_positive(cell_at(results_gp, 'ndof', {i_deg,1}));

    for i_tol = tol_ids
        ndof = [ndof, finite_positive(cell_at(results_interpolation, 'ndof', {i_deg,1,i_tol}))]; %#ok<AGROW>
        ndof = [ndof, finite_positive(cell_at(results_projection, 'ndof', {i_deg,1,i_tol}))]; %#ok<AGROW>
    end

    ndof = unique(ndof);
end


function value = finite_positive(value)
    value = value(:).';
    value = value(isfinite(value) & value > 0);
end


function value = safe_max(x)
    if isempty(x)
        value = NaN;
    else
        value = max(x);
    end
end


function [scale, unit_label] = memory_unit(unit)
    switch lower(unit)

        case 'bytes'
            scale = 1;
            unit_label = 'bytes';

        case 'gb'
            scale = 1e-9;
            unit_label = 'GB';

        otherwise
            scale = 1e-6;
            unit_label = 'MB';
    end
end


function label = pretty_metric(metric)
    switch metric

        case 'peak_memory'
            label = 'peak resident memory';

        case 'K_memory'
            label = 'stored operator memory';

        otherwise
            label = strrep(metric, '_', ' ');
    end
end


function print_plot_summary(out, cfg, policy, quad_kind)
    fprintf('\n--- stiffness peak-memory scaling, %s, %s auxiliary space, %s quadrature ---\n', ...
            cfg.name, policy, quad_kind);

    for k = 1:numel(out)

        rec = out(k);

        fprintf('degree %d:\n', rec.degree);

        fprintf('  GeoPDEs         slope %5.2f   completed up to ndof = %g\n', ...
                rec.slope_gp, rec.ndof_max_gp);

        fprintf('  interpolation   slope %s   completed up to ndof = %g\n', ...
                format_slopes(rec.slope_interpolation), rec.ndof_max_interpolation);

        fprintf('  projection      slope %s   completed up to ndof = %g\n', ...
                format_slopes(rec.slope_projection), rec.ndof_max_projection);

        if isfinite(rec.ndof_max_gp) && rec.ndof_max_projection > rec.ndof_max_gp
            fprintf('  -> the projection method reaches %.1fx more dofs than GeoPDEs before OOM.\n', ...
                    rec.ndof_max_projection / rec.ndof_max_gp);
        end
    end

    fprintf('\nA slope of one is linear growth in the number of dofs, a slope near zero is bounded memory.\n');
end


function label = format_slopes(values)
    label = ['[', strtrim(sprintf('%.2f ', values)), ']'];
end


function values = normalize_cellstr(values)
    if ischar(values) || isstring(values)
        values = cellstr(values);
    end

    values = cellfun(@(v) lower(char(v)), values, 'UniformOutput', false);
end


function validate_members(values, allowed, label)
    for k = 1:numel(values)

        if ~any(strcmp(values{k}, allowed))
            error('measure_stiffness_peak_memory:value', 'Unknown %s "%s".', label, values{k});
        end
    end
end


function require_nargs(args, required, action)
    if numel(args) < required
        error('measure_stiffness_peak_memory:arguments', ...
              'Action "%s" requires at least %d additional argument(s).', action, required);
    end
end


function print_usage()
    fprintf('\nStiffness peak-memory workflow\n');
    fprintf('------------------------------\n');
    fprintf('Use existing measurements:\n');
    fprintf('  results = measure_stiffness_peak_memory(''collect'', ''beam'', ''refined'', ''exact'');\n');
    fprintf('  results = measure_stiffness_peak_memory(''collect_all'');\n');
    fprintf('  measure_stiffness_peak_memory(''plot'', ''beam'', ''refined'', ''exact'');\n\n');
    fprintf('Run the sweep:\n');
    fprintf('  measure_stiffness_peak_memory(''prepare'', ''beam'');\n');
    fprintf('  %% then execute the printed .sh file in a terminal\n\n');
    fprintf('Inspect a configuration:\n');
    fprintf('  cfg = measure_stiffness_peak_memory(''config'', ''beam'');\n\n');
end
