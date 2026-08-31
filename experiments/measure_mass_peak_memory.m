function varargout = measure_mass_peak_memory(action, varargin)
%MEASURE_MASS_PEAK_MEMORY  Peak-memory workflow for the mass tensor.
%
% This file is the memory counterpart of experiments_paper_mass.m. That
% script measures assembly time, operator storage and accuracy; this one
% measures the transient memory required during the assembly itself.
%
% For each geometry, solution-space degree, quadrature rule, refinement
% level and rank tolerance the peak memory of
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
%   peak_out/<experiment>/peak_<method>_<ideg>_<iquad>_<insub>_<itol>.mat
%
% with method in {gp, interpolation, projection}.
%
% -------------------------------------------------------------------------
% WORKING WITH EXISTING MEASUREMENTS
% -------------------------------------------------------------------------
%
% Collect one experiment:
%
%   results = measure_mass_peak_memory('collect', 'flag');
%
% The returned structure contains
%
%   results.gp
%   results.interpolation
%   results.projection
%   results.config
%
% Collect every experiment for which measurements exist:
%
%   results = measure_mass_peak_memory('collect_all');
%
% Plot an already measured experiment:
%
%   measure_mass_peak_memory('plot', 'flag');
%
% Plot the stored operator memory instead of the peak memory:
%
%   measure_mass_peak_memory('plot', 'flag', 'metric', 'M_memory');
%
% -------------------------------------------------------------------------
% RUNNING THE SWEEP
% -------------------------------------------------------------------------
%
% Prepare the isolated jobs from an interactive session that has GeoPDEs,
% the TT toolbox and this repository on the path:
%
%   measure_mass_peak_memory('prepare', 'flag');
%
% This saves a snapshot of the search path in session_path.mat, so that
% every child process reproduces the present environment, and writes
% run_peak_jobs_flag.sh. Run that script in a terminal. The interactive
% session may be closed first, which frees the licence seat for the
% children. Completed cases are skipped, so the script may be interrupted
% and restarted.
%
% As in the original sweep, only quadrature index 1 is measured by default:
%
%   measure_mass_peak_memory('prepare', 'flag', 'quads', [1 2]);
%
% includes both quadrature rules.
%
% -------------------------------------------------------------------------
% WHAT IS MEASURED
% -------------------------------------------------------------------------
%
% For every isolated case the sequence is
%
%   1. build geometry, mesh and solution space,
%   2. record VmRSS as the baseline,
%   3. reset VmHWM through /proc/self/clear_refs,
%   4. execute exactly one assembly,
%   5. immediately read VmHWM,
%   6. only afterwards evaluate RecursiveSize and nnz,
%
% so that
%
%   peak_memory = VmHWM - VmRSS_baseline.
%
% The discretisation of step 1 is common to all three methods and is
% deliberately excluded. No timing, validation, conversion to a full matrix
% or error computation is inserted between steps 3 and 5.
%
% Supported actions
% -----------------
%   'collect'      collect existing case files for one geometry
%   'collect_all'  collect every experiment with existing case files
%   'plot'         collect and plot one geometry
%   'prepare'      generate the isolated shell jobs
%   'config'       return the experiment configuration
%   'worker'       internal mode used by the generated child jobs
%   'help'         print a short usage summary
%
% Supported experiments
% ---------------------
%   'beam', 'flag', 'rotor', 'singularity'.

    if nargin == 0
        action = 'help';
    end

    action = lower(char(action));

    switch action

        case 'help'
            print_usage();

        case 'config'
            require_nargs(varargin, 1, 'config');

            varargout{1} = peak_config(varargin{1});

        case 'collect'
            require_nargs(varargin, 1, 'collect');

            varargout{1} = collect_experiment(varargin{1});

        case 'collect_all'
            varargout{1} = collect_all_experiments();

        case 'plot'
            require_nargs(varargin, 1, 'plot');

            experiment = varargin{1};
            plot_options = varargin(2:end);

            cfg = peak_config(experiment);

            [results_gp, results_interpolation, results_projection] = collect_peaks(experiment);

            out = plot_peak_scaling(results_gp, results_interpolation, results_projection, cfg, plot_options{:});

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
            %   experiment, i_deg, i_quad, i_nsub, i_tol, method, outfile
            require_nargs(varargin, 7, 'worker');

            run_peak_case(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7});

        otherwise
            error('measure_mass_peak_memory:action', 'Unknown action "%s".', action);
    end
end


function results = collect_experiment(experiment)
    [results_gp, results_interpolation, results_projection] = collect_peaks(experiment);

    results = struct;
    results.config = peak_config(experiment);
    results.gp = results_gp;
    results.interpolation = results_interpolation;
    results.projection = results_projection;
end


function results = collect_all_experiments()
    % Every experiment is collected for which measurements exist. An
    % experiment that has never been run is skipped instead of producing an
    % empty entry.
    experiments = {'beam', 'flag', 'rotor', 'singularity'};

    results = struct;

    for i = 1:numel(experiments)

        experiment = experiments{i};
        cfg = peak_config(experiment);

        if ~exist(cfg.outdir, 'dir')
            fprintf('No measurements for %s, skipped.\n', experiment);
            continue;
        end

        results.(experiment) = collect_experiment(experiment);
    end
end


function script_file = write_peak_jobs(experiment, varargin)
%WRITE_PEAK_JOBS  Generate the terminal script for the isolated sweep.
%
% One MATLAB call per case. Existing result files are skipped by the shell
% script, so an interrupted sweep can simply be restarted.

    p = inputParser;
    p.addParameter('quads', 1);
    p.addParameter('methods', {'gp', 'interpolation', 'projection'});
    p.parse(varargin{:});

    quads = p.Results.quads;
    methods = p.Results.methods;

    if ischar(methods)
        methods = {methods};
    end

    cfg = peak_config(experiment);
    here = fileparts(mfilename('fullpath'));

    if ~exist(cfg.outdir, 'dir')
        mkdir(cfg.outdir);
    end

    % Snapshot the current MATLAB search path. Each child restores it before
    % constructing the geometry and the solution space.
    saved_path = path;
    save(fullfile(here, 'session_path.mat'), 'saved_path');

    matlab_executable = fullfile(matlabroot, 'bin', 'matlab');

    jobs = build_jobs(cfg, quads, methods);

    script_file = fullfile(here, sprintf('run_peak_jobs_%s.sh', cfg.experiment));

    fid = fopen(script_file, 'w');

    if fid < 0
        error('measure_mass_peak_memory:writeJobs', 'Could not create "%s".', script_file);
    end

    fprintf(fid, '#!/usr/bin/env bash\n');
    fprintf(fid, 'set -u\n');
    fprintf(fid, 'cd "%s" || { echo "cannot cd to code dir"; exit 1; }\n', here);
    fprintf(fid, 'ML="%s"\n', matlab_executable);
    fprintf(fid, 'if [ ! -x "$ML" ]; then echo "MATLAB not found at $ML"; exit 1; fi\n');
    fprintf(fid, 'mkdir -p "%s"\n\n', cfg.outdir);
    fprintf(fid, 'echo "Experiment %s: %d jobs"\n\n', cfg.experiment, numel(jobs));

    for k = 1:numel(jobs)

        job = jobs(k);
        tag = job_tag(job);

        outfile = fullfile(cfg.outdir, sprintf('peak_%s.mat', tag));
        logfile = fullfile(cfg.outdir, sprintf('peak_%s.log', tag));

        % The child enters the worker through this file.
        statement = sprintf('measure_mass_peak_memory(''worker'',''%s'',%d,%d,%d,%d,''%s'',''%s'');', ...
                            cfg.experiment, job.i_deg, job.i_quad, job.i_nsub, job.i_tol, job.method, outfile);

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

    fprintf('\nPrepared the mass peak-memory sweep.\n');
    fprintf('Experiment : %s\n', cfg.experiment);
    fprintf('Jobs       : %d\n', numel(jobs));
    fprintf('Shell file : %s\n', script_file);
    fprintf('Path file  : %s\n', fullfile(here, 'session_path.mat'));
    fprintf('\nRun in a terminal:\n');
    fprintf('  bash "%s"\n\n', script_file);
end


function jobs = build_jobs(cfg, quads, methods)
%BUILD_JOBS  Job grid of the isolated sweep.
%
% The GeoPDEs reference does not depend on the rank tolerance and is
% therefore measured once per degree, quadrature rule and refinement level.

    requested = @(method) any(strcmp(method, methods));

    jobs = repmat(make_job('', 1, 1, 1, 1), 0, 1);

    for i_deg = 1:size(cfg.degree, 1)

        available_quads = 1:size(cfg.n_quadrature{i_deg}, 1);
        selected_quads = intersect(quads, available_quads);
        selected_quads = selected_quads(:).';

        for i_quad = selected_quads

            for i_nsub = 1:size(cfg.n_sub{i_deg, i_quad}, 1)

                if requested('gp')
                    jobs(end+1,1) = make_job('gp', i_deg, i_quad, i_nsub, 1); %#ok<AGROW>
                end

                for i_tol = 1:numel(cfg.tol)

                    if requested('interpolation')
                        jobs(end+1,1) = make_job('interpolation', i_deg, i_quad, i_nsub, i_tol); %#ok<AGROW>
                    end

                    if requested('projection')
                        jobs(end+1,1) = make_job('projection', i_deg, i_quad, i_nsub, i_tol); %#ok<AGROW>
                    end
                end
            end
        end
    end
end


function job = make_job(method, i_deg, i_quad, i_nsub, i_tol)
    job = struct;
    job.method = method;
    job.i_deg = i_deg;
    job.i_quad = i_quad;
    job.i_nsub = i_nsub;
    job.i_tol = i_tol;
end


function tag = job_tag(job)
    tag = sprintf('%s_%d_%d_%d_%d', job.method, job.i_deg, job.i_quad, job.i_nsub, job.i_tol);
end


function run_peak_case(experiment, i_deg, i_quad, i_nsub, i_tol, method, outfile)
%RUN_PEAK_CASE  Measure one mass assembly in a fresh MATLAB process.
%
% Do not insert any operation between reset_peak_rss() and the final
% vmhwm_bytes() call. Everything that allocates in between enters the
% reported peak.

    restore_session_path();

    results = struct;
    results.ok = false;
    results.experiment = experiment;
    results.method = method;
    results.i_deg = i_deg;
    results.i_quad = i_quad;
    results.i_nsub = i_nsub;
    results.i_tol = i_tol;

    try
        % The scratch space of the multithreaded BLAS grows with the number
        % of computational threads and would otherwise enter the peak.
        maxNumCompThreads(1);

        cfg = peak_config(experiment);
        geometry = cfg.make_geometry();

        degree = cfg.degree(i_deg, :);
        n_sub = cfg.n_sub{i_deg, i_quad}(i_nsub, :);
        n_quad = cfg.n_quadrature{i_deg}(i_quad, :);
        tol = cfg.tol(i_tol);

        method_data = build_method_data(degree, n_sub, n_quad);

        [space, msh] = initialize_laplace(geometry, method_data);

        % -----------------------------------------------------------------
        % MEMORY MEASUREMENT
        % -----------------------------------------------------------------

        base = vmrss_bytes();
        reset_peak_rss();

        switch method

            case 'gp'
                M = op_u_v_tp(space, space, msh, cfg.c_diff);

            case 'interpolation'
                low_rank_data = build_low_rank_data(tol);

                H = cfg.interpolate_system(geometry, low_rank_data);

                M = cfg.assemble_interpolation(H, space, low_rank_data, method_data);

            case 'projection'
                M = cfg.assemble_projection(tol, geometry, space, n_quad);

            otherwise
                error('measure_mass_peak_memory:method', 'Unknown method "%s".', method);
        end

        hwm = vmhwm_bytes();

        % -----------------------------------------------------------------
        % POST-PROCESSING
        %
        % Only operations below this point may allocate additional memory;
        % they cannot influence the recorded VmHWM.
        % -----------------------------------------------------------------

        results.peak_memory = hwm - base;
        results.peak_abs = hwm;
        results.base_rss = base;
        results.M_memory = RecursiveSize(M);
        results.M_nnz = nnz(M);
        results.ndof = space.ndof;
        results.ok = true;

        fprintf('OK  %-11s %-13s deg=%d quad=%d nsub=%d tol=%d  peak=%.1f MB  ndof=%d\n', ...
                experiment, method, i_deg, i_quad, i_nsub, i_tol, results.peak_memory/1e6, results.ndof);

    catch ME
        results.ok = false;
        results.err = getReport(ME, 'extended', 'hyperlinks', 'off');

        fprintf(2, 'FAIL %s %s: %s\n', experiment, method, results.err);
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


function low_rank_data = build_low_rank_data(tol)
    low_rank_data = struct;
    low_rank_data.mass = 1;
    low_rank_data.stiffness = 0;
    low_rank_data.TT_interpolation = 1;
    low_rank_data.boundary_conditions = 'Dirichlet';
    low_rank_data.geometry_format = 'B-Splines';
    low_rank_data.rankTol = tol;
end


function bytes = vmrss_bytes()
    % Resident set size at this moment.
    bytes = status_field('VmRSS:');
end


function bytes = vmhwm_bytes()
    % Largest resident set size since the last reset.
    bytes = status_field('VmHWM:');
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


function reset_peak_rss()
%RESET_PEAK_RSS  Reset the high-water resident set size of this process.
%
% Writing 5 to /proc/self/clear_refs sets VmHWM to the current VmRSS. This
% is available from Linux 4.0 on and requires CONFIG_PROC_PAGE_MONITOR.
% Without the reset, VmHWM would still reflect the construction of the
% discretisation. A failing write is treated as harmless, as in the
% original implementation.

    fid = fopen('/proc/self/clear_refs', 'w');

    if fid > 0
        fprintf(fid, '5');
        fclose(fid);
    end
end


function cfg = peak_config(experiment)
%PEAK_CONFIG  Mass peak-memory experiment configuration.
%
% Single source of truth for the isolated worker, the job generator, the
% collector and the plotter. The grids agree with experiments_paper_mass.m,
% so that peak memory, assembly time and accuracy refer to the same cases.

    cfg.experiment = lower(experiment);
    cfg.tol = [1e-3, 1e-5, 1e-7];

    % Constant coefficient in the GeoPDEs mass form.
    cfg.c_diff = @(x, y, z) ones(size(x));

    % Assembly routines used in the comparison. These are the only names
    % that have to be adapted if the assembly functions are renamed.
    cfg.interpolate_system = @interpolate_system;
    cfg.assemble_interpolation = @assemble_mass_lowrank_interpolation;
    cfg.assemble_projection = @assemble_mass_lowrank_projection;

    % The output folder is anchored to this file, so that the collector
    % finds the case files independently of the working directory.
    cfg.outdir = fullfile(fileparts(mfilename('fullpath')), 'peak_out', cfg.experiment);

    switch cfg.experiment

        case 'beam'
            cfg.name = 'twisted beam';
            cfg.degree = [3, 3, 3; 5, 5, 5];
            cfg.n_quadrature = {[5, 5, 5; 6, 6, 6], [5, 5, 5; 8, 8, 8]};
            cfg.n_sub = broadcast_nsub({step_grid(0, 7, 1), step_grid(0, 5, 1)}, cfg.n_quadrature);
            cfg.make_geometry = @make_beam_geometry;

        case 'flag'
            cfg.name = 'thick flag';
            cfg.degree = [3, 3, 3; 5, 5, 5];
            cfg.n_quadrature = {[5, 5, 5; 6, 6, 6], [5, 5, 5; 8, 8, 8]};

            % The exact quadrature experiments use one fewer refinement
            % step because of their larger computational cost.
            cfg.n_sub = cell(2,2);
            cfg.n_sub{1,1} = step_grid(0, 18, 3);
            cfg.n_sub{2,1} = step_grid(0, 15, 3);
            cfg.n_sub{1,2} = step_grid(0, 15, 3);
            cfg.n_sub{2,2} = step_grid(0, 12, 3);

            cfg.make_geometry = @make_flag_geometry;

        case 'rotor'
            cfg.name = 'rotor blade';
            cfg.degree = [3, 3, 3; 5, 5, 5];
            cfg.n_quadrature = {[5, 5, 5; 6, 6, 6], [5, 5, 5; 8, 8, 8]};
            cfg.n_sub = broadcast_nsub({step_grid(0, 7, 1), step_grid(0, 5, 1)}, cfg.n_quadrature);
            cfg.make_geometry = @make_rotor_geometry;

        case 'singularity'
            cfg.name = 'singular geometry';
            cfg.degree = [3, 3, 3; 5, 5, 5];
            cfg.n_quadrature = {[5, 5, 5], [5, 5, 5; 7, 7, 7]};
            cfg.n_sub = broadcast_nsub({step_grid(4, 28, 4), step_grid(4, 20, 4)}, cfg.n_quadrature);
            cfg.make_geometry = @make_singularity_geometry;

        otherwise
            error('measure_mass_peak_memory:experiment', 'Unknown experiment "%s".', experiment);
    end

    cfg.n_sub_n = cellfun(@(x) size(x,1), cfg.n_sub);
end


function refinements = step_grid(first, last, step)
    % Rows of isotropically inserted knots, [n n n].
    levels = (first:step:last).';
    refinements = [levels, levels, levels];
end


function n_sub = broadcast_nsub(refinement_by_degree, n_quadrature)
%BROADCAST_NSUB  Use the same refinement grid for every quadrature rule.

    degree_n = numel(refinement_by_degree);
    n_quad_max = max(cellfun(@(q) size(q,1), n_quadrature));

    n_sub = cell(degree_n, n_quad_max);

    for i_deg = 1:degree_n

        for i_quad = 1:size(n_quadrature{i_deg}, 1)
            n_sub{i_deg, i_quad} = refinement_by_degree{i_deg};
        end
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


function [results_gp, results_interpolation, results_projection] = collect_peaks(experiment)
%COLLECT_PEAKS  Collect the isolated case files of one experiment.
%
% The output follows the layout of experiments_paper_mass.m:
%
%   results_gp.<field>{i_deg, i_quad}
%   results_interpolation.<field>{i_deg, i_quad, i_tol}
%   results_projection.<field>{i_deg, i_quad, i_tol}
%
% Each entry is a vector over the refinement levels. Cases that were not
% run, or that ran out of memory, stay NaN, so that a missing measurement
% never shifts the refinement indexing.

    cfg = peak_config(experiment);

    degree_n = size(cfg.degree, 1);
    tol_n = numel(cfg.tol);
    n_quad_max = max(cellfun(@(q) size(q,1), cfg.n_quadrature));

    fields = {'peak_memory', 'peak_abs', 'M_memory', 'M_nnz', 'ndof'};

    results_gp = allocate_results(fields, {degree_n, n_quad_max});
    results_interpolation = allocate_results(fields, {degree_n, n_quad_max, tol_n});
    results_projection = allocate_results(fields, {degree_n, n_quad_max, tol_n});

    for i_deg = 1:degree_n

        for i_quad = 1:size(cfg.n_quadrature{i_deg}, 1)

            refinement_n = size(cfg.n_sub{i_deg, i_quad}, 1);

            for i_field = 1:numel(fields)

                field = fields{i_field};

                results_gp.(field){i_deg,i_quad} = nan(1, refinement_n);

                for i_tol = 1:tol_n
                    results_interpolation.(field){i_deg,i_quad,i_tol} = nan(1, refinement_n);
                    results_projection.(field){i_deg,i_quad,i_tol} = nan(1, refinement_n);
                end
            end
        end
    end

    files = dir(fullfile(cfg.outdir, 'peak_*.mat'));

    expression = '^peak_(gp|interpolation|projection)_(\d+)_(\d+)_(\d+)_(\d+)\.mat$';

    n_ok = 0;
    n_fail = 0;

    for k = 1:numel(files)

        tokens = regexp(files(k).name, expression, 'tokens', 'once');

        if isempty(tokens)
            continue;
        end

        method = tokens{1};
        i_deg = str2double(tokens{2});
        i_quad = str2double(tokens{3});
        i_nsub = str2double(tokens{4});
        i_tol = str2double(tokens{5});

        case_results = load(fullfile(files(k).folder, files(k).name));

        if ~isfield(case_results, 'ok') || ~case_results.ok
            n_fail = n_fail + 1;
            continue;
        end

        n_ok = n_ok + 1;

        for i_field = 1:numel(fields)

            field = fields{i_field};
            value = NaN;

            if isfield(case_results, field)
                value = case_results.(field);
            end

            switch method

                case 'gp'
                    results_gp = put_result(results_gp, field, {i_deg,i_quad}, i_nsub, value);

                case 'interpolation'
                    results_interpolation = put_result(results_interpolation, field, {i_deg,i_quad,i_tol}, i_nsub, value);

                case 'projection'
                    results_projection = put_result(results_projection, field, {i_deg,i_quad,i_tol}, i_nsub, value);
            end
        end
    end

    fprintf('collect_peaks(%s): %d ok, %d failed, %d files in %s\n', ...
            cfg.experiment, n_ok, n_fail, numel(files), cfg.outdir);
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


function out = plot_peak_scaling(results_gp, results_interpolation, results_projection, cfg, varargin)
%PLOT_PEAK_SCALING  Peak memory against the number of degrees of freedom.
%
% One figure per degree, on doubly logarithmic axes, with an empirical
% slope fitted to each curve. The refinement range in which the GeoPDEs
% reference no longer completes is shaded.

    p = inputParser;
    p.addParameter('degrees', []);
    p.addParameter('quad_idx', []);
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

        if isempty(opt.quad_idx)
            i_quad = pick_quad(results_gp, results_projection, i_deg, cfg, opt.metric);
        else
            i_quad = min(opt.quad_idx, size(cfg.n_quadrature{i_deg}, 1));
        end

        fig = figure('Color', 'w', 'Name', sprintf('peak scaling, degree %d', degree_value));

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
        rec.ndof_max_interpolation = 0;
        rec.ndof_max_projection = 0;

        % --- GeoPDEs reference -------------------------------------------
        [x_gp, y_gp] = memory_curve(results_gp, {i_deg,i_quad}, opt.metric, scale);

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

            [x_int, y_int] = memory_curve(results_interpolation, {i_deg,i_quad,i_tol}, opt.metric, scale);

            slope_interpolation = loglog_slope(x_int, y_int);

            rec.slope_interpolation(i) = slope_interpolation;
            rec.ndof_max_interpolation = max([rec.ndof_max_interpolation, safe_max(x_int)]);

            if ~isempty(x_int)
                plot(ax, x_int, y_int, ['-', marker], 'Color', color_interpolation, ...
                     'LineWidth', 1.4, 'MarkerSize', 5, ...
                     'DisplayName', sprintf('interpolation, tol=%g  (slope %.2f)', tol, slope_interpolation));
            end

            [x_proj, y_proj] = memory_curve(results_projection, {i_deg,i_quad,i_tol}, opt.metric, scale);

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
        all_dofs = all_ndof(results_gp, results_interpolation, results_projection, i_deg, i_quad, opt.tols);

        shade_reference_limit(ax, rec.ndof_max_gp, safe_max(all_dofs));

        xlabel(ax, 'solution dofs  (space.ndof)');
        ylabel(ax, sprintf('%s  [%s]', pretty_metric(opt.metric), unit_label));
        title(ax, sprintf('%s, degree %d (quad idx %d)', cfg.name, degree_value, i_quad));
        legend(ax, 'Location', 'northwest', 'Box', 'off');

        if ~isempty(opt.savedir)

            if ~exist(opt.savedir, 'dir')
                mkdir(opt.savedir);
            end

            file_base = fullfile(opt.savedir, sprintf('%s_%s_deg%02d', cfg.experiment, opt.metric, degree_value));

            savefig(fig, [file_base, '.fig']);
            exportgraphics(ax, [file_base, '.png'], 'Resolution', 200);
        end

        out(end+1) = rec; %#ok<AGROW>
    end

    print_plot_summary(out, cfg);
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


function i_quad = pick_quad(results_gp, results_projection, i_deg, cfg, metric)
%PICK_QUAD  Quadrature rule that was actually measured.
%
% The sweep measures a single quadrature rule by default, so the index with
% the largest number of usable points is selected.

    n_quad = size(cfg.n_quadrature{i_deg}, 1);

    i_quad = 1;
    best = -1;

    for q = 1:n_quad

        count = count_finite(cell_at(results_gp, metric, {i_deg,q}));

        if count == 0
            count = count_finite(cell_at(results_projection, metric, {i_deg,q,1}));
        end

        if count > best
            best = count;
            i_quad = q;
        end
    end
end


function count = count_finite(value)
    value = value(:);
    count = sum(isfinite(value) & value > 0);
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


function ndof = all_ndof(results_gp, results_interpolation, results_projection, i_deg, i_quad, tol_ids)
    ndof = finite_positive(cell_at(results_gp, 'ndof', {i_deg,i_quad}));

    for i_tol = tol_ids
        ndof = [ndof, finite_positive(cell_at(results_interpolation, 'ndof', {i_deg,i_quad,i_tol}))]; %#ok<AGROW>
        ndof = [ndof, finite_positive(cell_at(results_projection, 'ndof', {i_deg,i_quad,i_tol}))]; %#ok<AGROW>
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

        case 'M_memory'
            label = 'stored matrix memory';

        otherwise
            label = strrep(metric, '_', ' ');
    end
end


function print_plot_summary(out, cfg)
    fprintf('\n--- peak-memory scaling summary, %s ---\n', cfg.name);

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


function require_nargs(args, required, action)
    if numel(args) < required
        error('measure_mass_peak_memory:arguments', ...
              'Action "%s" requires at least %d additional argument(s).', action, required);
    end
end


function print_usage()
    fprintf('\nMass peak-memory workflow\n');
    fprintf('-------------------------\n');
    fprintf('Use existing measurements:\n');
    fprintf('  results = measure_mass_peak_memory(''collect'', ''flag'');\n');
    fprintf('  results = measure_mass_peak_memory(''collect_all'');\n');
    fprintf('  measure_mass_peak_memory(''plot'', ''flag'');\n\n');
    fprintf('Run the sweep:\n');
    fprintf('  measure_mass_peak_memory(''prepare'', ''flag'');\n');
    fprintf('  %% then execute the printed .sh file in a terminal\n\n');
    fprintf('Inspect a configuration:\n');
    fprintf('  cfg = measure_mass_peak_memory(''config'', ''flag'');\n\n');
end
