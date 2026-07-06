function [results, testcases] = run_benchmark_batch(varargin)
% RUN_BENCHMARK_BATCH Execute Lambert solver benchmark on all test cases
%
% Usage:
%   % Run and save to file only
%   run_benchmark_batch('checkpoint_interval', 500, 'output_dir', 'benchmark_results')
%
%   % Run and return results for in-memory processing
%   [results, testcases] = run_benchmark_batch('output_dir', 'test_results')
%
%   % Resume from checkpoint
%   run_benchmark_batch('resume', true, 'output_dir', 'benchmark_results')
%
% Parameters:
%   resume              - Resume from latest checkpoint (default: false)
%   checkpoint_interval - Save checkpoint every N cases (default: 500)
%   output_dir          - Output directory for results (default: 'benchmark_results')
%   testcases_file      - Input test cases file (default: 'testcases_75k.mat')
%   debug               - Print detailed output on solver failures (default: false)
%   Tolr                - Position convergence tolerance in km (default: 1e-8)
%   AbsTol              - ODE absolute tolerance (default: 1e-10)
%   RelTol              - ODE relative tolerance (default: 1e-12)
%   N                   - Chebyshev nodes (default: 16)
%   maxIter             - Max solver iterations (default: 200)
%   delta               - Step-size divisor for orbital period (default: 32)
%   T_switch            - Residual threshold to switch to Broyden (default: Tolr*1000)
%   N_switch            - Max Phase-1 iterations before switch (default: floor(maxIter*0.05))
%   Cd                  - Drag coefficient (default: 2.2)
%
% Returns:
%   results   - Array of result structs (one per test case)
%   testcases - Array of test case structs
%
% Runs 7 Lambert solvers (T, MPS, B, TB, TRB, MPSB, MPSRB) on each test case with:
%   - Parallel execution using parfor
%   - Checkpoint recovery for long runs
%   - Error handling (record failures, continue)
%
% Estimated runtime: ~8-10 hours on 8-core system (75k cases)

    % Parse input arguments
    p = inputParser;
    p.addParameter('resume', false, @islogical);
    p.addParameter('checkpoint_interval', 500, @isnumeric);
    p.addParameter('output_dir', 'benchmark_results', @ischar);
    p.addParameter('testcases_file', 'testcases_75k.mat', @ischar);
    p.addParameter('debug', false, @islogical);
    p.addParameter('Cd',      2.2,                             @isnumeric);
    p.addParameter('Cr',      1.2,                             @isnumeric);
    p.addParameter('jdepoch', juliandate(2026,1,1,12,0,0),     @isnumeric);
    p.addParameter('Tolr', 1e-8, @isnumeric);
    p.addParameter('AbsTol', 1e-10, @isnumeric);
    p.addParameter('RelTol', 1e-12, @isnumeric);
    p.addParameter('N', 16, @isnumeric);
    p.addParameter('maxIter', 200, @isnumeric);
    p.addParameter('delta', 32, @isnumeric);
    p.addParameter('T_switch', [], @(x) isempty(x) || isnumeric(x));
    p.addParameter('N_switch', [], @(x) isempty(x) || isnumeric(x));
    p.parse(varargin{:});

    % Ensure vallado functions are on path for parfor workers
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'vallado', 'software', 'matlab')));

    % Suppress warnings during benchmark (restore at end)
    warning_state = warning('query', 'all');
    warning('off', 'all');

    % Create output directory if needed
    if ~exist(p.Results.output_dir, 'dir')
        mkdir(p.Results.output_dir);
    end

    Cd      = p.Results.Cd;
    Cr      = p.Results.Cr;
    jdepoch = p.Results.jdepoch;

    % Load test cases
    fprintf('Loading test cases from %s...\n', p.Results.testcases_file);
    load(p.Results.testcases_file, 'testcases');
    n_cases = length(testcases);
    fprintf('Loaded %d test cases\n\n', n_cases);

    active_solvers = {'T', 'MPS', 'B', 'TB', 'TRB', 'MPSB', 'MPSRB', 'BRB'};

    % Initialize or resume
    if p.Results.resume
        [results, start_idx] = resume_from_checkpoint(p.Results.output_dir, n_cases);
    else
        results = initialize_results(n_cases, active_solvers);
        start_idx = 1;
    end

    % Constants
    mu = 398600.5;
    Re = 6378.137;
    r_GEO = 42164;  % km (GEO radius from Earth center)
    r_max_allowed = 1.5 * r_GEO;  % km
    J = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, ...
         -0.2277161016e-6, 0.5396484906e-6];
    muM    = 4902.8;   % Moon GM (km^3/s^2)
    rhoSRP = 4.56e-3;  % Solar radiation pressure (N/km^2)
    %Cd = 0;%2.2;
    

    % ODE options
    options.AbsTol  = p.Results.AbsTol;
    options.RelTol  = p.Results.RelTol;
    options.N       = p.Results.N;
    options.maxIter = p.Results.maxIter;
    options.Tolr    = p.Results.Tolr;
    options.delta   = p.Results.delta;
    options.debug   = p.Results.debug;
    if ~isempty(p.Results.T_switch),  options.T_switch = p.Results.T_switch;  end
    if ~isempty(p.Results.N_switch),  options.N_switch = p.Results.N_switch;  end

    % Force pool restart so workers pick up the correct lambertb/lambhodograph
    delete(gcp('nocreate'));

    % Parallel execution with checkpoint batches
    fprintf('Starting benchmark from case %d to %d...\n', start_idx, n_cases);
    fprintf('Checkpoint interval: %d cases\n', p.Results.checkpoint_interval);
    fprintf('========================================\n\n');

    batch_size = p.Results.checkpoint_interval;
    total_time_start = tic;

    for batch_start = start_idx:batch_size:n_cases
        batch_end = min(batch_start + batch_size - 1, n_cases);
        batch_indices = batch_start:batch_end;

        fprintf('Processing batch %d-%d (%d cases)...\n', ...
                batch_start, batch_end, length(batch_indices));
        batch_tic = tic;

        % Parallel loop over batch (use regular for loop in debug mode)
        results_batch = cell(length(batch_indices), 1);

        if p.Results.debug
            % Serial execution in debug mode (allows keyboard)
            for i = 1:length(batch_indices)
                case_idx = batch_indices(i);
                tc = testcases(case_idx);

                % if case_idx == 13
                %     keyboard
                % end


                % Run all 4 solvers on this test case
                results_batch{i} = run_all_solvers(tc, mu, Re, J, Cd, jdepoch, muM, Cr, rhoSRP, options, p.Results.debug, active_solvers);
            end
        else
            % Parallel execution in normal mode
            parfor i = 1:length(batch_indices)
                warning('off', 'all');
                case_idx = batch_indices(i);
                tc = testcases(case_idx);

                % Run all solvers on this test case
                results_batch{i} = run_all_solvers(tc, mu, Re, J, Cd, jdepoch, muM, Cr, rhoSRP, options, false, active_solvers);
            end
        end

        % Consolidate results
        for i = 1:length(batch_indices)
            results(batch_indices(i)) = results_batch{i};
        end

        batch_time = toc(batch_tic);
        avg_time_per_case = batch_time / length(batch_indices);

        fprintf('  Batch completed in %.1f sec (%.2f sec/case)\n', ...
                batch_time, avg_time_per_case);

        % Save checkpoint
        checkpoint_file = sprintf('%s/checkpoint_batch_%06d.mat', ...
                                  p.Results.output_dir, batch_start);
        save(checkpoint_file, 'results', 'batch_start', 'batch_end', '-v7.3');
        fprintf('  Checkpoint saved: %s\n', checkpoint_file);

        % Progress estimate
        cases_completed = batch_end;
        cases_remaining = n_cases - cases_completed;
        time_elapsed = toc(total_time_start);
        time_per_case = time_elapsed / cases_completed;
        time_remaining = time_per_case * cases_remaining;

        fprintf('  Progress: %d/%d cases (%.1f%%)\n', ...
                cases_completed, n_cases, 100*cases_completed/n_cases);
        fprintf('  Estimated time remaining: %.1f hours\n\n', time_remaining/3600);
    end

    total_time = toc(total_time_start);

    % Save final consolidated results
    fprintf('========================================\n');
    fprintf('Benchmark complete!\n');
    fprintf('Total time: %.2f hours\n', total_time/3600);
    fprintf('Average time per case: %.2f sec\n\n', total_time/n_cases);

    % Print quick summary
    print_quick_summary(results, active_solvers);

    % Restore warning state
    warning(warning_state);
end

function result = run_all_solvers(tc, mu, Re, J, Cd, jdepoch, muM, Cr, rhoSRP, options, debug_mode, solvers)
% Run specified perturbed Lambert solvers on a single test case

    % Constants for orbital period calculation
    r_GEO = 42164;  % km (GEO radius from Earth center)
    r_max_allowed = 1.5 * r_GEO;  % km
    options.debug = debug_mode;  % Sync debug flag from parameter into options struct

    % Initialize result with test case data
    result = struct('test_id', tc.id, 'testcase', tc);

    % Positions and velocities from test case (transpose to row vectors for solvers)
    r1 = tc.r1;  % Row vector
    r2 = tc.r2;  % Row vector
    v1 = [0 0 0];  % Zero initial guess

    % Force model
    forcemodel = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, mu, Cd, tc.A, tc.m, Re, J, jdepoch, rhoSRP, Cr, tc.A, muM);

    % Run each solver
    for s = 1:length(solvers)
        solver_name = solvers{s};

        % Call appropriate perturbed Lambert solver
        switch solver_name
            case 'T'
                [v1t, v2t, slv_stats] = prtlambertT(...
                    forcemodel, @odeMPCI, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);

            case 'MPS'
                [v1t, v2t, slv_stats] = prtlambertMPS(...
                    forcemodel, @odeMPCI, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);

            case 'B'
                [v1t, v2t, slv_stats] = prtlambertB(...
                    forcemodel, @odeMPCI, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);

            case 'TB'
                [v1t, v2t, slv_stats] = prtlambertTB(...
                    forcemodel, @odeMPCI, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);

            case 'TRB'
                [v1t, v2t, slv_stats] = prtlambertTRB(...
                    forcemodel, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);

            case 'MPSB'
                [v1t, v2t, slv_stats] = prtlambertMPSB(...
                    forcemodel, @odeMPCI, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);

            case 'MPSRB'
                [v1t, v2t, slv_stats] = prtlambertMPSRB(...
                    forcemodel, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);

            case 'BR'
                [v1t, v2t, slv_stats] = prtlambertBR(...
                    forcemodel, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);

            case 'BRB'
                [v1t, v2t, slv_stats] = prtlambertBRB(...
                    forcemodel, ...
                    r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
        end

            % Check convergence
            converged_check = all(isreal(v1t)) && all(~isnan(v1t)) && all(isfinite(v1t));

            % Debug: log suspicious cases
            if debug_mode && ~converged_check
                fprintf('DEBUG: %s returned v1t=[%.6f, %.6f, %.6f], isreal=%d, hasNaN=%d, isfinite=[%d %d %d]\n', ...
                    solver_name, v1t, all(isreal(v1t)), any(isnan(v1t)), isfinite(v1t));
            end

            if converged_check
                % Compute accuracy: propagate with achieved v1 (row vectors)
                % Set Sec parameter based on orbital period
                oe_verify = kep_elements(r1, v1t, mu);
                a_verify = oe_verify(1);
                if a_verify <= 0 || ~isreal(a_verify) || isnan(a_verify) || a_verify > 2*r_max_allowed
                    % Hyperbolic or very large orbit - use max altitude approximation
                    T_orb_verify = 2*pi*sqrt((2*r_max_allowed)^3/mu);
                else
                    T_orb_verify = sqrt(a_verify^3/mu) * 2*pi;
                end
                options_verify = options;
                options_verify.Sec = T_orb_verify / options.delta;

                [~, x_verify] = odeMPCI(forcemodel, [0 tc.tm], [r1, v1t], options_verify);

                % Check if verification propagation succeeded
                if all(isfinite(x_verify(end,:))) && ~any(isnan(x_verify(end,:)))
                    r2_achieved = x_verify(end, 1:3);  % Keep as row vector
                    dr2 = norm(r2 - r2_achieved);

                    % Achieved v2 from propagation
                    v2_achieved = x_verify(end, 4:6);  % Keep as row vector

                    result.(solver_name) = struct(...
                        'v1', v1t', ...
                        'v2', v2t', ...
                        'v2_achieved', v2_achieved', ...
                        'dr2', dr2, ...
                        'runtime', slv_stats.runtime_preamble + slv_stats.runtime_phase1 + slv_stats.runtime_phase2, ...
                        'integration_number',   slv_stats.integration_number, ...
                        'integration_preamble', slv_stats.integration_preamble, ...
                        'integration_phase1',   slv_stats.integration_phase1, ...
                        'integration_phase2',   slv_stats.integration_phase2, ...
                        'runtime_preamble',     slv_stats.runtime_preamble, ...
                        'runtime_phase1',       slv_stats.runtime_phase1, ...
                        'runtime_phase2',       slv_stats.runtime_phase2, ...
                        'direction_chosen',     slv_stats.direction_chosen, ...
                        'hitearthvar',          slv_stats.hitearthvar, ...
                        'hitrad',               slv_stats.hitrad, ...
                        'converged', true);
                else
                    % Solver converged but verification failed (e.g., satellite decay)
                    result.(solver_name) = struct(...
                        'v1', v1t', ...
                        'v2', v2t', ...
                        'v2_achieved', [NaN; NaN; NaN], ...
                        'dr2', NaN, ...
                        'runtime', slv_stats.runtime_preamble + slv_stats.runtime_phase1 + slv_stats.runtime_phase2, ...
                        'integration_number',   slv_stats.integration_number, ...
                        'integration_preamble', slv_stats.integration_preamble, ...
                        'integration_phase1',   slv_stats.integration_phase1, ...
                        'integration_phase2',   slv_stats.integration_phase2, ...
                        'runtime_preamble',     slv_stats.runtime_preamble, ...
                        'runtime_phase1',       slv_stats.runtime_phase1, ...
                        'runtime_phase2',       slv_stats.runtime_phase2, ...
                        'direction_chosen',     slv_stats.direction_chosen, ...
                        'hitearthvar',          slv_stats.hitearthvar, ...
                        'hitrad',               slv_stats.hitrad, ...
                        'converged', false);

                    if debug_mode
                        fprintf('\n*** VERIFICATION FAILED ***\n');
                        fprintf('Test case ID: %d\n', tc.id);
                        fprintf('Solver: %s returned valid v1t but propagation failed\n', solver_name);
                        fprintf('x_verify(end,:) = [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f]\n\n', x_verify(end,:));
                    end
                end
            else
                % Convergence failure
                result.(solver_name) = struct(...
                    'v1', [NaN; NaN; NaN], ...
                    'v2', [NaN; NaN; NaN], ...
                    'v2_achieved', [NaN; NaN; NaN], ...
                    'dr2', NaN, ...
                    'runtime', slv_stats.runtime_preamble + slv_stats.runtime_phase1 + slv_stats.runtime_phase2, ...
                    'integration_number',   slv_stats.integration_number, ...
                    'integration_preamble', slv_stats.integration_preamble, ...
                    'integration_phase1',   slv_stats.integration_phase1, ...
                    'integration_phase2',   slv_stats.integration_phase2, ...
                    'runtime_preamble',     slv_stats.runtime_preamble, ...
                    'runtime_phase1',       slv_stats.runtime_phase1, ...
                    'runtime_phase2',       slv_stats.runtime_phase2, ...
                    'direction_chosen',     slv_stats.direction_chosen, ...
                    'hitearthvar',          slv_stats.hitearthvar, ...
                    'hitrad',               slv_stats.hitrad, ...
                    'converged', false);

                if debug_mode
                    fprintf('\n*** CONVERGENCE FAILURE DETECTED ***\n');
                    fprintf('Test case ID: %d\n', tc.id);
                    fprintf('Solver: %s\n', solver_name);
                    fprintf('Returned v1t: [%.6f, %.6f, %.6f]\n', v1t);
                    fprintf('Integrations: %d\n\n', slv_stats.integration_number);
                end
            end
    end
end

function results = initialize_results(n_cases, solvers)
% Initialize empty results array

    fprintf('Initializing results structure for %d cases...\n', n_cases);

    % Create template
    template = struct('test_id', 0, 'testcase', struct());

    for s = 1:length(solvers)
        template.(solvers{s}) = struct(...
            'v1', [NaN; NaN; NaN], ...
            'v2', [NaN; NaN; NaN], ...
            'v2_achieved', [NaN; NaN; NaN], ...
            'dr2', NaN, ...
            'runtime', NaN, ...
            'integration_number',   NaN, ...
            'integration_preamble', NaN, ...
            'integration_phase1',   NaN, ...
            'integration_phase2',   NaN, ...
            'runtime_preamble',     NaN, ...
            'runtime_phase1',       NaN, ...
            'runtime_phase2',       NaN, ...
            'direction_chosen',     [], ...
            'hitearthvar',          NaN, ...
            'hitrad',               NaN, ...
            'converged', false);
    end

    results = repmat(template, n_cases, 1);
end

function [results, start_idx] = resume_from_checkpoint(output_dir, n_cases)
% Resume from the latest checkpoint

    % Find all checkpoint files
    checkpoint_files = dir(fullfile(output_dir, 'checkpoint_batch_*.mat'));

    if isempty(checkpoint_files)
        error('No checkpoint files found in %s', output_dir);
    end

    % Sort by date and load latest
    [~, latest_idx] = max([checkpoint_files.datenum]);
    latest_file = fullfile(output_dir, checkpoint_files(latest_idx).name);

    fprintf('Resuming from checkpoint: %s\n', latest_file);
    checkpoint_data = load(latest_file, 'results', 'batch_end');

    results = checkpoint_data.results;
    start_idx = checkpoint_data.batch_end + 1;

    fprintf('Resuming from case %d\n\n', start_idx);

    % Verify results array size matches
    if length(results) ~= n_cases
        error('Checkpoint has %d cases but testcases file has %d', ...
              length(results), n_cases);
    end
end

function print_quick_summary(results, solvers)
% Print quick summary of results

    fprintf('\n========================================\n');
    fprintf('QUICK SUMMARY\n');
    fprintf('========================================\n\n');
    n_cases = length(results);

    fprintf('%-10s %10s %15s %15s\n', 'Solver', 'Converged', 'Mean Runtime', 'Median dr2');
    fprintf('%s\n', repmat('-', 1, 55));

    for s = 1:length(solvers)
        solver = solvers{s};

        % Convergence rate - extract using arrayfun
        converged = arrayfun(@(x) x.(solver).converged, results);
        n_converged = sum(converged);

        % Runtime statistics - extract using arrayfun
        runtimes = arrayfun(@(x) x.(solver).runtime, results);
        runtimes_valid = runtimes(~isnan(runtimes));
        if ~isempty(runtimes_valid)
            mean_runtime = mean(runtimes_valid);
        else
            mean_runtime = NaN;
        end

        % Accuracy statistics - extract using arrayfun
        dr2_vals = arrayfun(@(x) x.(solver).dr2, results);
        dr2_valid = dr2_vals(converged);
        if ~isempty(dr2_valid)
            median_dr2 = median(dr2_valid);
        else
            median_dr2 = NaN;
        end

        fprintf('%-10s %6d/%6d %12.3f s %12.2e km\n', ...
                solver, n_converged, n_cases, mean_runtime, median_dr2);
    end

    fprintf('\n');
end
