%% RUN_CONVERGENCE_ANALYSIS  Run all 9 solvers (debug mode) on selected cases.
%
% Loads convergence_cases.mat (from select_convergence_cases.m) and the
% original testcases file, then runs all 9 solvers with options.debug=true
% so that each solver returns a history struct with per-integration normdr2.
%
% TR and MPSR (backward-only intermediates) are included here — they are
% NOT part of the 75k benchmark but are needed for convergence graphs.
%
% Usage:
%   run_convergence_analysis()
%   run_convergence_analysis('convergence_cases.mat', 'testcases_75k.mat', 'convergence_history.mat')

function convergence_history = run_convergence_analysis(testcases_or_file, selected_indices_or_file, output_file, Cr, jdepoch)

if nargin < 1 || isempty(testcases_or_file)
    testcases_or_file = 'testcases_75k.mat';
end
if nargin < 2 || isempty(selected_indices_or_file)
    selected_indices_or_file = 'convergence_cases.mat';
end
if nargin < 3 || isempty(output_file)
    output_file = 'convergence_history.mat';
end
if nargin < 4 || isempty(Cr),      Cr      = 1.2;                            end
if nargin < 5 || isempty(jdepoch), jdepoch = juliandate(2026,1,1,12,0,0);    end

addpath(genpath('vallado/software/matlab'));

if ischar(testcases_or_file)
    fprintf('Loading testcases from %s...\n', testcases_or_file);
    load(testcases_or_file, 'testcases');
else
    testcases = testcases_or_file;
end

if ischar(selected_indices_or_file)
    fprintf('Loading convergence cases from %s...\n', selected_indices_or_file);
    load(selected_indices_or_file, 'selected_indices');
else
    selected_indices = selected_indices_or_file;
end

%% Configuration (match benchmark settings)
mu = 398600.5;
Re = 6378.137;
J  = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, ...
      -0.2277161016e-6, 0.5396484906e-6];
Cd     = 2.2;
muM    = 4902.8;   % Moon GM (km^3/s^2)
rhoSRP = 4.56e-3;  % Solar radiation pressure (N/km^2)

options.AbsTol  = 1e-10;
options.RelTol  = 1e-12;
options.N       = 16;
options.maxIter = 200;
options.Tolr    = 1e-8;
options.delta   = 64;
options.debug   = false;

fprintf('\nTotal unique cases to run: %d\n\n', length(selected_indices));

%% All 11 solvers (9 benchmark + 2 backward-only intermediates)
active_solvers = {'T', 'MPS', 'B', 'TR', 'MPSR', 'TB', 'TRB', 'MPSB', 'MPSRB', 'BR', 'BRB'};

%% Pre-extract selected subset so parfor doesn't broadcast full testcases array
tc_selected = testcases(selected_indices);

%% Run
n_cases = length(selected_indices);
convergence_history = cell(n_cases, 1);

parfor idx = 1:n_cases
    warning('off', 'all');
    tc = tc_selected(idx);

    r1 = tc.r1;
    r2 = tc.r2;
    v1 = [0 0 0];
    forcemodel = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, mu, Cd, tc.A, tc.m, Re, J, jdepoch, rhoSRP, Cr, tc.A, muM);

    entry = struct('test_id', tc.id, 'testcase', tc);

    for s = 1:length(active_solvers)
        sname = active_solvers{s};
        try
            switch sname
                case 'T'
                    [~,~,slv_stats,hist] = prtlambertT( ...
                        forcemodel, @odeMPCI, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);
                case 'MPS'
                    [~,~,slv_stats,hist] = prtlambertMPS( ...
                        forcemodel, @odeMPCI, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
                case 'B'
                    [~,~,slv_stats,hist] = prtlambertB( ...
                        forcemodel, @odeMPCI, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
                case 'TR'
                    [~,~,slv_stats,hist] = prtlambertTR( ...
                        forcemodel, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);
                case 'MPSR'
                    [~,~,slv_stats,hist] = prtlambertMPSR( ...
                        forcemodel, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
                case 'TB'
                    [~,~,slv_stats,hist] = prtlambertTB( ...
                        forcemodel, @odeMPCI, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);
                case 'TRB'
                    [~,~,slv_stats,hist] = prtlambertTRB( ...
                        forcemodel, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, 1, options);
                case 'MPSB'
                    [~,~,slv_stats,hist] = prtlambertMPSB( ...
                        forcemodel, @odeMPCI, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
                case 'MPSRB'
                    [~,~,slv_stats,hist] = prtlambertMPSRB( ...
                        forcemodel, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
                case 'BR'
                    [~,~,slv_stats,hist] = prtlambertBR( ...
                        forcemodel, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
                case 'BRB'
                    [~,~,slv_stats,hist] = prtlambertBRB( ...
                        forcemodel, ...
                        r1, r2, v1, tc.dm, tc.de, tc.nrev, tc.tm, options);
            end
            if isfield(hist, 'integration_number')
                int_count_hist = hist.integration_number;  % MPS-family: true cumulative ODE count
            else
                int_count_hist = hist.iterations;          % others: cumulative count (1 per iter)
            end
            entry.(sname) = struct( ...
                'converged',              hist.converged, ...
                'final_error',            hist.final_error, ...
                'normdr2_history',        hist.normdr2, ...
                'integration_number',     slv_stats.integration_number, ...
                'integration_count_history', int_count_hist, ...
                'phase',                  [], ...
                'direction_chosen',       slv_stats.direction_chosen);
            if isfield(hist, 'phase')
                entry.(sname).phase = hist.phase;
            end
        catch ME
            fprintf('  %s ERROR: %s\n', sname, ME.message);
            for k = 1:length(ME.stack)
                fprintf('    at %s > %s (line %d)\n', ME.stack(k).file, ME.stack(k).name, ME.stack(k).line);
            end
            slv_stats = struct('integration_number', 0, 'direction_chosen', []);
            entry.(sname) = struct( ...
                'converged',         false, ...
                'final_error',       NaN, ...
                'normdr2_history',   [], ...
                'integration_number', 0, ...
                'phase',             [], ...
                'direction_chosen',  []);
        end

    end

    convergence_history{idx} = entry;
end

save(output_file, 'convergence_history', 'selected_indices', '-v7.3');
fprintf('Saved convergence history to: %s\n', output_file);
end
