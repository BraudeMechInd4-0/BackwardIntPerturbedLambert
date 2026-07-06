function sensitivity = run_switch_sensitivity(testcases_file, output_dir, n_sample)
%RUN_SWITCH_SENSITIVITY  Sensitivity analysis on T_switch and N_switch for TRB and MPSRB.
%
% Runs TRB and MPSRB on a stratified random sample over a grid of (T_switch, N_switch)
% values. Saves raw per-case results — all aggregation is done in
% visualize_switch_sensitivity.
%
% Saved variables:
%   conv_per_case         (n_sampled x nT x nN x nS, logical)
%   runtime_per_case      (n_sampled x nT x nN x nS, double, NaN if not converged)
%   integrations_per_case (n_sampled x nT x nN x nS, double, NaN if not converged)
%   is_hard_sampled       (n_sampled x 1, logical)
%   sampled_indices, n_sampled, T_switch_vals, N_switch_vals, solvers_tested
%
% Usage:
%   run_switch_sensitivity()
%   run_switch_sensitivity('benchmark_results/testcases_75k.mat', 'benchmark_results', 8000)

if nargin < 1 || isempty(testcases_file)
    testcases_file = fullfile('benchmark_results', 'testcases_75k.mat');
end
if nargin < 2 || isempty(output_dir)
    output_dir = 'benchmark_results';
end
if nargin < 3 || isempty(n_sample)
    n_sample = 8000;
end

fprintf('========================================\n');
fprintf('Switch Sensitivity Analysis\n');
fprintf('  Solvers  : TRB, MPSRB, BRB\n');
fprintf('  n_sample : %d\n', n_sample);
fprintf('========================================\n\n');

if ischar(testcases_file)
    fprintf('Loading testcases from %s...\n', testcases_file);
    load(testcases_file, 'testcases');
else
    testcases = testcases_file;
end
n_total = length(testcases);
fprintf('Loaded %d test cases.\n\n', n_total);

%% Stratified sample: floor(n_sample/25) per nrev×regime cell
rng(42);
regimes    = {'LEO','MEO','GEO','HEO','hard'};
n_cells    = 5 * 5;
n_per_cell = floor(n_sample / n_cells);
sampled_indices = [];
for nrev = 0:4
    for ri = 1:5
        pool = find([testcases.nrev] == nrev & ...
                    strcmp({testcases.regime}, regimes{ri}));
        n    = min(n_per_cell, length(pool));
        perm = pool(randperm(length(pool), n));
        sampled_indices = [sampled_indices, perm]; %#ok<AGROW>
    end
end
n_sampled       = length(sampled_indices);
is_hard_sampled = logical([testcases(sampled_indices).is_hard])';
regime_sampled  = {testcases(sampled_indices).regime}';

%% Pre-extract sampled subset so parfor doesn't broadcast the full testcases array
tc_sampled = testcases(sampled_indices);
fprintf('Sampled %d cases (%d per nrev×regime cell, %d hard).\n\n', ...
    n_sampled, n_per_cell, sum(is_hard_sampled));

%% Constants — match run_benchmark_batch.m exactly
mu      = 398600.5;
Re      = 6378.137;
J       = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, ...
           -0.2277161016e-6, 0.5396484906e-6];
Cd      = 2.2;
Cr      = 1.2;
muM     = 4902.8;
rhoSRP  = 4.56e-3;
jdepoch = juliandate(2026, 1, 1, 12, 0, 0);

make_fm = @(A, m) @(t,x) orbit_eq_J6_drag_SRP_moon( ...
    t, x, mu, Cd, A, m, Re, J, jdepoch, rhoSRP, Cr, A, muM);

base_opts.AbsTol  = 1e-10;
base_opts.RelTol  = 1e-12;
base_opts.N       = 16;
base_opts.maxIter = 200;
base_opts.Tolr    = 1e-8;
base_opts.delta   = 32;
base_opts.debug   = false;

%% Parameter grids
T_switch_vals  = [1e-6, 1e-5, 1e-4];
N_switch_vals  = [5, 10, 15];
solvers_tested = {'TRB', 'MPSRB', 'BRB'};

nT = length(T_switch_vals);
nN = length(N_switch_vals);
nS = length(solvers_tested);

%% Per-case result arrays — initialise fresh, then overwrite from checkpoint if found
if ~exist(output_dir, 'dir'),  mkdir(output_dir);  end
solver_tag      = strjoin(sort(solvers_tested), '_');
checkpoint_file = fullfile(output_dir, ['switch_sensitivity_checkpoint_' solver_tag '.mat']);

conv_per_case         = false(n_sampled, nT, nN, nS);
runtime_per_case      = NaN(n_sampled,   nT, nN, nS);
integrations_per_case = NaN(n_sampled,   nT, nN, nS);
completed             = false(nT, nN);

if exist(checkpoint_file, 'file')
    fprintf('Checkpoint found — resuming from %s\n', checkpoint_file);
    ck = load(checkpoint_file, ...
        'conv_per_case','runtime_per_case','integrations_per_case','completed');
    conv_per_case         = ck.conv_per_case;
    runtime_per_case      = ck.runtime_per_case;
    integrations_per_case = ck.integrations_per_case;
    completed             = ck.completed;
    fprintf('  %d / %d grid cells already done.\n\n', sum(completed(:)), nT*nN);
end

warning('off', 'all');

for ti = 1:nT
    for ni = 1:nN
        if completed(ti, ni)
            fprintf('SKIP  T_switch=%.0e  N_switch=%d  (checkpoint)\n', ...
                T_switch_vals(ti), N_switch_vals(ni));
            continue;
        end

        opts          = base_opts;
        opts.T_switch = T_switch_vals(ti);
        opts.N_switch = N_switch_vals(ni);

        fprintf('T_switch=%.0e  N_switch=%d\n', T_switch_vals(ti), N_switch_vals(ni));

        for si = 1:nS
            solver = solvers_tested{si};

            conv_k         = false(n_sampled, 1);
            runtime_k      = NaN(n_sampled, 1);
            integrations_k = NaN(n_sampled, 1);

            parfor k = 1:n_sampled
                tc_k = tc_sampled(k);
                if ~isempty(tc_k.dm)
                    fm_k = make_fm(tc_k.A, tc_k.m);
                    try
                        switch solver
                            case 'TRB'
                                [v1t, ~, st] = prtlambertTRB(fm_k, tc_k.r1, tc_k.r2, [0 0 0], ...
                                    tc_k.dm, tc_k.de, tc_k.nrev, tc_k.tm, 1, opts);
                            case 'MPSRB'
                                [v1t, ~, st] = prtlambertMPSRB(fm_k, tc_k.r1, tc_k.r2, [0 0 0], ...
                                    tc_k.dm, tc_k.de, tc_k.nrev, tc_k.tm, opts);
                            case 'BRB'
                                [v1t, ~, st] = prtlambertBRB(fm_k, tc_k.r1, tc_k.r2, [0 0 0], ...
                                    tc_k.dm, tc_k.de, tc_k.nrev, tc_k.tm, opts);
                        end
                        if all(isreal(v1t)) && all(isfinite(v1t)) && ~any(isnan(v1t))
                            conv_k(k)         = true;
                            runtime_k(k)      = st.runtime_phase1 + st.runtime_phase2;
                            integrations_k(k) = st.integration_number;
                        end
                    catch
                    end
                end
            end

            conv_per_case(:, ti, ni, si)         = conv_k;
            runtime_per_case(:, ti, ni, si)      = runtime_k;
            integrations_per_case(:, ti, ni, si) = integrations_k;

            fprintf('  %s: %d / %d converged,  mean_runtime=%.4fs\n', ...
                solver, sum(conv_k), n_sampled, mean(runtime_k(isfinite(runtime_k))));
        end

        %% Checkpoint after each completed grid cell
        completed(ti, ni) = true;
        save(checkpoint_file, ...
            'conv_per_case', 'runtime_per_case', 'integrations_per_case', 'completed', ...
            'is_hard_sampled', 'regime_sampled', 'sampled_indices', 'n_sampled', ...
            'T_switch_vals', 'N_switch_vals', 'solvers_tested', '-v7.3');
        fprintf('  [checkpoint saved: %d/%d cells done]\n', sum(completed(:)), nT*nN);
    end
end

warning('on', 'all');

%% Assemble return struct and remove checkpoint
sensitivity.conv_per_case         = conv_per_case;
sensitivity.runtime_per_case      = runtime_per_case;
sensitivity.integrations_per_case = integrations_per_case;
sensitivity.is_hard_sampled       = is_hard_sampled;
sensitivity.regime_sampled        = regime_sampled;
sensitivity.sampled_indices       = sampled_indices;
sensitivity.n_sampled             = n_sampled;
sensitivity.T_switch_vals         = T_switch_vals;
sensitivity.N_switch_vals         = N_switch_vals;
sensitivity.solvers_tested        = solvers_tested;
if exist(checkpoint_file, 'file'),  delete(checkpoint_file);  end
end
