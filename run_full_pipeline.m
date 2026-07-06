%% RUN_FULL_PIPELINE  Complete 75k benchmark + convergence analysis + paper figures.
%
% Stages (each can be skipped by setting the corresponding flag to false):
%
%   Stage 1 — Generate 75k test cases  (skip if testcases_75k.mat already exists)
%   Stage 2 — Run 7-solver benchmark   (skip to load a previous consolidated result)
%   Stage 3 — Summary statistics table
%   Stage 4 — Select convergence cases
%   Stage 5 — Run convergence analysis (9 solvers, debug mode, selected cases)
%   Stage 6 — Switch sensitivity analysis (T_switch / N_switch grid)
%   Stage 7 — Generate paper figures (benchmark + convergence + sensitivity + direction)
%
% Outputs (all saved under output_dir):
%   results_75k_consolidated.mat   — benchmark results
%   convergence_cases.mat          — selected case indices
%   convergence_history.mat        — per-integration error history
%   switch_sensitivity.mat         — raw per-case sensitivity grid data
%   direction_analysis*.mat        — direction-selection false-positive stats
%   paper_figures/                 — tikz + png for all paper figures
%   convergence_analysis/          — tikz + png convergence graphs
%
% Vallado submodule changes (vallado/software/matlab/):
%   lambhodograph.m — replaced constastro script call with direct constants
%                     (mu1=398600.4415, twopi=2*pi) and renamed mu→mu1 throughout.
%                     Reason: constastro is a workspace-injection script; in parfor
%                     workers it failed silently and mu collided with MATLAB's built-in
%                     mu() function, producing v1t=[0;0;0] (wrong shape and value).
%   lambertb.m      — constmath/constastro calls were already commented out;
%                     no functional change, kept in sync with lambhodograph fix.

clc; close all;
warning('off', 'all')
addpath(genpath('vallado/software/matlab'))
addpath(genpath('matlab2tikz'))

diary_file = sprintf('pipeline_%s.log', datestr(now, 'yyyymmdd_HHMMSS'));
diary(diary_file);
diary on;
fprintf('Pipeline started: %s\n\n', datestr(now));

%% =========================================================================
%%  CONFIGURATION — edit these flags and paths as needed
%% =========================================================================

% Stages to run
run_generation         = true;  % Stage 1: generate 75k test cases
run_benchmark          = true;  % Stage 2: run 8-solver benchmark
run_analysis           = true;  % Stage 3: summary statistics table
run_select_cases       = true;  % Stage 4: select convergence cases
run_conv_analysis      = true;  % Stage 5: convergence analysis (11 solvers, debug mode)
run_sensitivity        = true;  % Stage 6: switch sensitivity (T_switch/N_switch grid)
run_figures            = true;  % Stage 7: generate paper figures (master switch)
run_figures_benchmark  = true;  % Stage 7a: benchmark bar/box figures
run_figures_conv       = true;  % Stage 7b: convergence graphs
run_figures_sensitivity= true;  % Stage 7c: switch sensitivity figures
run_figures_direction  = true;  % Stage 7d: direction-selection analysis

% Paths
output_dir       = 'benchmark_results';
testcases_file   = fullfile(output_dir, 'testcases_75k.mat');
results_file     = fullfile(output_dir, 'results_75k_consolidated.mat');
cases_file       = fullfile(output_dir, 'convergence_cases.mat');
history_file     = fullfile(output_dir, 'convergence_history.mat');
sensitivity_file = fullfile(output_dir, 'switch_sensitivity.mat');
figures_dir      = fullfile(output_dir, 'paper_figures');
conv_figures_dir = fullfile(output_dir, 'convergence_analysis');

% Switch sensitivity settings
switch_n_sample   = 8000;   % cases to sample for T/N_switch sensitivity grid

% Benchmark settings
Tolr              = 1e-8;
maxIter           = 200;
AbsTol            = 1e-10;
RelTol            = 1e-12;
N_cheb            = 16;
delta_cheb        = 32;
T_switch          = Tolr * 1000;
N_switch          = floor(maxIter * 0.05);
Cd                = 2.2;
Cr                = 1.2;
jdepoch           = juliandate(2026, 1, 1, 12, 0, 0);
checkpoint_interval = 5000;

% Convergence case selection
N_per_combo = 8;   % cases per (nrev × regime) cell → 25 cells × 8 = 200 total

% Ensure output directory exists before any stage writes to it
if ~exist(output_dir, 'dir'),  mkdir(output_dir);  end

%% =========================================================================
%%  Stage 1 — Generate test cases
%% =========================================================================
if run_generation
    fprintf('========================================\n');
    fprintf('Stage 1: Generating 75k test cases...\n');
    fprintf('========================================\n\n');
    generate_benchmark_testcases(75000, testcases_file,AbsTol,RelTol,N_cheb,delta_cheb,maxIter);
    fprintf('Done.\n\n');
else
    fprintf('Stage 1: SKIPPED (using %s)\n\n', testcases_file);
end

%% =========================================================================
%%  Stage 2 — Run 75k benchmark
%% =========================================================================
if run_benchmark
    fprintf('========================================\n');
    fprintf('Stage 2: Running 75k benchmark...\n');
    fprintf('  Solvers: T, MPS, B, TB, TRB, MPSB, MPSRB, BRB\n');
    fprintf('  Expected runtime: ~8-10 h on 8-core machine\n');
    fprintf('========================================\n\n');
    tic;
    [results, testcases] = run_benchmark_batch( ...
        'testcases_file',      testcases_file, ...
        'checkpoint_interval', checkpoint_interval, ...
        'output_dir',          output_dir, ...
        'debug',               false, ...
        'Tolr',                Tolr, ...
        'maxIter',             maxIter, ...
        'AbsTol',              AbsTol, ...
        'RelTol',              RelTol, ...
        'N',                   N_cheb, ...
        'delta',               delta_cheb, ...
        'T_switch',            T_switch, ...
        'N_switch',            N_switch, ...
        'Cd',                  Cd, ...
        'Cr',                  Cr, ...
        'jdepoch',             jdepoch);
    elapsed = toc;
    fprintf('\nBenchmark completed in %.2f hours.\n\n', elapsed/3600);
    save(results_file, 'results', 'testcases', '-v7.3');
    fprintf('Saved to %s\n\n', results_file);
else
    if exist('results', 'var') && exist('testcases', 'var')
        fprintf('Stage 2: SKIPPED — results/testcases already in workspace\n\n');
    else
        fprintf('Stage 2: SKIPPED — loading %s\n\n', results_file);
        load(results_file, 'results', 'testcases');
    end
end

%% =========================================================================
%%  Stage 3 — Summary statistics table
%% =========================================================================
if run_analysis
    fprintf('========================================\n');
    fprintf('Stage 3: Generating summary statistics...\n');
    fprintf('========================================\n\n');
    analyze_benchmark_results(results, testcases, output_dir, 'all',     maxIter);
    analyze_benchmark_results(results, testcases, output_dir, 'hard',    maxIter);
    analyze_benchmark_results(results, testcases, output_dir, 'nonhard', maxIter);
    fprintf('\n  Failure mechanism check...\n');
    analyze_mps_failure_mechanism(results, figures_dir, 'all',     maxIter);
    analyze_mps_failure_mechanism(results, figures_dir, 'hard',    maxIter);
    analyze_mps_failure_mechanism(results, figures_dir, 'nonhard', maxIter);
    fprintf('\n  Failure transition matrices...\n');
    compute_failure_transitions(results, figures_dir, 'all',     maxIter);
    compute_failure_transitions(results, figures_dir, 'hard',    maxIter);
    compute_failure_transitions(results, figures_dir, 'nonhard', maxIter);
    fprintf('Done.\n\n');
else
    fprintf('Stage 3: SKIPPED\n\n');
end

%% =========================================================================
%%  Stage 4 — Select convergence cases
%% =========================================================================
if run_select_cases
    fprintf('========================================\n');
    fprintf('Stage 4: Selecting convergence cases...\n');
    fprintf('========================================\n\n');
    selected_indices = select_convergence_cases(testcases, cases_file, N_per_combo);
    fprintf('\n');
else
    fprintf('Stage 4: SKIPPED — loading %s\n\n', cases_file);
    tmp = load(cases_file, 'selected_indices');
    selected_indices = tmp.selected_indices;
end

%% =========================================================================
%%  Stage 5 — Convergence analysis (debug runs)
%% =========================================================================
if run_conv_analysis
    fprintf('========================================\n');
    fprintf('Stage 5: Running convergence analysis...\n');
    fprintf('  (all 11 solvers, debug mode, selected cases)\n');
    fprintf('========================================\n\n');
    convergence_history = run_convergence_analysis(testcases, selected_indices, history_file, Cr, jdepoch);
    fprintf('\n');
else
    fprintf('Stage 5: SKIPPED — loading %s\n\n', history_file);
    if isfile(history_file)
        tmp = load(history_file, 'convergence_history');
        convergence_history = tmp.convergence_history;
    else
        fprintf('  (convergence_history.mat not found — Stage 7b will be skipped)\n\n');
        convergence_history = [];
    end
end

%% =========================================================================
%%  Stage 6 — Switch sensitivity analysis
%% =========================================================================
if run_sensitivity
    fprintf('========================================\n');
    fprintf('Stage 6: Running switch sensitivity analysis...\n');
    fprintf('  (TRB, MPSRB, BRB — %d sampled cases, T/N_switch grid)\n', switch_n_sample);
    fprintf('========================================\n\n');
    sensitivity = run_switch_sensitivity(testcases, output_dir, switch_n_sample);
    save(sensitivity_file, '-struct', 'sensitivity', '-v7.3');
    fprintf('Sensitivity saved to %s\n\n', sensitivity_file);
else
    fprintf('Stage 6: SKIPPED — loading %s\n\n', sensitivity_file);
    sensitivity = load(sensitivity_file);
end

%% =========================================================================
%%  Stage 7 — Paper figures
%% =========================================================================
if run_figures
    fprintf('========================================\n');
    fprintf('Stage 7: Generating paper figures...\n');
    fprintf('========================================\n\n');

    if run_figures_benchmark
        fprintf('  7a. Benchmark figures...\n');
        visualize_benchmark_paper(results, testcases, figures_dir, Tolr);
        visualize_benchmark_paper(results, testcases, figures_dir, Tolr, 'hard');
        visualize_benchmark_paper(results, testcases, figures_dir, Tolr, 'nonhard');
    else
        fprintf('  7a. SKIPPED\n');
    end

    if run_figures_conv
        if ~isempty(convergence_history)
            fprintf('\n  7b. Convergence graphs...\n');
            visualize_convergence_paper(convergence_history, conv_figures_dir);
        else
            fprintf('\n  7b. Skipped (convergence_history not available)\n');
        end
    else
        fprintf('  7b. SKIPPED\n');
    end

    if run_figures_sensitivity
        fprintf('\n  7c. Switch sensitivity figures...\n');
        visualize_switch_sensitivity(sensitivity, figures_dir);
        visualize_switch_sensitivity(sensitivity, figures_dir, 'hard');
        visualize_switch_sensitivity(sensitivity, figures_dir, 'nonhard');
    else
        fprintf('  7c. SKIPPED\n');
    end

    if run_figures_direction
        fprintf('\n  7d. Direction-selection analysis...\n');
        visualize_direction_analysis(results, figures_dir);
        visualize_direction_analysis(results, figures_dir, 'hard');
        visualize_direction_analysis(results, figures_dir, 'nonhard');
    else
        fprintf('  7d. SKIPPED\n');
    end
else
    fprintf('Stage 7: SKIPPED\n\n');
end

%% =========================================================================
%%  Done
%% =========================================================================
fprintf('\n========================================\n');
fprintf('PIPELINE COMPLETE\n');
fprintf('  Benchmark figures : %s\n', figures_dir);
fprintf('  Convergence graphs: %s\n', conv_figures_dir);
fprintf('========================================\n\n');
fprintf('Pipeline finished: %s\n', datestr(now));
diary off;
