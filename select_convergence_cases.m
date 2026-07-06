%% SELECT_CONVERGENCE_CASES  Pick representative cases for convergence-graph analysis.
%
% Selects cases stratified by nrev × regime (25 cells):
%   N_per_cell cases per (nrev × regime) combination — 5 nrevs × 5 regimes.
%   Default N_per_cell=8 → 200 total cases.
%
% Outputs:  convergence_cases.mat  with variables:
%   selected_indices  (vector of indices into testcases)
%   nrev_counts       (1x5 vector: total count per nrev 0-4)
%   n_hard_selected   (scalar: number of hard-regime cases)
%
% Usage:
%   select_convergence_cases('benchmark_results/results_75k_consolidated.mat')
%   select_convergence_cases(results_file, output_file, N_per_cell)

function selected_indices = select_convergence_cases(testcases_or_file, output_file, N_per_cell)

if nargin < 1 || isempty(testcases_or_file)
    testcases_or_file = 'benchmark_results/results_75k_consolidated.mat';
end
if nargin < 2 || isempty(output_file)
    output_file = 'convergence_cases.mat';
end
if nargin < 3 || isempty(N_per_cell)
    N_per_cell = 8;
end

if ischar(testcases_or_file)
    fprintf('Loading testcases from %s...\n', testcases_or_file);
    load(testcases_or_file, 'testcases');
else
    testcases = testcases_or_file;
end
n_cases = length(testcases);
fprintf('Loaded %d cases.\n\n', n_cases);

regimes  = {'LEO','MEO','GEO','HEO','hard'};
nrev_set = 0:4;

rng(42);

selected_indices = [];
nrev_counts      = zeros(1, 5);
n_hard_selected  = 0;

for ni = 1:5
    nrev = nrev_set(ni);
    for ri = 1:length(regimes)
        pool = find([testcases.nrev] == nrev & ...
                    strcmp({testcases.regime}, regimes{ri}));
        n    = min(N_per_cell, length(pool));
        perm = pool(randperm(length(pool), n));
        selected_indices      = [selected_indices, perm]; %#ok<AGROW>
        nrev_counts(ni)       = nrev_counts(ni) + n;
        if strcmp(regimes{ri}, 'hard')
            n_hard_selected = n_hard_selected + n;
        end
    end
end

%% Summary
fprintf('========== SELECTION SUMMARY ==========\n');
fprintf('  N_per_cell : %d\n', N_per_cell);
fprintf('  Cases per nrev (all regimes):\n');
for nrev = 0:4
    fprintf('    nrev=%d : %d cases\n', nrev, nrev_counts(nrev+1));
end
fprintf('  Hard cases  : %d\n', n_hard_selected);
fprintf('  Total       : %d\n', length(selected_indices));
fprintf('========================================\n\n');

%% Save
save(output_file, 'selected_indices', 'nrev_counts', 'n_hard_selected');
fprintf('Saved: %s\n', output_file);
end
