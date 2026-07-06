function analyze_mps_failure_mechanism(results_or_file, output_dir, filter_type, maxIter)
%ANALYZE_MPS_FAILURE_MECHANISM  Classify failure mode of backward-chosen FP cases.
%
% For each backward-capable solver (TRB, MPSRB, BRB), identifies genuine
% false-positive cases (chose backward, hybrid failed, >=1 forward would have
% converged) and classifies the failure as:
%   integration_fail : integration_number < maxIter  (diverged before corrector finished)
%   corrector_fail   : integration_number >= maxIter (hit iteration cap)
%
% Usage:
%   analyze_mps_failure_mechanism(results, 'benchmark_results/paper_figures')
%   analyze_mps_failure_mechanism(results, output_dir, 'hard', 200)
%   analyze_mps_failure_mechanism('benchmark_results/results_75k_consolidated.mat', ...)

if nargin < 1 || isempty(results_or_file)
    results_or_file = fullfile('benchmark_results', 'results_75k_consolidated.mat');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('benchmark_results', 'paper_figures');
end
if nargin < 3 || isempty(filter_type)
    filter_type = 'all';
end
if nargin < 4 || isempty(maxIter)
    maxIter = 200;
end

if ischar(results_or_file)
    fprintf('Loading %s ...\n', results_or_file);
    tmp = load(results_or_file, 'results');
    results = tmp.results;
    clear tmp;
    fprintf('Loaded (%d cases).\n', length(results));
else
    results = results_or_file;
end

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

switch filter_type
    case 'hard'
        is_hard = arrayfun(@(r) isfield(r,'testcase') && isfield(r.testcase,'is_hard') && ...
            r.testcase.is_hard, results);
        results = results(is_hard);
        suffix  = '_hard';
    case 'nonhard'
        is_hard = arrayfun(@(r) isfield(r,'testcase') && isfield(r.testcase,'is_hard') && ...
            r.testcase.is_hard, results);
        results = results(~is_hard);
        suffix  = '_nonhard';
    otherwise
        suffix = '';
end

n_results = length(results);

pairs = {
    'TRB',   'T',   'TB';
    'MPSRB', 'MPS', 'MPSB';
    'BRB',   'B',   ''      % BRB has only one forward counterpart
};

stats_all = struct();

for pi = 1:size(pairs, 1)
    hyb  = pairs{pi, 1};
    fwd1 = pairs{pi, 2};
    fwd2 = pairs{pi, 3};

    n_integ_fail    = 0;
    n_corrector_fail = 0;
    n_fp = 0;

    for k = 1:n_results
        r = results(k);
        if ~isfield(r, hyb) || ~isfield(r.(hyb), 'direction_chosen')
            continue;
        end
        chose_bwd = strcmp(r.(hyb).direction_chosen, 'backward');
        if ~chose_bwd, continue; end

        hyb_conv  = r.(hyb).converged;
        fwd1_conv = isfield(r, fwd1) && r.(fwd1).converged;
        fwd2_conv = ~isempty(fwd2) && isfield(r, fwd2) && r.(fwd2).converged;
        fwd_would_succeed = fwd1_conv || fwd2_conv;

        % Genuine FP: chose backward, hybrid failed, forward would have worked
        if ~hyb_conv && fwd_would_succeed
            n_fp = n_fp + 1;
            n_iter = r.(hyb).integration_number;
            if n_iter < maxIter
                n_integ_fail = n_integ_fail + 1;
            else
                n_corrector_fail = n_corrector_fail + 1;
            end
        end
    end

    assert(n_integ_fail + n_corrector_fail == n_fp, ...
        'Count mismatch for %s: %d+%d ~= %d', hyb, n_integ_fail, n_corrector_fail, n_fp);

    pct_integ = 100 * n_integ_fail    / max(n_fp, 1);
    pct_corr  = 100 * n_corrector_fail / max(n_fp, 1);

    stats_all.(hyb) = struct( ...
        'n_fp',             n_fp, ...
        'n_integ_fail',     n_integ_fail, ...
        'n_corrector_fail', n_corrector_fail, ...
        'pct_integ_fail',   pct_integ, ...
        'pct_corrector_fail', pct_corr, ...
        'maxIter_used',     maxIter);

    fprintf('%s FP cases (filter=%s): N=%d\n', hyb, filter_type, n_fp);
    fprintf('  integration_fail : %d (%.1f%%)\n', n_integ_fail, pct_integ);
    fprintf('  corrector_fail   : %d (%.1f%%)\n', n_corrector_fail, pct_corr);
    fprintf('  (classification threshold: integration_number < %d)\n\n', maxIter);
end

%% Save .mat
mat_out = fullfile(output_dir, ['mps_failure_mechanism' suffix '.mat']);
save(mat_out, 'stats_all', '-v7.3');
fprintf('Stats saved to %s\n', mat_out);

%% LaTeX table
tex_out = fullfile(output_dir, ['mps_failure_mechanism_table' suffix '.tex']);
fid = fopen(tex_out, 'w');
fprintf(fid, '%% FP failure mechanism analysis (filter=%s, maxIter=%d, n=%d)\n', filter_type, maxIter, n_results);
fprintf(fid, '\\begin{tabular}{l r r r r r}\n\\hline\n');
fprintf(fid, 'Solver & \\#FP & Integ-fail & Integ-fail\\%% & Corr-fail & Corr-fail\\%% \\\\\n\\hline\n');
for pi = 1:size(pairs, 1)
    hyb = pairs{pi, 1};
    s   = stats_all.(hyb);
    fprintf(fid, '%s & %d & %d & %.1f\\%% & %d & %.1f\\%% \\\\\n', ...
        hyb, s.n_fp, s.n_integ_fail, s.pct_integ_fail, s.n_corrector_fail, s.pct_corrector_fail);
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fclose(fid);
fprintf('LaTeX table saved to %s\n\n', tex_out);
end
