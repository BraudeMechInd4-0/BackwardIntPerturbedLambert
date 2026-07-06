function visualize_direction_analysis(results_or_file, output_dir, filter_type)
%VISUALIZE_DIRECTION_ANALYSIS  Direction-selection (fwd/bwd) false-positive analysis.
%
% Loads stage-2 benchmark results and classifies each backward-chosen case for
% TRB (vs T + TB) and MPSRB (vs MPS + MPSB) as:
%   True positive  : hybrid converged AND at least one forward counterpart failed
%   False positive : hybrid failed OR all three converged (backward not needed)
%   Undecided      : hybrid converged AND exactly one forward also converged
%
% For the all-converged subset, compares hybrid runtime vs min(fwd1, fwd2) runtime.
%
% filter_type: 'all' (default), 'hard', or 'nonhard'. Output filenames get the
% corresponding suffix ('', '_hard', '_nonhard').
%
% Produces:
%   direction_analysis[_hard|_nonhard].tex   — LaTeX table of stats
%   direction_analysis[_hard|_nonhard].mat   — aggregate stats struct
%
% Usage:
%   visualize_direction_analysis()
%   visualize_direction_analysis('benchmark_results/results_75k_consolidated.mat', 'benchmark_results')
%   visualize_direction_analysis('benchmark_results/results_75k_consolidated.mat', 'benchmark_results', 'hard')
%   visualize_direction_analysis('benchmark_results/results_75k_consolidated.mat', 'benchmark_results', 'nonhard')

if nargin < 1 || isempty(results_or_file)
    results_or_file = fullfile('benchmark_results', 'results_75k_consolidated.mat');
end
if nargin < 2 || isempty(output_dir)
    output_dir = 'benchmark_results';
end
if nargin < 3 || isempty(filter_type)
    filter_type = 'all';
end
if islogical(filter_type) && filter_type
    filter_type = 'hard';
elseif islogical(filter_type)
    filter_type = 'all';
end

if ischar(results_or_file)
    fprintf('Loading %s ...\n', results_or_file);
    tmp = load(results_or_file, 'results');
    results_struct = tmp.results;
    clear tmp;
    fprintf('Done loading (%d cases).\n', length(results_struct));
    % Load once, run all three filters without reloading
    visualize_direction_analysis(results_struct, output_dir, 'all');
    visualize_direction_analysis(results_struct, output_dir, 'hard');
    visualize_direction_analysis(results_struct, output_dir, 'nonhard');
    fprintf('DONE\n');
    return;
else
    results = results_or_file;
end
if ~exist(output_dir, 'dir'),  mkdir(output_dir);  end

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

%% Analyse each hybrid solver pair
pairs = { ...
    'TRB',   'T',   'TB';   ...
    'MPSRB', 'MPS', 'MPSB'; ...
    'BRB',   'B',   ''      ...  % BRB has only one forward counterpart
};

stats_all = struct();

for pi = 1:size(pairs, 1)
    hyb  = pairs{pi, 1};
    fwd1 = pairs{pi, 2};
    fwd2 = pairs{pi, 3};

    % Cases where hybrid chose backward
    bwd_mask = false(n_results, 1);
    for k = 1:n_results
        r = results(k);
        if isfield(r, hyb) && isfield(r.(hyb), 'direction_chosen')
            bwd_mask(k) = strcmp(r.(hyb).direction_chosen, 'backward');
        end
    end
    bwd_idx = find(bwd_mask);
    n_bwd   = length(bwd_idx);

    n_true_pos   = 0;
    n_false_pos  = 0;    % TRUE FP: backward failed, >=1 forward would have converged
    n_immaterial = 0;   % both converged, or both failed — choice didn't matter
    n_both_conv  = 0;   % immaterial sub-case: hybrid + ALL forwards converged
    n_both_fail  = 0;   % immaterial sub-case: hybrid + ALL forwards failed
    n_mixed      = 0;   % immaterial sub-case: hybrid converged, only some forwards converged (TRB/MPSRB only)
    n_fwd_all    = 1 + ~isempty(fwd2);  % number of forward counterparts (1 for BRB, 2 for TRB/MPSRB)

    % All-converged subset: runtime comparison
    n_all_conv     = 0;
    n_hyb_faster   = 0;
    rt_hyb_sum     = 0;
    rt_best_fwd_sum = 0;

    for k = bwd_idx'
        r = results(k);
        hyb_conv  = isfield(r, hyb)  && r.(hyb).converged;
        fwd1_conv = isfield(r, fwd1) && r.(fwd1).converged;
        fwd2_conv = ~isempty(fwd2) && isfield(r, fwd2) && r.(fwd2).converged;

        n_fwd_conv = fwd1_conv + fwd2_conv;

        fwd_would_succeed = (n_fwd_conv >= 1);
        if hyb_conv && ~fwd_would_succeed
            n_true_pos   = n_true_pos   + 1;   % backward converged, both forwards failed
        elseif ~hyb_conv && fwd_would_succeed
            n_false_pos  = n_false_pos  + 1;   % backward failed, a forward would have worked
        else
            n_immaterial = n_immaterial + 1;   % both worked, or both failed
            if hyb_conv && n_fwd_conv == n_fwd_all
                n_both_conv = n_both_conv + 1;   % choice didn't matter — all converged
            elseif ~hyb_conv && n_fwd_conv == 0
                n_both_fail = n_both_fail + 1;   % choice didn't matter — nothing works
            else
                n_mixed = n_mixed + 1;           % hybrid converged, partial forward coverage
            end
        end

        % Speed comparison when hybrid and all available forwards converged
        all_fwd_conv = fwd1_conv && (isempty(fwd2) || fwd2_conv);
        if hyb_conv && all_fwd_conv
            rt_hyb  = r.(hyb).runtime;
            rt_fwd1 = r.(fwd1).runtime;
            if ~isempty(fwd2) && isfield(r, fwd2)
                rt_best = min(rt_fwd1, r.(fwd2).runtime);
            else
                rt_best = rt_fwd1;
            end
            if isfinite(rt_hyb) && isfinite(rt_best)
                n_all_conv      = n_all_conv + 1;
                rt_hyb_sum      = rt_hyb_sum + rt_hyb;
                rt_best_fwd_sum = rt_best_fwd_sum + rt_best;
                if rt_hyb < rt_best
                    n_hyb_faster = n_hyb_faster + 1;
                end
            end
        end
    end

    mean_rt_hyb     = NaN;
    mean_rt_best    = NaN;
    if n_all_conv > 0
        mean_rt_hyb  = rt_hyb_sum  / n_all_conv;
        mean_rt_best = rt_best_fwd_sum / n_all_conv;
    end

    % Backward rate over all filtered cases
    bwd_rate = 100 * n_bwd / n_results;

    n_decisive = n_true_pos + n_false_pos;
    if n_bwd > 0
        pct_tp         = 100 * n_true_pos   / n_bwd;
        pct_fp_raw     = 100 * n_false_pos  / n_bwd;
        pct_immaterial = 100 * n_immaterial / n_bwd;
        fp_rate        = 100 * n_false_pos  / max(n_decisive, 1);
        pct_faster     = 100 * n_hyb_faster / max(n_all_conv, 1);
    else
        pct_tp = NaN;  pct_fp_raw = NaN;  pct_immaterial = NaN;
        fp_rate = NaN;  pct_faster = NaN;
    end

    assert(n_true_pos + n_false_pos + n_immaterial == n_bwd, ...
        'Classification mismatch for %s: %d+%d+%d ~= %d', hyb, n_true_pos, n_false_pos, n_immaterial, n_bwd);
    assert(n_both_conv + n_both_fail + n_mixed == n_immaterial, ...
        'Immaterial split mismatch for %s: %d+%d+%d ~= %d', hyb, n_both_conv, n_both_fail, n_mixed, n_immaterial);

    stats_all.(hyb) = struct( ...
        'n_results',       n_results, ...
        'n_bwd',           n_bwd, ...
        'bwd_rate_pct',    bwd_rate, ...
        'n_true_pos',      n_true_pos, ...
        'n_false_pos',     n_false_pos, ...
        'n_immaterial',    n_immaterial, ...
        'n_both_conv',     n_both_conv, ...
        'n_both_fail',     n_both_fail, ...
        'n_mixed',         n_mixed, ...
        'n_decisive',      n_decisive, ...
        'pct_true_pos',    pct_tp, ...
        'pct_false_pos',   pct_fp_raw, ...
        'pct_immaterial',  pct_immaterial, ...
        'fp_rate',         fp_rate, ...
        'n_all_conv',      n_all_conv, ...
        'n_hyb_faster',    n_hyb_faster, ...
        'pct_hyb_faster',  pct_faster, ...
        'mean_rt_hyb',     mean_rt_hyb, ...
        'mean_rt_best_fwd',mean_rt_best);

    fprintf('%s: %d/%d backward (%.1f%%) | TP=%d FP=%d | fp_rate=%.1f%% (of decisive=%d) | immaterial=%d (both-conv=%d mixed=%d both-fail=%d) | faster=%.1f%%\n', ...
        hyb, n_bwd, n_results, bwd_rate, n_true_pos, n_false_pos, fp_rate, n_decisive, n_immaterial, n_both_conv, n_mixed, n_both_fail, pct_faster);
end

%% Save .mat
mat_out = fullfile(output_dir, ['direction_analysis' suffix '.mat']);
save(mat_out, 'stats_all', '-v7.3');
fprintf('Stats saved to %s\n', mat_out);

%% LaTeX table
tex_out = fullfile(output_dir, ['direction_analysis_table' suffix '.tex']);
fid = fopen(tex_out, 'w');
fprintf(fid, '%% Direction-selection analysis (filter=%s, n=%d)\n', filter_type, n_results);
fprintf(fid, '\\begin{tabular}{l r r r r r r r r r}\n\\hline\n');
fprintf(fid, 'Solver & bwd rate & TP & FP & FP rate & \\#immaterial & \\#both-conv & \\#mixed & \\#both-fail & hyb faster \\\\\n\\hline\n');
for pi = 1:size(pairs, 1)
    hyb = pairs{pi, 1};
    s   = stats_all.(hyb);
    fprintf(fid, '%s & %.1f\\%% & %d & %d & %.1f\\%% & %d & %d & %d & %d & %.1f\\%% \\\\\n', ...
        hyb, s.bwd_rate_pct, s.n_true_pos, s.n_false_pos, s.fp_rate, ...
        s.n_immaterial, s.n_both_conv, s.n_mixed, s.n_both_fail, s.pct_hyb_faster);
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fclose(fid);
fprintf('LaTeX table saved to %s\n', tex_out);
end
