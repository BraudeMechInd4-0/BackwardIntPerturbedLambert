function compute_failure_transitions(results_or_file, output_dir, filter_type, maxIter)
%COMPUTE_FAILURE_TRANSITIONS  Full-path state table for solver family chains.
%
% For every case that failed on at least one solver in a family chain,
% shows how that same case performed on ALL solvers in the chain.
%
%   Families:
%     Thompson : T, TB, TRB   (3 states each -> 27 paths, show 26 non-all-success)
%     MPS      : MPS, MPSB, MPSRB
%     Broyden  : B, BRB        (9 paths, show 8 non-all-success)
%
%   States:
%     success (1) : converged == true
%     explode (2) : ~converged && integration_number <  maxIter
%     stall   (3) : ~converged && integration_number >= maxIter
%
% Outputs (to output_dir):
%   failure_transitions[suffix].csv           -- one row per path (all 27/9)
%   failure_transitions[suffix].mat           -- path_counts structs
%   failure_transitions_[family][suffix].tex  -- table of non-all-success paths

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
    fprintf('Done loading (%d cases).\n', length(results));
else
    results = results_or_file;
end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

fprintf('\n[compute_failure_transitions] fields: converged, integration_number (maxIter=%d)\n', maxIter);

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

n_cases     = length(results);
STATE_NAMES = {'success', 'explode', 'stall'};

families = {
    'thompson', {'T',   'TB',   'TRB'  };
    'mps',      {'MPS', 'MPSB', 'MPSRB'};
    'broyden',  {'B',   'BRB'          }
};

all_data = struct();
csv_rows = {'family,filter,path,count,pct'};

for fi = 1:size(families, 1)
    fam_name = families{fi, 1};
    solvers  = families{fi, 2};
    n_s      = length(solvers);

    % Pre-compute state for each solver (0 = field missing)
    states = zeros(n_cases, n_s);
    for si = 1:n_s
        sol = solvers{si};
        for k = 1:n_cases
            r = results(k);
            if isfield(r, sol) && isfield(r.(sol), 'converged') && isfield(r.(sol), 'integration_number')
                states(k, si) = classify_state(r.(sol).converged, r.(sol).integration_number, maxIter);
            end
        end
    end

    % Only cases where every solver has a valid state
    valid_mask   = all(states > 0, 2);
    n_valid      = sum(valid_mask);
    states_valid = states(valid_mask, :);

    fprintf('\n=== %s (filter=%s, n=%d, valid=%d) ===\n', upper(fam_name), filter_type, n_cases, n_valid);
    if n_valid < n_cases
        fprintf('  NOTE: %d cases skipped (missing fields)\n', n_cases - n_valid);
    end

    % Build path count array
    if n_s == 3
        path_counts = zeros(3, 3, 3);
        for k = 1:n_valid
            s = states_valid(k, :);
            path_counts(s(1), s(2), s(3)) = path_counts(s(1), s(2), s(3)) + 1;
        end
    else  % n_s == 2
        path_counts = zeros(3, 3);
        for k = 1:n_valid
            s = states_valid(k, :);
            path_counts(s(1), s(2)) = path_counts(s(1), s(2)) + 1;
        end
    end
    assert(sum(path_counts(:)) == n_valid, 'path_counts sum mismatch for %s', fam_name);

    % Print: header
    sol_header = sprintf('%-9s', solvers{:});
    fprintf('\n  %-s  Count    %%\n', sol_header);
    fprintf('  %s\n', repmat('-', 1, 9*n_s + 18));

    % Enumerate all paths, sorted by count descending, skip all-success
    if n_s == 3
        paths = zeros(27, 3);
        counts = zeros(27, 1);
        idx = 0;
        for s1 = 1:3
            for s2 = 1:3
                for s3 = 1:3
                    idx = idx + 1;
                    paths(idx,:) = [s1 s2 s3];
                    counts(idx)  = path_counts(s1, s2, s3);
                end
            end
        end
    else
        paths = zeros(9, 2);
        counts = zeros(9, 1);
        idx = 0;
        for s1 = 1:3
            for s2 = 1:3
                idx = idx + 1;
                paths(idx,:) = [s1 s2];
                counts(idx)  = path_counts(s1, s2);
            end
        end
    end

    % Sort descending by count
    [counts_sorted, ord] = sort(counts, 'descend');
    paths_sorted = paths(ord, :);

    % All-success path index: all states == 1
    is_all_success = all(paths_sorted == 1, 2);

    % Print all-success line first for context, then failures
    for pass = 1:2
        if pass == 1
            mask = is_all_success;
            if ~any(mask), continue; end
            fprintf('  [all succeeded]\n');
        else
            mask = ~is_all_success;
            if ~any(mask), continue; end
            fprintf('  [at least one failure]\n');
        end
        sel_paths  = paths_sorted(mask, :);
        sel_counts = counts_sorted(mask);
        for row = 1:size(sel_paths, 1)
            if sel_counts(row) == 0, continue; end
            path_str = '';
            for si = 1:n_s
                path_str = [path_str sprintf('%-9s', STATE_NAMES{sel_paths(row,si)})]; %#ok<AGROW>
            end
            pct = 100 * sel_counts(row) / n_valid;
            fprintf('  %s  %6d   %5.1f%%\n', path_str, sel_counts(row), pct);

            % CSV
            path_lbl = strjoin(arrayfun(@(x) STATE_NAMES{x}, sel_paths(row,:), 'UniformOutput', false), '-');
            csv_rows{end+1} = sprintf('%s,%s,%s,%d,%.2f', fam_name, filter_type, path_lbl, sel_counts(row), pct); %#ok<AGROW>
        end
    end

    all_data.(fam_name) = struct('path_counts', path_counts, 'n_valid', n_valid);

    % TEX: non-all-success paths only, sorted by count descending
    tex_out = fullfile(output_dir, sprintf('failure_transitions_%s%s.tex', fam_name, suffix));
    fid = fopen(tex_out, 'w');
    col_fmt = repmat('l ', 1, n_s);
    fprintf(fid, '%% Failure paths -- %s (filter=%s, n_valid=%d, maxIter=%d)\n', fam_name, filter_type, n_valid, maxIter);
    fprintf(fid, '\\begin{tabular}{%srr}\n\\hline\n', col_fmt);
    % Header
    hdr = '';
    for si = 1:n_s
        hdr = [hdr solvers{si} ' & ']; %#ok<AGROW>
    end
    fprintf(fid, '%sCount & \\%% \\\\\n\\hline\n', hdr);
    % Rows: non-all-success with count > 0
    for row = 1:size(paths_sorted, 1)
        if is_all_success(row) || counts_sorted(row) == 0, continue; end
        row_str = '';
        for si = 1:n_s
            row_str = [row_str STATE_NAMES{paths_sorted(row,si)} ' & ']; %#ok<AGROW>
        end
        pct = 100 * counts_sorted(row) / n_valid;
        fprintf(fid, '%s%d & %.1f\\%% \\\\\n', row_str, counts_sorted(row), pct);
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n');
    fclose(fid);
    fprintf('  LaTeX saved: %s\n', tex_out);
end

% MAT
mat_out = fullfile(output_dir, sprintf('failure_transitions%s.mat', suffix));
save(mat_out, 'all_data', 'filter_type', 'n_cases', 'maxIter', '-v7.3');
fprintf('\nMAT saved: %s\n', mat_out);

% CSV
csv_out = fullfile(output_dir, sprintf('failure_transitions%s.csv', suffix));
fid = fopen(csv_out, 'w');
for i = 1:length(csv_rows)
    fprintf(fid, '%s\n', csv_rows{i});
end
fclose(fid);
fprintf('CSV saved: %s\n', csv_out);

end

%% -------------------------------------------------------------------------
function s = classify_state(converged, integration_number, maxIter)
if converged
    s = 1;
elseif integration_number < maxIter
    s = 2;
else
    s = 3;
end
end
