%% VISUALIZE_CONVERGENCE_PAPER  Generate per-solver convergence graphs for the paper.
%
% For each of the 9 solvers, produces one figure showing normdr2 vs integration
% number for ALL selected cases overlaid as thin lines, with a bold average line.
%
% Output: 9 figures total (one per solver)
%   conv_<SOLVER>.png / .tikz
%
% Usage:
%   visualize_convergence_paper()
%   visualize_convergence_paper('convergence_history.mat', 'convergence_analysis')

function visualize_convergence_paper(history_or_file, output_dir)

if nargin < 1 || isempty(history_or_file)
    history_or_file = 'convergence_history.mat';
end
if nargin < 2 || isempty(output_dir)
    output_dir = 'convergence_analysis';
end

addpath(genpath('matlab2tikz'));

if ischar(history_or_file)
    fprintf('Loading convergence history from %s...\n', history_or_file);
    load(history_or_file, 'convergence_history');
else
    convergence_history = history_or_file;
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

Tolr = 1e-8;

%% Solver display settings
solver_info = struct();
solver_info.T     = struct('label', 'T (Thompson)',           'color', [0.00 0.45 0.70]);
solver_info.MPS   = struct('label', 'MPS',                    'color', [0.85 0.33 0.10]);
solver_info.B     = struct('label', 'B (Broyden)',            'color', [0.47 0.67 0.19]);
solver_info.TR    = struct('label', 'TR (Thompson Rear)',      'color', [0.00 0.45 0.70]);
solver_info.MPSR  = struct('label', 'MPSR',                   'color', [0.85 0.33 0.10]);
solver_info.TB    = struct('label', 'TB (T+Broyden)',         'color', [0.49 0.18 0.56]);
solver_info.TRB   = struct('label', 'TRB (TR+Broyden)',       'color', [0.30 0.75 0.93]);
solver_info.MPSB  = struct('label', 'MPSB (MPS+Broyden)',    'color', [0.93 0.69 0.13]);
solver_info.MPSRB = struct('label', 'MPSRB (MPSR+Broyden)',  'color', [0.64 0.08 0.18]);
solver_info.BR    = struct('label', 'BR (Broyden Rear)',      'color', [0.47 0.67 0.19]);
solver_info.BRB   = struct('label', 'BRB (BR+Broyden)',       'color', [0.30 0.60 0.20]);

all_solvers = {'T','MPS','B','TR','MPSR','TB','TRB','MPSB','MPSRB','BR','BRB'};

n_cases = length(convergence_history);
fprintf('Generating figures for %d cases, %d solvers...\n\n', n_cases, length(all_solvers));

for s = 1:length(all_solvers)
    sname = all_solvers{s};
    info  = solver_info.(sname);

    fig = figure('Position', [100, 100, 1000, 700]);
    hold on;

    all_iters = {};
    all_errs  = {};

    for idx = 1:n_cases
        entry = convergence_history{idx};
        if ~isfield(entry, sname)
            continue;
        end
        data = entry.(sname);
        if isempty(data.normdr2_history)
            continue;
        end

        errs  = data.normdr2_history(:)';
        if isfield(data, 'integration_count_history') && length(data.integration_count_history) == length(errs)
            iters = data.integration_count_history(:)';
        else
            iters = 1:length(errs);  % fallback for old data without the vector
        end
        valid = errs > 0 & isfinite(errs) & iters > 0;
        if sum(valid) < 2
            continue;
        end

        xv = iters(valid);
        yv = errs(valid);

        semilogy(xv, yv, '-', ...
                 'Color', [info.color, 0.3], 'LineWidth', 0.5, ...
                 'HandleVisibility', 'off');

        all_iters{end+1} = xv; %#ok<AGROW>
        all_errs{end+1}  = yv; %#ok<AGROW>
    end

    %% Bold average line
    if ~isempty(all_iters)
        max_iter = max(cellfun(@max, all_iters));

        % Extend each case to max_iter using its last error value.
        % Prevents survivorship bias: converged (short) cases would otherwise
        % drop out at later iterations, pulling the mean toward slow cases only.
        for k = 1:length(all_iters)
            last_iter = all_iters{k}(end);
            if last_iter < max_iter
                ext          = (last_iter+1):max_iter;
                all_iters{k} = [all_iters{k}, ext];                    %#ok<AGROW>
                all_errs{k}  = [all_errs{k},  repmat(all_errs{k}(end), 1, numel(ext))]; %#ok<AGROW>
            end
        end

        n_cases_ext = length(all_iters);
        err_mat = NaN(n_cases_ext, max_iter);
        for k = 1:n_cases_ext
            for j = 1:length(all_iters{k})
                err_mat(k, all_iters{k}(j)) = all_errs{k}(j);
            end
        end

        med_e   = median(err_mat, 1, 'omitnan');
        has_med = ~isnan(med_e) & med_e > 0;
        semilogy(find(has_med), med_e(has_med), '-', ...
                 'Color', info.color, 'LineWidth', 2.5, ...
                 'DisplayName', [info.label ' (median)']);
    end

    %% Tolerance line
    yline(Tolr, 'k--', 'LineWidth', 1.5, ...
          'DisplayName', sprintf('Tolerance (%.0e km)', Tolr));

    %% Dynamic y-axis limits
    if ~isempty(all_errs)
        all_vals = [all_errs{:}];
        y_min = max(1e-10, min(all_vals) / 10);
        y_max = max(all_vals) * 10;
        ylim([y_min, y_max]);
    end

    hold off;
    set(gca, 'YScale', 'log', 'FontSize', 11);
    xlabel('Cumulative integration count', 'FontSize', 12);
    ylabel('Position error ||r_2|| (km)', 'FontSize', 12);
    title(sprintf('%s — Convergence (all cases)', info.label), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    grid on;

    %% Save
    base_name = sprintf('conv_%s', sname);
    png_file  = fullfile(output_dir, [base_name '.png']);
    tikz_file = fullfile(output_dir, [base_name '.tikz']);
    saveas(fig, png_file);
    if exist('matlab2tikz', 'file')
        try
            matlab2tikz(tikz_file, 'figurehandle', fig, ...
                        'width',  '\figwidth', ...
                        'height', '\figheight', ...
                        'showInfo', false, 'parseStrings', false);
        catch ME
            fprintf('  matlab2tikz failed for %s: %s\n', base_name, ME.message);
        end
    end
    fprintf('  Saved: %s\n', base_name);
    close(fig);
end

fprintf('\nDone. %d figures saved to: %s\n', length(all_solvers), output_dir);
end
