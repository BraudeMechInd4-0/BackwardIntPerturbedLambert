function visualize_switch_sensitivity(sensitivity_or_file, output_dir, filter_type)
%VISUALIZE_SWITCH_SENSITIVITY  LaTeX tables and TikZ heatmaps for T/N_switch sensitivity.
%
% Loads raw per-case data saved by run_switch_sensitivity and aggregates here.
% filter_type: 'all' (default), 'hard', or 'nonhard'. Output filenames get the
% corresponding suffix ('', '_hard', '_nonhard').
%
% Produces:
%   switch_sensitivity_table[_hard|_nonhard].tex      — convergence count + mean runtime tables
%   switch_sensitivity_conv[_hard|_nonhard].tikz      — heatmap: n_conv per grid point
%   switch_sensitivity_runtime[_hard|_nonhard].tikz   — heatmap: mean runtime per grid point
%
% Usage:
%   visualize_switch_sensitivity()
%   visualize_switch_sensitivity('benchmark_results/switch_sensitivity.mat', 'benchmark_results')
%   visualize_switch_sensitivity('benchmark_results/switch_sensitivity.mat', 'benchmark_results', 'hard')
%   visualize_switch_sensitivity('benchmark_results/switch_sensitivity.mat', 'benchmark_results', 'nonhard')

if nargin < 1 || isempty(sensitivity_or_file)
    sensitivity_or_file = fullfile('benchmark_results', 'switch_sensitivity.mat');
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

if ischar(sensitivity_or_file)
    s = load(sensitivity_or_file);
    conv_per_case         = s.conv_per_case;
    runtime_per_case      = s.runtime_per_case;
    integrations_per_case = s.integrations_per_case;
    is_hard_sampled       = s.is_hard_sampled;
    n_sampled             = s.n_sampled;
    T_switch_vals         = s.T_switch_vals;
    N_switch_vals         = s.N_switch_vals;
    solvers_tested        = s.solvers_tested;
    if isfield(s, 'regime_sampled')
        regime_sampled = s.regime_sampled;
    else
        regime_sampled = {};
    end
else
    sn = sensitivity_or_file;
    conv_per_case         = sn.conv_per_case;
    runtime_per_case      = sn.runtime_per_case;
    integrations_per_case = sn.integrations_per_case;
    is_hard_sampled       = sn.is_hard_sampled;
    n_sampled             = sn.n_sampled;
    T_switch_vals         = sn.T_switch_vals;
    N_switch_vals         = sn.N_switch_vals;
    solvers_tested        = sn.solvers_tested;
    if isfield(sn, 'regime_sampled')
        regime_sampled = sn.regime_sampled;
    else
        regime_sampled = {};
    end
end

if ~exist(output_dir, 'dir'),  mkdir(output_dir);  end

%% Apply filter mask
switch filter_type
    case 'hard'
        mask   = is_hard_sampled;
        suffix = '_hard';
        n_used = sum(mask);
    case 'nonhard'
        mask   = ~is_hard_sampled;
        suffix = '_nonhard';
        n_used = sum(mask);
    otherwise
        mask   = true(n_sampled, 1);
        suffix = '';
        n_used = n_sampled;
end

conv_use   = conv_per_case(mask, :, :, :);
rt_use     = runtime_per_case(mask, :, :, :);
integ_use  = integrations_per_case(mask, :, :, :);

nT = length(T_switch_vals);
nN = length(N_switch_vals);
nS = length(solvers_tested);

%% Aggregate
n_conv            = zeros(nT, nN, nS);
mean_runtime      = NaN(nT, nN, nS);
mean_integrations = NaN(nT, nN, nS);

for ti = 1:nT
    for ni = 1:nN
        for si = 1:nS
            c = conv_use(:, ti, ni, si);
            r = rt_use(:, ti, ni, si);
            g = integ_use(:, ti, ni, si);
            n_conv(ti, ni, si) = sum(c);
            r_ok = r(isfinite(r));
            g_ok = g(isfinite(g));
            if ~isempty(r_ok),  mean_runtime(ti, ni, si) = mean(r_ok);      end
            if ~isempty(g_ok),  mean_integrations(ti, ni, si) = mean(g_ok); end
        end
    end
end

fmt_T = @(v) sprintf('$10^{%d}$', round(log10(v)));

%% ---- LaTeX tables ----
tex_file = fullfile(output_dir, ['switch_sensitivity_table' suffix '.tex']);
fid = fopen(tex_file, 'w');

% Convergence count table
fprintf(fid, '%% Convergence count (n_used=%d, filter=%s)\n', n_used, filter_type);
for si = 1:nS
    fprintf(fid, '\\subsection*{%s — Converged cases}\n', solvers_tested{si});
    fprintf(fid, '\\begin{tabular}{l');
    for ni = 1:nN,  fprintf(fid, ' r');  end
    fprintf(fid, '}\n\\hline\n');
    fprintf(fid, '$T_{\\mathrm{sw}}$ \\textbackslash{} $N_{\\mathrm{sw}}$');
    for ni = 1:nN,  fprintf(fid, ' & %d', N_switch_vals(ni));  end
    fprintf(fid, ' \\\\\n\\hline\n');
    for ti = 1:nT
        fprintf(fid, '%s', fmt_T(T_switch_vals(ti)));
        for ni = 1:nN,  fprintf(fid, ' & %d', n_conv(ti, ni, si));  end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n\n');
end

% Mean runtime table
fprintf(fid, '%% Mean runtime (seconds, converged cases only)\n');
for si = 1:nS
    fprintf(fid, '\\subsection*{%s — Mean runtime (s)}\n', solvers_tested{si});
    fprintf(fid, '\\begin{tabular}{l');
    for ni = 1:nN,  fprintf(fid, ' r');  end
    fprintf(fid, '}\n\\hline\n');
    fprintf(fid, '$T_{\\mathrm{sw}}$ \\textbackslash{} $N_{\\mathrm{sw}}$');
    for ni = 1:nN,  fprintf(fid, ' & %d', N_switch_vals(ni));  end
    fprintf(fid, ' \\\\\n\\hline\n');
    for ti = 1:nT
        fprintf(fid, '%s', fmt_T(T_switch_vals(ti)));
        for ni = 1:nN
            if isnan(mean_runtime(ti, ni, si))
                fprintf(fid, ' & ---');
            else
                fprintf(fid, ' & %.4f', mean_runtime(ti, ni, si));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n\n');
end

% Mean integrations table
fprintf(fid, '%% Mean integration count (converged cases only)\n');
for si = 1:nS
    fprintf(fid, '\\subsection*{%s — Mean integrations}\n', solvers_tested{si});
    fprintf(fid, '\\begin{tabular}{l');
    for ni = 1:nN,  fprintf(fid, ' r');  end
    fprintf(fid, '}\n\\hline\n');
    fprintf(fid, '$T_{\\mathrm{sw}}$ \\textbackslash{} $N_{\\mathrm{sw}}$');
    for ni = 1:nN,  fprintf(fid, ' & %d', N_switch_vals(ni));  end
    fprintf(fid, ' \\\\\n\\hline\n');
    for ti = 1:nT
        fprintf(fid, '%s', fmt_T(T_switch_vals(ti)));
        for ni = 1:nN
            if isnan(mean_integrations(ti, ni, si))
                fprintf(fid, ' & ---');
            else
                fprintf(fid, ' & %.1f', mean_integrations(ti, ni, si));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n\n');
end

% Per-regime convergence summary (averaged over the full T/N grid)
% Skipped for hard-only (single regime, trivial)
if ~isempty(regime_sampled) && ~strcmp(filter_type, 'hard')
    all_regimes = {'LEO','MEO','GEO','HEO','hard'};
    fprintf(fid, '%% Per-regime convergence rate (averaged over T/N grid)\n');
    fprintf(fid, '\\subsection*{Regime breakdown — convergence rate (\\%%)}\n');
    fprintf(fid, '\\begin{tabular}{l');
    for si = 1:nS,  fprintf(fid, ' r');  end
    fprintf(fid, '}\n\\hline\n');
    fprintf(fid, 'Regime');
    for si = 1:nS,  fprintf(fid, ' & %s', solvers_tested{si});  end
    fprintf(fid, ' \\\\\n\\hline\n');
    for ri = 1:length(all_regimes)
        rg    = all_regimes{ri};
        rmask = mask & strcmp(regime_sampled, rg);
        n_rg  = sum(rmask);
        if n_rg == 0,  continue;  end
        fprintf(fid, '%s', rg);
        for si = 1:nS
            c_rg = conv_per_case(rmask, :, :, si);
            rate = 100 * sum(c_rg(:)) / (n_rg * nT * nN);
            fprintf(fid, ' & %.1f', rate);
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n\n');
end

fclose(fid);
fprintf('LaTeX table saved to %s\n', tex_file);

%% ---- TikZ heatmap helper ----
    function write_heatmap(tikz_path, data, data_fmt, title_prefix, cbar_label)
        all_vals = data(isfinite(data(:)));
        if isempty(all_vals)
            cmin = 0;  cmax = 1;
        else
            cmin = min(all_vals);  cmax = max(all_vals);
        end
        if cmin == cmax,  cmax = cmin + 1;  end

        fh = fopen(tikz_path, 'w');
        fprintf(fh, '\\begin{tikzpicture}\n');

        for sii = 1:nS
            y_shift = (nS - sii) * (nT + 3);
            fprintf(fh, '\\begin{scope}[yshift=%dcm]\n', y_shift);
            fprintf(fh, '\\begin{axis}[\n');
            fprintf(fh, '  title={%s: %s},\n', title_prefix, solvers_tested{sii});
            fprintf(fh, '  xlabel={$N_{\\mathrm{sw}}$},\n');
            fprintf(fh, '  ylabel={$T_{\\mathrm{sw}}$},\n');
            fprintf(fh, '  xtick={1,...,%d},\n', nN);
            fprintf(fh, '  xticklabels={');
            for nii = 1:nN
                if nii > 1,  fprintf(fh, ',');  end
                fprintf(fh, '%d', N_switch_vals(nii));
            end
            fprintf(fh, '},\n');
            fprintf(fh, '  ytick={1,...,%d},\n', nT);
            fprintf(fh, '  yticklabels={');
            for tii = 1:nT
                if tii > 1,  fprintf(fh, ',');  end
                fprintf(fh, '$10^{%d}$', round(log10(T_switch_vals(tii))));
            end
            fprintf(fh, '},\n');
            fprintf(fh, '  colormap/hot,\n');
            fprintf(fh, '  colorbar,\n');
            fprintf(fh, '  colorbar style={ylabel={%s}},\n', cbar_label);
            fprintf(fh, '  point meta min=%.6g, point meta max=%.6g,\n', cmin, cmax);
            fprintf(fh, '  matrix plot*,\n');
            fprintf(fh, '  mesh/rows=%d, mesh/cols=%d,\n', nT, nN);
            fprintf(fh, '  nodes near coords,\n');
            fprintf(fh, '  nodes near coords align={center},\n');
            fprintf(fh, '  every node near coord/.append style={font=\\scriptsize,white},\n');
            fprintf(fh, ']\n');
            fprintf(fh, '\\addplot[matrix plot*,point meta=explicit] coordinates {\n');
            for tii = 1:nT
                for nii = 1:nN
                    val = data(tii, nii, sii);
                    if isnan(val),  val = cmin;  end
                    fprintf(fh, ['  (%d,%d) [' data_fmt ']\n'], nii, tii, val);
                end
            end
            fprintf(fh, '};\n\\end{axis}\n\\end{scope}\n\n');
        end

        fprintf(fh, '\\end{tikzpicture}\n');
        fclose(fh);
        fprintf('TikZ heatmap saved to %s\n', tikz_path);
    end

%% Convergence count heatmap
conv_tikz = fullfile(output_dir, ['switch_sensitivity_conv' suffix '.tikz']);
write_heatmap(conv_tikz, double(n_conv), '%d', 'Converged cases', 'count');
save_heatmap_png(fullfile(output_dir, ['switch_sensitivity_conv' suffix '.png']), ...
    double(n_conv), 'Converged Cases', N_switch_vals, T_switch_vals, solvers_tested);

%% Mean runtime heatmap
rt_tikz = fullfile(output_dir, ['switch_sensitivity_runtime' suffix '.tikz']);
write_heatmap(rt_tikz, mean_runtime, '%.4f', 'Mean runtime (s)', 's');
save_heatmap_png(fullfile(output_dir, ['switch_sensitivity_runtime' suffix '.png']), ...
    mean_runtime, 'Mean Runtime (s)', N_switch_vals, T_switch_vals, solvers_tested);

end

%% ---- PNG heatmap helper ----
function save_heatmap_png(png_path, data, cbar_label, N_vals, T_vals, solver_names)
    nS = length(solver_names);
    nT = length(T_vals);
    nN = length(N_vals);
    figure('Visible', 'off');
    for sii = 1:nS
        subplot(1, nS, sii);
        imagesc(data(:, :, sii));
        colorbar;
        colormap hot;
        set(gca, 'XTick', 1:nN, 'XTickLabel', arrayfun(@num2str, N_vals, 'UniformOutput', false));
        set(gca, 'YTick', 1:nT, 'YTickLabel', ...
            arrayfun(@(v) sprintf('1e%d', round(log10(v))), T_vals, 'UniformOutput', false));
        xlabel('N_{sw}');  ylabel('T_{sw}');
        title(sprintf('%s: %s', cbar_label, solver_names{sii}));
        % annotate cells
        for ti = 1:nT
            for ni = 1:nN
                val = data(ti, ni, sii);
                if isfinite(val)
                    text(ni, ti, sprintf('%.3g', val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'w');
                end
            end
        end
    end
    saveas(gcf, png_path);
    close(gcf);
    fprintf('PNG heatmap saved to %s\n', png_path);
end
