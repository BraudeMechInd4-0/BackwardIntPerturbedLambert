function visualize_benchmark_paper(results_or_file, varargin)
% VISUALIZE_BENCHMARK_PAPER  Generate comparison figures for 8 solvers.
%
% Usage:
%   visualize_benchmark_paper(results_file, output_dir, Tolr)
%   visualize_benchmark_paper(results_file, output_dir, Tolr, filter_type)
%   visualize_benchmark_paper(results, testcases, output_dir, Tolr)
%   visualize_benchmark_paper(results, testcases, output_dir, Tolr, filter_type)
%
%   filter_type  : 'all' (default) | 'hard' | 'nonhard'

    addpath(genpath('matlab2tikz'));

    if ischar(results_or_file)
        % file mode: (results_file, output_dir, Tolr, filter_type)
        fprintf('Loading results from %s...\n', results_or_file);
        load(results_or_file, 'results', 'testcases');
        if numel(varargin) >= 1, output_dir  = varargin{1}; else output_dir  = '.';   end
        if numel(varargin) >= 2, Tolr        = varargin{2}; else Tolr        = 1e-8;  end
        if numel(varargin) >= 3, filter_type = varargin{3}; else filter_type = 'all'; end
    else
        % in-memory mode: (results, testcases, output_dir, Tolr, filter_type)
        results   = results_or_file;
        testcases = varargin{1};
        if numel(varargin) >= 2, output_dir  = varargin{2}; else output_dir  = '.';   end
        if numel(varargin) >= 3, Tolr        = varargin{3}; else Tolr        = 1e-8;  end
        if numel(varargin) >= 4, filter_type = varargin{4}; else filter_type = 'all'; end
    end

    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    %% Apply filter
    switch filter_type
        case 'hard'
            fmask  = arrayfun(@(tc) tc.is_hard, testcases);
            suffix = '_hard';
        case 'nonhard'
            fmask  = ~arrayfun(@(tc) tc.is_hard, testcases);
            suffix = '_nonhard';
        otherwise
            fmask  = true(length(testcases), 1);
            suffix = '';
    end
    results   = results(fmask);
    testcases = testcases(fmask);

    fprintf('Output directory: %s  [filter=%s, n=%d]\n', output_dir, filter_type, length(results));
    fprintf('Generating visualizations...\n\n');

    solvers = {'T', 'MPS', 'B', 'TB', 'TRB', 'MPSB', 'MPSRB', 'BRB'};
    colors  = lines(8);
    Re      = 6378.137;  % km

    %% Figure 1: Runtime Distribution
    figure;
    runtime_data = [];
    labels = {};
    for s = 1:length(solvers)
        runtimes = arrayfun(@(x) x.(solvers{s}).runtime, results);
        runtimes = runtimes(~isnan(runtimes));
        runtime_data = [runtime_data; runtimes(:)]; %#ok<AGROW>
        labels = [labels; repmat(solvers(s), length(runtimes), 1)]; %#ok<AGROW>
    end
    boxplot(runtime_data, labels, 'Colors', colors);
    ylabel('Runtime (seconds)', 'FontSize', 12);
    title('Runtime Distribution by Solver', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    saveas(gcf, fullfile(output_dir, ['runtime_distribution' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['runtime_distribution' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: runtime_distribution%s\n', suffix);
    close(gcf);

    %% Figure 2: Convergence by Regime
    figure;
    regimes    = {'LEO','MEO','GEO','HEO','hard'};
    conv_rates = zeros(length(solvers), length(regimes));
    n_per_regime = zeros(1, length(regimes));
    for ri = 1:length(regimes)
        rg_mask = strcmp({testcases.regime}, regimes{ri});
        n_per_regime(ri) = sum(rg_mask);
        if n_per_regime(ri) == 0, continue; end
        for s = 1:length(solvers)
            conv = arrayfun(@(x) x.(solvers{s}).converged, results(rg_mask));
            conv_rates(s, ri) = 100 * sum(conv) / n_per_regime(ri);
        end
    end
    % Only plot regimes that have cases
    present = n_per_regime > 0;
    bar_handle = bar(find(present), conv_rates(:, present)', 'grouped');
    set(gca, 'XTick', find(present), 'XTickLabel', regimes(present), 'FontSize', 11);
    ylabel('Convergence Rate (%)', 'FontSize', 12);
    legend(solvers, 'Location', 'best', 'FontSize', 10);
    title('Convergence Rate by Regime', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    ylim([0 105]);
    for s = 1:length(solvers), bar_handle(s).FaceColor = colors(s,:); end
    saveas(gcf, fullfile(output_dir, ['convergence_by_regime' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['convergence_by_regime' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: convergence_by_regime%s\n', suffix);
    close(gcf);

    %% Figures 3-9: Accuracy vs Transfer Time (one per solver)
    for s = 1:length(solvers)
        figure;
        converged = arrayfun(@(x) x.(solvers{s}).converged, results);
        dr2_vals  = arrayfun(@(x) x.(solvers{s}).dr2,       results);
        tm_vals   = [testcases.tm] / 3600;  % hours
        dr2_plot  = dr2_vals(converged);
        tm_plot   = tm_vals(converged);
        valid     = dr2_plot > 0;
        dr2_plot  = dr2_plot(valid);
        tm_plot   = tm_plot(valid);
        if ~isempty(dr2_plot)
            scatter(tm_plot, log10(dr2_plot), 20, colors(s,:), 'filled', 'MarkerFaceAlpha', 0.3);
        end
        xlabel('Transfer Time (hours)', 'FontSize', 11);
        ylabel('log_{10}(Position Error) [km]', 'FontSize', 11);
        title(sprintf('%s: Accuracy vs Transfer Time', solvers{s}), 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca, 'FontSize', 10);
        hold on;
        yline(log10(0.001), 'g--', '1 m',  'LineWidth', 1.5, 'FontSize', 9);
        yline(log10(0.1),   'y--', '100 m', 'LineWidth', 1.5, 'FontSize', 9);
        yline(log10(1.0),   'r--', '1 km',  'LineWidth', 1.5, 'FontSize', 9);
        yline(log10(Tolr),  'b--', sprintf('Tolr=%.0e km', Tolr), 'LineWidth', 1.5, 'FontSize', 9);
        hold off;
        fname = sprintf('accuracy_vs_time_%s%s', solvers{s}, suffix);
        saveas(gcf, fullfile(output_dir, [fname '.png']));
        if exist('matlab2tikz', 'file')
            try; matlab2tikz(fullfile(output_dir, [fname '.tikz'])); catch; end
        end
        fprintf('Saved: %s\n', fname);
        close(gcf);
    end

    %% Figure 10: Success Rate vs Revolution Count
    figure;
    nrevs = [0, 1, 2, 3, 4];
    success_by_nrev = zeros(length(solvers), length(nrevs));
    for nr = 1:length(nrevs)
        mask = [testcases.nrev] == nrevs(nr);
        if sum(mask) == 0, continue; end
        for s = 1:length(solvers)
            converged = arrayfun(@(x) x.(solvers{s}).converged, results(mask));
            success_by_nrev(s, nr) = 100 * sum(converged) / sum(mask);
        end
    end
    bar_handle = bar(nrevs, success_by_nrev', 'grouped');
    xlabel('Number of Revolutions', 'FontSize', 12);
    ylabel('Success Rate (%)', 'FontSize', 12);
    legend(solvers, 'Location', 'best', 'FontSize', 11);
    title('Success Rate vs Multi-Revolution Transfers', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    ylim([0 105]);
    set(gca, 'FontSize', 11);
    for s = 1:length(solvers), bar_handle(s).FaceColor = colors(s,:); end
    saveas(gcf, fullfile(output_dir, ['success_vs_revolutions' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['success_vs_revolutions' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: success_vs_revolutions%s\n', suffix);
    close(gcf);

    %% Figure 11: Speedup Comparison
    figure;
    runtime_T     = arrayfun(@(x) x.T.runtime,    results);
    converged_T   = arrayfun(@(x) x.T.converged,  results);
    runtime_MPS   = arrayfun(@(x) x.MPS.runtime,  results);
    converged_MPS = arrayfun(@(x) x.MPS.converged,results);
    runtime_TB    = arrayfun(@(x) x.TB.runtime,   results);
    converged_TB  = arrayfun(@(x) x.TB.converged, results);
    runtime_TRB   = arrayfun(@(x) x.TRB.runtime,  results);
    converged_TRB = arrayfun(@(x) x.TRB.converged,results);
    runtime_MPSB  = arrayfun(@(x) x.MPSB.runtime, results);
    converged_MPSB= arrayfun(@(x) x.MPSB.converged,results);
    runtime_MPSRB = arrayfun(@(x) x.MPSRB.runtime,results);
    converged_MPSRB=arrayfun(@(x) x.MPSRB.converged,results);
    runtime_B     = arrayfun(@(x) x.B.runtime,    results);
    converged_B   = arrayfun(@(x) x.B.converged,  results);
    runtime_BRB   = arrayfun(@(x) x.BRB.runtime,  results);
    converged_BRB = arrayfun(@(x) x.BRB.converged,results);

    comparisons = {
        'T\rightarrowTB',     runtime_T,   converged_T,   runtime_TB,    converged_TB;
        'T\rightarrowTRB',    runtime_T,   converged_T,   runtime_TRB,   converged_TRB;
        'MPS\rightarrowMPSB', runtime_MPS, converged_MPS, runtime_MPSB,  converged_MPSB;
        'MPS\rightarrowMPSRB',runtime_MPS, converged_MPS, runtime_MPSRB, converged_MPSRB;
        'B\rightarrowBRB',    runtime_B,   converged_B,   runtime_BRB,   converged_BRB;
    };
    speedup_data = zeros(5,1);
    speedup_std  = zeros(5,1);
    comp_colors  = [colors(4,:); colors(5,:); colors(6,:); colors(7,:); colors(3,:)];
    for c = 1:5
        rt_base = comparisons{c,2};  cv_base = comparisons{c,3};
        rt_hyb  = comparisons{c,4};  cv_hyb  = comparisons{c,5};
        valid   = cv_base & cv_hyb & ~isnan(rt_base) & ~isnan(rt_hyb) & rt_hyb > 0;
        if sum(valid) > 0
            sp = rt_base(valid) ./ rt_hyb(valid);
            speedup_data(c) = mean(sp);
            speedup_std(c)  = std(sp);
        end
    end
    bar_handle = bar(speedup_data);
    bar_handle.FaceColor = 'flat';
    for c = 1:5, bar_handle.CData(c,:) = comp_colors(c,:); end
    hold on;
    lower_err = min(speedup_std, speedup_data);
    errorbar(1:5, speedup_data, lower_err, speedup_std, 'k.', 'LineWidth', 1.5);
    hold off;
    set(gca, 'XTickLabel', comparisons(:,1), 'FontSize', 10);
    ylabel('Mean Speedup Factor', 'FontSize', 12);
    title('Speedup: Hybrid vs Base Solver', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    yline(1.0, 'r--', 'No Speedup', 'LineWidth', 1.5, 'FontSize', 10);
    saveas(gcf, fullfile(output_dir, ['speedup_comparison' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['speedup_comparison' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: speedup_comparison%s\n', suffix);
    close(gcf);

    %% Figure 12: Iteration Distribution (box plot)
    figure;
    iter_data   = [];
    labels_iter = {};
    for s = 1:length(solvers)
        iters = arrayfun(@(x) x.(solvers{s}).integration_number, results);
        iters = iters(~isnan(iters) & iters > 0);
        iter_data   = [iter_data;   iters(:)]; %#ok<AGROW>
        labels_iter = [labels_iter; repmat(solvers(s), length(iters), 1)]; %#ok<AGROW>
    end
    boxplot(iter_data, labels_iter, 'Colors', colors);
    ylabel('Iteration Count', 'FontSize', 12);
    title('Iteration Count Distribution', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    saveas(gcf, fullfile(output_dir, ['iteration_distribution' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['iteration_distribution' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: iteration_distribution%s\n', suffix);
    close(gcf);

    %% Figure 13: Mean Iterations (bar chart)
    figure;
    mean_iters = zeros(length(solvers), 1);
    std_iters  = zeros(length(solvers), 1);
    for s = 1:length(solvers)
        iters = arrayfun(@(x) x.(solvers{s}).integration_number, results);
        iters = iters(~isnan(iters) & iters > 0);
        if ~isempty(iters)
            mean_iters(s) = mean(iters);
            std_iters(s)  = std(iters);
        end
    end
    bar_handle = bar(mean_iters);
    bar_handle.FaceColor = 'flat';
    for s = 1:length(solvers), bar_handle.CData(s,:) = colors(s,:); end
    hold on;
    lower_err = min(std_iters, mean_iters);
    errorbar(1:length(solvers), mean_iters, lower_err, std_iters, 'k.', 'LineWidth', 1.5);
    hold off;
    set(gca, 'XTickLabel', solvers, 'FontSize', 11);
    ylabel('Mean Iteration Count', 'FontSize', 12);
    title('Average Iterations to Convergence', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    saveas(gcf, fullfile(output_dir, ['mean_iterations' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['mean_iterations' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: mean_iterations%s\n', suffix);
    close(gcf);

    %% Figure 14: Convergence vs Minimum Altitude
    figure;
    min_altitudes = zeros(length(testcases), 1);
    for i = 1:length(testcases)
        r1_alt = norm(testcases(i).r1) - Re;
        r2_alt = norm(testcases(i).r2) - Re;
        min_altitudes(i) = min(r1_alt, r2_alt);
    end
    alt_bins   = [200, 400, 600, 800, 1000, 1500, 2000];
    n_bins     = length(alt_bins) - 1;
    bin_labels = cell(n_bins, 1);
    for b = 1:n_bins
        bin_labels{b} = sprintf('%d-%d', alt_bins(b), alt_bins(b+1));
    end
    conv_by_alt = zeros(length(solvers), n_bins);
    for b = 1:n_bins
        mask     = min_altitudes >= alt_bins(b) & min_altitudes < alt_bins(b+1);
        n_in_bin = sum(mask);
        if n_in_bin == 0, continue; end
        for s = 1:length(solvers)
            converged = arrayfun(@(x) x.(solvers{s}).converged, results(mask));
            conv_by_alt(s, b) = 100 * sum(converged) / n_in_bin;
        end
    end
    bar_handle = bar(1:n_bins, conv_by_alt', 'grouped');
    set(gca, 'XTickLabel', bin_labels, 'FontSize', 10);
    xlabel('Minimum Altitude (km)', 'FontSize', 12);
    ylabel('Convergence Rate (%)', 'FontSize', 12);
    legend(solvers, 'Location', 'best', 'FontSize', 10);
    title('Convergence Rate vs Minimum Altitude', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    ylim([0 105]);
    set(gca, 'FontSize', 11);
    for s = 1:length(solvers), bar_handle(s).FaceColor = colors(s,:); end
    saveas(gcf, fullfile(output_dir, ['convergence_vs_min_altitude' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['convergence_vs_min_altitude' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: convergence_vs_min_altitude%s\n', suffix);
    close(gcf);

    %% Normalized figures — all metrics relative to T (cases where T converged)
    t_conv  = arrayfun(@(x) x.T.converged,          results);
    rt_T    = arrayfun(@(x) x.T.runtime,             results);
    it_T    = arrayfun(@(x) x.T.integration_number, results);
    norm_mask = t_conv & ~isnan(rt_T) & rt_T > 0 & ~isnan(it_T) & it_T > 0;
    res_norm  = results(norm_mask);
    rt_T_n    = rt_T(norm_mask);
    it_T_n    = it_T(norm_mask);

    %% Norm-Fig 1: Runtime Distribution (normalised to T)
    figure;
    rt_norm_data   = [];
    rt_norm_labels = {};
    for s = 1:length(solvers)
        rt_s = arrayfun(@(x) x.(solvers{s}).runtime, res_norm);
        cv_s = arrayfun(@(x) x.(solvers{s}).converged, res_norm);
        ratio = rt_s(cv_s) ./ rt_T_n(cv_s);
        ratio = ratio(isfinite(ratio));
        rt_norm_data   = [rt_norm_data;   ratio(:)];   %#ok<AGROW>
        rt_norm_labels = [rt_norm_labels; repmat(solvers(s), length(ratio), 1)]; %#ok<AGROW>
    end
    boxplot(rt_norm_data, rt_norm_labels, 'Colors', colors);
    ylabel('Runtime / T Runtime', 'FontSize', 12);
    yline(1.0, 'r--', 'T baseline', 'LineWidth', 1.5, 'FontSize', 10);
    title('Normalised Runtime Distribution (relative to T)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    saveas(gcf, fullfile(output_dir, ['runtime_distribution_norm' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['runtime_distribution_norm' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: runtime_distribution_norm%s\n', suffix);
    close(gcf);

    %% Norm-Fig 2: Iteration Distribution (normalised to T)
    figure;
    it_norm_data   = [];
    it_norm_labels = {};
    for s = 1:length(solvers)
        it_s = arrayfun(@(x) x.(solvers{s}).integration_number, res_norm);
        cv_s = arrayfun(@(x) x.(solvers{s}).converged, res_norm);
        ratio = it_s(cv_s) ./ it_T_n(cv_s);
        ratio = ratio(isfinite(ratio));
        it_norm_data   = [it_norm_data;   ratio(:)];   %#ok<AGROW>
        it_norm_labels = [it_norm_labels; repmat(solvers(s), length(ratio), 1)]; %#ok<AGROW>
    end
    boxplot(it_norm_data, it_norm_labels, 'Colors', colors);
    ylabel('Iterations / T Iterations', 'FontSize', 12);
    yline(1.0, 'r--', 'T baseline', 'LineWidth', 1.5, 'FontSize', 10);
    title('Normalised Iteration Distribution (relative to T)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);
    saveas(gcf, fullfile(output_dir, ['iteration_distribution_norm' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['iteration_distribution_norm' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: iteration_distribution_norm%s\n', suffix);
    close(gcf);

    %% Norm-Fig 3: Mean Iterations (normalised to T)
    figure;
    mean_it_norm = zeros(length(solvers), 1);
    std_it_norm  = zeros(length(solvers), 1);
    for s = 1:length(solvers)
        it_s = arrayfun(@(x) x.(solvers{s}).integration_number, res_norm);
        cv_s = arrayfun(@(x) x.(solvers{s}).converged, res_norm);
        ratio = it_s(cv_s) ./ it_T_n(cv_s);
        ratio = ratio(isfinite(ratio));
        if ~isempty(ratio)
            mean_it_norm(s) = mean(ratio);
            std_it_norm(s)  = std(ratio);
        end
    end
    bar_handle = bar(mean_it_norm);
    bar_handle.FaceColor = 'flat';
    for s = 1:length(solvers), bar_handle.CData(s,:) = colors(s,:); end
    hold on;
    lower_err = min(std_it_norm, mean_it_norm);
    errorbar(1:length(solvers), mean_it_norm, lower_err, std_it_norm, 'k.', 'LineWidth', 1.5);
    yline(1.0, 'r--', 'T baseline', 'LineWidth', 1.5, 'FontSize', 10);
    hold off;
    set(gca, 'XTickLabel', solvers, 'FontSize', 11);
    ylabel('Mean Iterations / T Mean Iterations', 'FontSize', 12);
    title('Normalised Mean Iterations (relative to T)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    saveas(gcf, fullfile(output_dir, ['mean_iterations_norm' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['mean_iterations_norm' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: mean_iterations_norm%s\n', suffix);
    close(gcf);

    %% Norm-Fig 4: Convergence by Regime (normalised to T)
    figure;
    conv_rates_norm = zeros(length(solvers), length(regimes));
    for ri = 1:length(regimes)
        rg_mask = strcmp({testcases.regime}, regimes{ri});
        if sum(rg_mask) == 0, continue; end
        t_rate = 100 * sum(arrayfun(@(x) x.T.converged, results(rg_mask))) / sum(rg_mask);
        if t_rate == 0, continue; end
        for s = 1:length(solvers)
            s_rate = 100 * sum(arrayfun(@(x) x.(solvers{s}).converged, results(rg_mask))) / sum(rg_mask);
            conv_rates_norm(s, ri) = s_rate / t_rate;
        end
    end
    present_n = any(conv_rates_norm ~= 0, 1);
    bar_handle = bar(find(present_n), conv_rates_norm(:, present_n)', 'grouped');
    set(gca, 'XTick', find(present_n), 'XTickLabel', regimes(present_n), 'FontSize', 11);
    ylabel('Conv. Rate / T Conv. Rate', 'FontSize', 12);
    legend(solvers, 'Location', 'best', 'FontSize', 10);
    title('Normalised Convergence Rate by Regime (relative to T)', 'FontSize', 14, 'FontWeight', 'bold');
    yline(1.0, 'r--', 'T baseline', 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    for s = 1:length(solvers), bar_handle(s).FaceColor = colors(s,:); end
    saveas(gcf, fullfile(output_dir, ['convergence_by_regime_norm' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['convergence_by_regime_norm' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: convergence_by_regime_norm%s\n', suffix);
    close(gcf);

    %% Norm-Fig 5: Success vs Revolutions (normalised to T)
    figure;
    success_norm = zeros(length(solvers), length(nrevs));
    for nr = 1:length(nrevs)
        mask = [testcases.nrev] == nrevs(nr);
        if sum(mask) == 0, continue; end
        t_rate = 100 * sum(arrayfun(@(x) x.T.converged, results(mask))) / sum(mask);
        if t_rate == 0, continue; end
        for s = 1:length(solvers)
            s_rate = 100 * sum(arrayfun(@(x) x.(solvers{s}).converged, results(mask))) / sum(mask);
            success_norm(s, nr) = s_rate / t_rate;
        end
    end
    bar_handle = bar(nrevs, success_norm', 'grouped');
    xlabel('Number of Revolutions', 'FontSize', 12);
    ylabel('Conv. Rate / T Conv. Rate', 'FontSize', 12);
    legend(solvers, 'Location', 'best', 'FontSize', 11);
    title('Normalised Success Rate vs Revolutions (relative to T)', 'FontSize', 14, 'FontWeight', 'bold');
    yline(1.0, 'r--', 'T baseline', 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
    for s = 1:length(solvers), bar_handle(s).FaceColor = colors(s,:); end
    saveas(gcf, fullfile(output_dir, ['success_vs_revolutions_norm' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['success_vs_revolutions_norm' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: success_vs_revolutions_norm%s\n', suffix);
    close(gcf);

    %% Norm-Fig 6: Convergence vs Min Altitude (normalised to T)
    figure;
    conv_alt_norm = zeros(length(solvers), n_bins);
    for b = 1:n_bins
        mask     = min_altitudes >= alt_bins(b) & min_altitudes < alt_bins(b+1);
        n_in_bin = sum(mask);
        if n_in_bin == 0, continue; end
        t_rate = 100 * sum(arrayfun(@(x) x.T.converged, results(mask))) / n_in_bin;
        if t_rate == 0, continue; end
        for s = 1:length(solvers)
            s_rate = 100 * sum(arrayfun(@(x) x.(solvers{s}).converged, results(mask))) / n_in_bin;
            conv_alt_norm(s, b) = s_rate / t_rate;
        end
    end
    bar_handle = bar(1:n_bins, conv_alt_norm', 'grouped');
    set(gca, 'XTickLabel', bin_labels, 'FontSize', 10);
    xlabel('Minimum Altitude (km)', 'FontSize', 12);
    ylabel('Conv. Rate / T Conv. Rate', 'FontSize', 12);
    legend(solvers, 'Location', 'best', 'FontSize', 10);
    title('Normalised Convergence Rate vs Min Altitude (relative to T)', 'FontSize', 14, 'FontWeight', 'bold');
    yline(1.0, 'r--', 'T baseline', 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
    for s = 1:length(solvers), bar_handle(s).FaceColor = colors(s,:); end
    saveas(gcf, fullfile(output_dir, ['convergence_vs_min_altitude_norm' suffix '.png']));
    if exist('matlab2tikz', 'file')
        try; matlab2tikz(fullfile(output_dir, ['convergence_vs_min_altitude_norm' suffix '.tikz'])); catch; end
    end
    fprintf('Saved: convergence_vs_min_altitude_norm%s\n', suffix);
    close(gcf);

    fprintf('\nVisualization complete! Generated 14 + 6 normalised PNG figures [filter=%s] in %s\n', ...
        filter_type, output_dir);
end
