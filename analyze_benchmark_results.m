function analyze_benchmark_results(varargin)
% ANALYZE_BENCHMARK_RESULTS Comprehensive statistical analysis of benchmark
%
% Usage:
%   analyze_benchmark_results('results_file.mat')
%   analyze_benchmark_results(results, testcases)
%   analyze_benchmark_results(results, testcases, 'output_dir')
%   analyze_benchmark_results(results, testcases, 'output_dir', filter_type)
%
% filter_type : 'all' (default) | 'hard' | 'nonhard'
%
% Generates:
%   - Convergence rates per solver
%   - Failure mode analysis (max iter hit vs integration failure)
%   - Runtime statistics (mean, std, median, 95th percentile)
%   - Accuracy statistics (dr2 mean, std, median, 95th percentile)
%   - Success rates by regime  [skipped for filter_type='hard']
%   - Success rates by revolution count (0-4)
%   - Iteration count statistics
%   - Phase timing breakdown (preamble / phase1 / phase2)
%   - Speedup vs Thompson baseline
%   - Accuracy comparison (all solvers converged)
%   - Solver pair convergence comparison
%
% Output files (suffix = '' | '_hard' | '_nonhard'):
%   benchmark_analysis<suffix>.mat
%   benchmark_analysis<suffix>.csv
%   benchmark_analysis<suffix>.tex

    % Parse input: file path or direct data
    if nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1}))
        results_file = varargin{1};
        fprintf('Loading results from %s...\n', results_file);
        load(results_file, 'results', 'testcases');
        [output_dir, ~, ~] = fileparts(results_file);
        if isempty(output_dir), output_dir = '.'; end
        filter_type = 'all';
    elseif nargin >= 2
        results   = varargin{1};
        testcases = varargin{2};
        output_dir  = 'benchmark_results';
        if nargin >= 3, output_dir  = varargin{3}; end
        filter_type = 'all';
        if nargin >= 4, filter_type = varargin{4}; end
        fprintf('Using data from memory...\n');
    else
        error('Usage: analyze_benchmark_results(file) or analyze_benchmark_results(results, testcases, [output_dir], [filter_type])');
    end

    n_cases = length(results);

    %% Apply filter
    switch filter_type
        case 'hard'
            fmask = strcmp({testcases.regime}, 'hard');
            suffix = '_hard';
        case 'nonhard'
            fmask = ~strcmp({testcases.regime}, 'hard');
            suffix = '_nonhard';
        otherwise
            fmask = true(n_cases, 1);
            suffix = '';
    end
    results   = results(fmask);
    testcases = testcases(fmask);
    n_cases   = sum(fmask);

    solvers = {'T', 'MPS', 'B', 'TB', 'TRB', 'MPSB', 'MPSRB', 'BRB'};
    if nargin >= 5
        maxIter_val = varargin{5};
    else
        maxIter_val = 200;
    end
    maxIter_map = struct('T', maxIter_val, 'MPS', maxIter_val, 'B', maxIter_val, 'TB', maxIter_val, ...
                         'TRB', maxIter_val, 'MPSB', maxIter_val, 'MPSRB', maxIter_val, 'BRB', maxIter_val);

    % Open output files
    fid_csv = fopen(fullfile(output_dir, ['benchmark_analysis' suffix '.csv']), 'w');
    fid_tex = fopen(fullfile(output_dir, ['benchmark_analysis' suffix '.tex']), 'w');

    fprintf('\n========================================\n');
    fprintf('BENCHMARK ANALYSIS: %d Test Cases%s\n', n_cases, ...
        struct('all','','hard',' [Hard only]','nonhard',' [Non-hard only]').(filter_type));
    fprintf('========================================\n\n');

    %% 1. OVERALL CONVERGENCE RATES
    fprintf('--- Convergence Rates ---\n');
    fprintf('%-10s %15s %10s\n', 'Solver', 'Converged', 'Rate');
    fprintf('%s\n', repmat('-', 1, 40));
    fprintf(fid_csv, 'Convergence Rates\n');
    fprintf(fid_csv, 'Solver,Converged,Total,Rate\n');
    fprintf(fid_tex, '%% Convergence Rates\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Converged & Total & Rate \\\\\n\\hline\n');
    for s = 1:length(solvers)
        solver = solvers{s};
        converged = sum(arrayfun(@(x) x.(solver).converged, results));
        rate = 100 * converged / n_cases;
        fprintf('%-10s %6d / %6d   %6.2f%%\n', solver, converged, n_cases, rate);
        fprintf(fid_csv, '%s,%d,%d,%.2f\n', solver, converged, n_cases, rate);
        fprintf(fid_tex, '%s & %d & %d & %.2f\\%% \\\\\n', solver, converged, n_cases, rate);
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 1.5. FAILURE MODE ANALYSIS
    fprintf('--- Failure Mode Analysis ---\n');
    fprintf('%-10s %15s %18s\n', 'Solver', 'Max Iter Hit', 'Integration Fail');
    fprintf('%s\n', repmat('-', 1, 48));
    fprintf(fid_csv, 'Failure Mode Analysis\n');
    fprintf(fid_csv, 'Solver,MaxIterHit,IntegrationFail,TotalFailed\n');
    fprintf(fid_tex, '%% Failure Mode Analysis\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Max Iter Hit & Integration Fail & Total Failed \\\\\n\\hline\n');
    for s = 1:length(solvers)
        solver = solvers{s};
        max_iter = maxIter_map.(solver);
        converged_arr = arrayfun(@(x) x.(solver).converged, results);
        iters_arr     = arrayfun(@(x) x.(solver).integration_number, results);
        failed_mask   = ~converged_arr;
        n_failed      = sum(failed_mask);
        max_iter_hit      = sum(failed_mask & (iters_arr >= max_iter));
        integration_fail  = sum(failed_mask & (iters_arr <  max_iter));
        fprintf('%-10s %15d %18d\n', solver, max_iter_hit, integration_fail);
        fprintf(fid_csv, '%s,%d,%d,%d\n', solver, max_iter_hit, integration_fail, n_failed);
        fprintf(fid_tex, '%s & %d & %d & %d \\\\\n', solver, max_iter_hit, integration_fail, n_failed);
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 2. RUNTIME STATISTICS
    fprintf('--- Runtime Statistics (seconds) ---\n');
    fprintf('%-10s %10s %10s %10s %10s\n', 'Solver', 'Mean', 'Std', 'Median', '95th%');
    fprintf('%s\n', repmat('-', 1, 55));
    fprintf(fid_csv, 'Runtime Statistics (seconds)\n');
    fprintf(fid_csv, 'Solver,Mean,Std,Median,95th\n');
    fprintf(fid_tex, '%% Runtime Statistics (seconds)\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Mean & Std & Median & 95th\\%% \\\\\n\\hline\n');
    for s = 1:length(solvers)
        solver = solvers{s};
        runtimes = arrayfun(@(x) x.(solver).runtime, results);
        runtimes = runtimes(~isnan(runtimes));
        if ~isempty(runtimes)
            m = mean(runtimes); st = std(runtimes);
            med = median(runtimes); p95 = prctile(runtimes, 95);
            fprintf('%-10s %10.3f %10.3f %10.3f %10.3f\n', solver, m, st, med, p95);
            fprintf(fid_csv, '%s,%.4f,%.4f,%.4f,%.4f\n', solver, m, st, med, p95);
            fprintf(fid_tex, '%s & %.3f & %.3f & %.3f & %.3f \\\\\n', solver, m, st, med, p95);
        else
            fprintf('%-10s %10s %10s %10s %10s\n', solver, 'N/A', 'N/A', 'N/A', 'N/A');
            fprintf(fid_csv, '%s,N/A,N/A,N/A,N/A\n', solver);
            fprintf(fid_tex, '%s & N/A & N/A & N/A & N/A \\\\\n', solver);
        end
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 3. ACCURACY STATISTICS (dr2)
    fprintf('--- Position Error Statistics (dr2, km) ---\n');
    fprintf('%-10s %12s %12s %12s %12s\n', 'Solver', 'Mean', 'Std', 'Median', '95th%');
    fprintf('%s\n', repmat('-', 1, 60));
    fprintf(fid_csv, 'Position Error Statistics (dr2 km)\n');
    fprintf(fid_csv, 'Solver,Mean,Std,Median,95th\n');
    fprintf(fid_tex, '%% Position Error Statistics (dr2, km)\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Mean & Std & Median & 95th\\%% \\\\\n\\hline\n');
    for s = 1:length(solvers)
        solver = solvers{s};
        converged = arrayfun(@(x) x.(solver).converged, results);
        dr2_vals  = arrayfun(@(x) x.(solver).dr2, results);
        dr2_vals  = dr2_vals(converged);
        if ~isempty(dr2_vals)
            m = mean(dr2_vals); st = std(dr2_vals);
            med = median(dr2_vals); p95 = prctile(dr2_vals, 95);
            fprintf('%-10s %12.4e %12.4e %12.4e %12.4e\n', solver, m, st, med, p95);
            fprintf(fid_csv, '%s,%.4e,%.4e,%.4e,%.4e\n', solver, m, st, med, p95);
            fprintf(fid_tex, '%s & %.2e & %.2e & %.2e & %.2e \\\\\n', solver, m, st, med, p95);
        else
            fprintf('%-10s %12s %12s %12s %12s\n', solver, 'N/A', 'N/A', 'N/A', 'N/A');
            fprintf(fid_csv, '%s,N/A,N/A,N/A,N/A\n', solver);
            fprintf(fid_tex, '%s & N/A & N/A & N/A & N/A \\\\\n', solver);
        end
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 4. SUCCESS RATE BY REGIME  (skipped for hard-only: single regime, trivial)
    if ~strcmp(filter_type, 'hard')
        fprintf('--- Success Rate by Regime ---\n');
        regimes = {'LEO','MEO','GEO','HEO','hard'};
        fprintf(fid_csv, 'Success Rate by Regime\n');
        fprintf(fid_csv, 'Regime,Solver,SuccessRate,Count\n');
        fprintf(fid_tex, '%% Success Rate by Regime\n');
        fprintf(fid_tex, '\\begin{tabular}{l l r r}\n\\hline\n');
        fprintf(fid_tex, 'Regime & Solver & Success Rate & Count \\\\\n\\hline\n');
        for ri = 1:length(regimes)
            regime  = regimes{ri};
            rg_mask = strcmp({testcases.regime}, regime);
            n_type  = sum(rg_mask);
            if n_type == 0, continue; end
            fprintf('\n%s (%d cases):\n', regime, n_type);
            fprintf('%-10s %10s\n', 'Solver', 'Success Rate');
            fprintf('%s\n', repmat('-', 1, 25));
            for s = 1:length(solvers)
                solver    = solvers{s};
                converged = arrayfun(@(x) x.(solver).converged, results(rg_mask));
                rate      = 100 * sum(converged) / n_type;
                fprintf('%-10s %9.2f%%\n', solver, rate);
                fprintf(fid_csv, '%s,%s,%.2f,%d\n', regime, solver, rate, n_type);
                fprintf(fid_tex, '%s & %s & %.2f\\%% & %d \\\\\n', regime, solver, rate, n_type);
            end
        end
        fprintf('\n');
        fprintf(fid_csv, '\n');
        fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');
    end

    %% 5. PERFORMANCE BY REVOLUTION COUNT
    fprintf('--- Success Rate by Revolution Count ---\n');
    nrevs = [0, 1, 2, 3, 4];
    fprintf(fid_csv, 'Success Rate by Revolution Count\n');
    fprintf(fid_csv, 'Nrev,Solver,SuccessRate,Count\n');
    fprintf(fid_tex, '%% Success Rate by Revolution Count\n');
    fprintf(fid_tex, '\\begin{tabular}{l l r r}\n\\hline\n');
    fprintf(fid_tex, 'Revolutions & Solver & Success Rate & Count \\\\\n\\hline\n');
    for nr = 1:length(nrevs)
        nrev    = nrevs(nr);
        mask    = [testcases.nrev] == nrev;
        n_nrev  = sum(mask);
        if n_nrev == 0, continue; end
        fprintf('\n%d Revolutions (%d cases):\n', nrev, n_nrev);
        fprintf('%-10s %10s\n', 'Solver', 'Success Rate');
        fprintf('%s\n', repmat('-', 1, 25));
        for s = 1:length(solvers)
            solver    = solvers{s};
            converged = arrayfun(@(x) x.(solver).converged, results(mask));
            rate      = 100 * sum(converged) / n_nrev;
            fprintf('%-10s %9.2f%%\n', solver, rate);
            fprintf(fid_csv, '%d,%s,%.2f,%d\n', nrev, solver, rate, n_nrev);
            fprintf(fid_tex, '%d & %s & %.2f\\%% & %d \\\\\n', nrev, solver, rate, n_nrev);
        end
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 6. ITERATION COUNT COMPARISON
    fprintf('--- Iteration Counts (Mean +/- Std) ---\n');
    fprintf('%-10s %20s\n', 'Solver', 'Iterations');
    fprintf('%s\n', repmat('-', 1, 35));
    fprintf(fid_csv, 'Iteration Counts\n');
    fprintf(fid_csv, 'Solver,Mean,Std\n');
    fprintf(fid_tex, '%% Iteration Counts\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Mean & Std \\\\\n\\hline\n');
    for s = 1:length(solvers)
        solver = solvers{s};
        iters  = arrayfun(@(x) x.(solver).integration_number, results);
        iters  = iters(~isnan(iters) & iters > 0);
        if ~isempty(iters)
            m = mean(iters); st = std(iters);
            fprintf('%-10s %10.1f +/- %.1f\n', solver, m, st);
            fprintf(fid_csv, '%s,%.2f,%.2f\n', solver, m, st);
            fprintf(fid_tex, '%s & %.1f & %.1f \\\\\n', solver, m, st);
        else
            fprintf('%-10s %20s\n', solver, 'N/A');
            fprintf(fid_csv, '%s,N/A,N/A\n', solver);
            fprintf(fid_tex, '%s & N/A & N/A \\\\\n', solver);
        end
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 6.5. PHASE TIMING BREAKDOWN
    fprintf('--- Phase Timing Breakdown (converged cases, mean) ---\n');
    fprintf('%-10s %10s %10s %10s %10s %10s %10s\n', ...
        'Solver', 'Pre integ', 'Ph1 integ', 'Ph2 integ', 'Pre t(s)', 'Ph1 t(s)', 'Ph2 t(s)');
    fprintf('%s\n', repmat('-', 1, 75));
    fprintf(fid_csv, 'Phase Timing Breakdown (mean over converged cases)\n');
    fprintf(fid_csv, 'Solver,MeanPreInteg,MeanPh1Integ,MeanPh2Integ,MeanPreRT,MeanPh1RT,MeanPh2RT\n');
    fprintf(fid_tex, '%% Phase Timing Breakdown\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r r r r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Pre integ & Ph1 integ & Ph2 integ & Pre t (s) & Ph1 t (s) & Ph2 t (s) \\\\\n\\hline\n');
    for s = 1:length(solvers)
        solver    = solvers{s};
        converged = arrayfun(@(x) x.(solver).converged, results);
        res_conv  = results(converged);
        if isempty(res_conv)
            fprintf('%-10s %10s %10s %10s %10s %10s %10s\n', solver, 'N/A','N/A','N/A','N/A','N/A','N/A');
            fprintf(fid_csv, '%s,N/A,N/A,N/A,N/A,N/A,N/A\n', solver);
            fprintf(fid_tex, '%s & N/A & N/A & N/A & N/A & N/A & N/A \\\\\n', solver);
            continue;
        end
        mi_pre  = mean(arrayfun(@(x) x.(solver).integration_preamble, res_conv));
        mi_ph1  = mean(arrayfun(@(x) x.(solver).integration_phase1,   res_conv));
        mi_ph2  = mean(arrayfun(@(x) x.(solver).integration_phase2,   res_conv));
        rt_pre  = mean(arrayfun(@(x) x.(solver).runtime_preamble,      res_conv));
        rt_ph1  = mean(arrayfun(@(x) x.(solver).runtime_phase1,        res_conv));
        rt_ph2  = mean(arrayfun(@(x) x.(solver).runtime_phase2,        res_conv));
        fprintf('%-10s %10.1f %10.1f %10.1f %10.4f %10.4f %10.4f\n', ...
            solver, mi_pre, mi_ph1, mi_ph2, rt_pre, rt_ph1, rt_ph2);
        fprintf(fid_csv, '%s,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f\n', ...
            solver, mi_pre, mi_ph1, mi_ph2, rt_pre, rt_ph1, rt_ph2);
        fprintf(fid_tex, '%s & %.1f & %.1f & %.1f & %.4f & %.4f & %.4f \\\\\n', ...
            solver, mi_pre, mi_ph1, mi_ph2, rt_pre, rt_ph1, rt_ph2);
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 7. COMPARATIVE SPEEDUP WITHIN METHOD FAMILIES
    fprintf('--- Speedup: Hybrid vs Base Solver ---\n');
    fprintf('%-15s %-10s %15s %15s %15s\n', 'Comparison', 'Baseline', 'Mean Speedup', 'Median Speedup', 'Std Speedup');
    fprintf('%s\n', repmat('-', 1, 75));
    fprintf(fid_csv, 'Speedup Hybrid vs Base Solver\n');
    fprintf(fid_csv, 'Comparison,Baseline,MeanSpeedup,MedianSpeedup,StdSpeedup\n');
    fprintf(fid_tex, '%% Speedup: Hybrid vs Base Solver\n');
    fprintf(fid_tex, '\\begin{tabular}{l l r r r}\n\\hline\n');
    fprintf(fid_tex, 'Comparison & Baseline & Mean Speedup & Median Speedup & Std Speedup \\\\\n\\hline\n');

    runtime_T    = arrayfun(@(x) x.T.runtime,    results);
    converged_T  = arrayfun(@(x) x.T.converged,  results);
    runtime_MPS  = arrayfun(@(x) x.MPS.runtime,  results);
    converged_MPS= arrayfun(@(x) x.MPS.converged,results);
    runtime_B    = arrayfun(@(x) x.B.runtime,    results);
    converged_B  = arrayfun(@(x) x.B.converged,  results);
    runtime_TB   = arrayfun(@(x) x.TB.runtime,   results);
    converged_TB = arrayfun(@(x) x.TB.converged, results);
    runtime_TRB  = arrayfun(@(x) x.TRB.runtime,  results);
    converged_TRB= arrayfun(@(x) x.TRB.converged,results);
    runtime_MPSB = arrayfun(@(x) x.MPSB.runtime, results);
    converged_MPSB=arrayfun(@(x) x.MPSB.converged,results);
    runtime_MPSRB= arrayfun(@(x) x.MPSRB.runtime,results);
    converged_MPSRB=arrayfun(@(x) x.MPSRB.converged,results);
    runtime_BRB  = arrayfun(@(x) x.BRB.runtime,  results);
    converged_BRB= arrayfun(@(x) x.BRB.converged,results);

    comparisons = {
        'T vs TB',     runtime_T,   converged_T,   runtime_TB,    converged_TB;
        'T vs TRB',    runtime_T,   converged_T,   runtime_TRB,   converged_TRB;
        'MPS vs MPSB', runtime_MPS, converged_MPS, runtime_MPSB,  converged_MPSB;
        'MPS vs MPSRB',runtime_MPS, converged_MPS, runtime_MPSRB, converged_MPSRB;
        'B vs BRB',    runtime_B,   converged_B,   runtime_BRB,   converged_BRB;
    };

    for c = 1:size(comparisons, 1)
        label   = comparisons{c,1};
        rt_base = comparisons{c,2};  cv_base = comparisons{c,3};
        rt_hyb  = comparisons{c,4};  cv_hyb  = comparisons{c,5};
        parts   = strsplit(label, ' vs ');
        base_name = parts{1};
        valid_mask = cv_base & cv_hyb & ~isnan(rt_base) & ~isnan(rt_hyb);
        if sum(valid_mask) > 0
            speedup = rt_base(valid_mask) ./ rt_hyb(valid_mask);
            m_sp = mean(speedup); med_sp = median(speedup); std_sp = std(speedup);
            fprintf('%-15s %-10s %15.3fx %15.3fx %15.3f\n', label, base_name, m_sp, med_sp, std_sp);
            fprintf(fid_csv, '%s,%s,%.3f,%.3f,%.3f\n', label, base_name, m_sp, med_sp, std_sp);
            fprintf(fid_tex, '%s & %s & %.3fx & %.3fx & %.3f \\\\\n', label, base_name, m_sp, med_sp, std_sp);
        else
            fprintf('%-15s %-10s %15s %15s %15s\n', label, base_name, 'N/A', 'N/A', 'N/A');
            fprintf(fid_csv, '%s,%s,N/A,N/A,N/A\n', label, base_name);
            fprintf(fid_tex, '%s & %s & N/A & N/A & N/A \\\\\n', label, base_name);
        end
    end
    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 8. ACCURACY COMPARISON (Cases where all solvers converged)
    fprintf('--- Accuracy Comparison (All Converged) ---\n');
    all_converged_mask = arrayfun(@(x) x.TB.converged,   results) & ...
                         arrayfun(@(x) x.TRB.converged,  results) & ...
                         arrayfun(@(x) x.MPSB.converged, results) & ...
                         arrayfun(@(x) x.MPSRB.converged,results);
    n_all_converged = sum(all_converged_mask);
    fprintf('Cases where all 4 hybrid solvers converged: %d (%.1f%%)\n\n', ...
            n_all_converged, 100*n_all_converged/n_cases);
    fprintf(fid_csv, 'Accuracy Comparison (All Converged)\n');
    fprintf(fid_csv, 'AllConvergedCases,%d,%.1f%%\n', n_all_converged, 100*n_all_converged/n_cases);
    fprintf(fid_csv, 'Solver,MeanDr2,MedianDr2\n');
    fprintf(fid_tex, '%% Accuracy Comparison (All Converged)\n');
    fprintf(fid_tex, '%% Cases where all 4 hybrid solvers converged: %d (%.1f%%)\n', n_all_converged, 100*n_all_converged/n_cases);
    fprintf(fid_tex, '\\begin{tabular}{l r r}\n\\hline\n');
    fprintf(fid_tex, 'Solver & Mean dr2 (km) & Median dr2 (km) \\\\\n\\hline\n');
    if n_all_converged > 0
        fprintf('%-10s %15s %15s\n', 'Solver', 'Mean dr2 (km)', 'Median dr2 (km)');
        fprintf('%s\n', repmat('-', 1, 45));
        for s = 1:length(solvers)
            solver  = solvers{s};
            dr2_all = arrayfun(@(x) x.(solver).dr2, results(all_converged_mask));
            m = mean(dr2_all); med = median(dr2_all);
            fprintf('%-10s %15.4e %15.4e\n', solver, m, med);
            fprintf(fid_csv, '%s,%.4e,%.4e\n', solver, m, med);
            fprintf(fid_tex, '%s & %.2e & %.2e \\\\\n', solver, m, med);
        end
        fprintf('\n');
    end
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 8.5. SOLVER PAIR CONVERGENCE COMPARISON
    fprintf('--- Solver Pair Convergence Comparison ---\n\n');
    fprintf(fid_csv, 'Solver Pair Convergence Comparison\n');
    fprintf(fid_csv, 'Comparison,BothConverged,FirstOnly,SecondOnly\n');
    fprintf(fid_tex, '%% Solver Pair Convergence Comparison\n');
    fprintf(fid_tex, '\\begin{tabular}{l r r r}\n\\hline\n');
    fprintf(fid_tex, 'Comparison & Both Converged & First Only & Second Only \\\\\n\\hline\n');

    all_converged_ids = find(all_converged_mask);
    fprintf('Cases where ALL 4 hybrid solvers converged: %d cases\n', length(all_converged_ids));
    fprintf('\n');

    conv_TB  = arrayfun(@(x) x.TB.converged,  results);
    conv_TRB = arrayfun(@(x) x.TRB.converged, results);
    tb_only   = find(conv_TB  & ~conv_TRB);
    trb_only  = find(~conv_TB &  conv_TRB);
    tb_trb_both = find(conv_TB & conv_TRB);
    fprintf('Thompson Family (TB vs TRB):\n');
    fprintf('  Both converged: %d  |  TB only: %d  |  TRB only: %d\n', ...
        length(tb_trb_both), length(tb_only), length(trb_only));
    fprintf(fid_csv, 'Thompson (TB vs TRB),%d,%d,%d\n', length(tb_trb_both), length(tb_only), length(trb_only));
    fprintf(fid_tex, 'Thompson (TB vs TRB) & %d & %d & %d \\\\\n', length(tb_trb_both), length(tb_only), length(trb_only));

    conv_MPSB  = arrayfun(@(x) x.MPSB.converged,  results);
    conv_MPSRB = arrayfun(@(x) x.MPSRB.converged, results);
    mpsb_only   = find(conv_MPSB  & ~conv_MPSRB);
    mpsrb_only  = find(~conv_MPSB &  conv_MPSRB);
    mpsb_mpsrb_both = find(conv_MPSB & conv_MPSRB);
    fprintf('MPS Family (MPSB vs MPSRB):\n');
    fprintf('  Both converged: %d  |  MPSB only: %d  |  MPSRB only: %d\n', ...
        length(mpsb_mpsrb_both), length(mpsb_only), length(mpsrb_only));
    fprintf(fid_csv, 'MPS (MPSB vs MPSRB),%d,%d,%d\n', length(mpsb_mpsrb_both), length(mpsb_only), length(mpsrb_only));
    fprintf(fid_tex, 'MPS (MPSB vs MPSRB) & %d & %d & %d \\\\\n', length(mpsb_mpsrb_both), length(mpsb_only), length(mpsrb_only));

    fprintf('\n');
    fprintf(fid_csv, '\n');
    fprintf(fid_tex, '\\hline\n\\end{tabular}\n\n');

    %% 9. SAVE SUMMARY DATA
    fprintf('Saving analysis summary...\n');
    summary = struct();
    summary.n_cases    = n_cases;
    summary.solvers    = solvers;
    summary.filter_type = filter_type;

    for s = 1:length(solvers)
        solver = solvers{s};
        summary.convergence.(solver) = sum(arrayfun(@(x) x.(solver).converged, results)) / n_cases;
    end

    for s = 1:length(solvers)
        solver = solvers{s};
        max_iter      = maxIter_map.(solver);
        converged_arr = arrayfun(@(x) x.(solver).converged, results);
        iters_arr     = arrayfun(@(x) x.(solver).integration_number, results);
        failed_mask   = ~converged_arr;
        summary.failure_mode.(solver).max_iter_hit        = sum(failed_mask & (iters_arr >= max_iter));
        summary.failure_mode.(solver).integration_fail    = sum(failed_mask & (iters_arr <  max_iter));
    end

    for s = 1:length(solvers)
        solver   = solvers{s};
        runtimes = arrayfun(@(x) x.(solver).runtime, results);
        runtimes = runtimes(~isnan(runtimes));
        summary.runtime.(solver) = struct('mean', mean(runtimes), 'std', std(runtimes), ...
            'median', median(runtimes), 'p95', prctile(runtimes, 95));
    end

    for s = 1:length(solvers)
        solver   = solvers{s};
        converged = arrayfun(@(x) x.(solver).converged, results);
        dr2_vals  = arrayfun(@(x) x.(solver).dr2, results);
        dr2_vals  = dr2_vals(converged);
        if ~isempty(dr2_vals)
            summary.accuracy.(solver) = struct('mean', mean(dr2_vals), 'std', std(dr2_vals), ...
                'median', median(dr2_vals), 'p95', prctile(dr2_vals, 95));
        end
    end

    if ~strcmp(filter_type, 'hard')
        regimes = {'LEO','MEO','GEO','HEO','hard'};
        for ri = 1:length(regimes)
            rg      = regimes{ri};
            rg_mask = strcmp({testcases.regime}, rg);
            if sum(rg_mask) == 0, continue; end
            summary.by_regime.(rg) = struct();
            for s = 1:length(solvers)
                solver    = solvers{s};
                converged = arrayfun(@(x) x.(solver).converged, results(rg_mask));
                summary.by_regime.(rg).(solver) = sum(converged) / max(sum(rg_mask), 1);
            end
        end
    end

    for nr = 1:length(nrevs)
        nrev = nrevs(nr);
        mask = [testcases.nrev] == nrev;
        if sum(mask) == 0, continue; end
        field_name = sprintf('nrev_%d', nrev);
        summary.by_nrev.(field_name) = struct();
        for s = 1:length(solvers)
            solver    = solvers{s};
            converged = arrayfun(@(x) x.(solver).converged, results(mask));
            summary.by_nrev.(field_name).(solver) = sum(converged) / sum(mask);
        end
    end

    save(fullfile(output_dir, ['benchmark_analysis_summary' suffix '.mat']), 'summary');
    fprintf('Analysis summary saved to %s\n', fullfile(output_dir, ['benchmark_analysis_summary' suffix '.mat']));

    fclose(fid_csv);
    fclose(fid_tex);
    fprintf('CSV saved to %s\n', fullfile(output_dir, ['benchmark_analysis' suffix '.csv']));
    fprintf('LaTeX saved to %s\n\n', fullfile(output_dir, ['benchmark_analysis' suffix '.tex']));

    fprintf('========================================\n');
    fprintf('ANALYSIS COMPLETE\n');
    fprintf('========================================\n');
end
