function generate_benchmark_testcases(n_total, output_file, AbsTol, RelTol, N, delta, maxIter, use_parallel)
% GENERATE_BENCHMARK_TESTCASES  Generate Lambert problem test cases.
%
% Usage:
%   generate_benchmark_testcases()
%   generate_benchmark_testcases(75000, 'testcases_75k.mat')
%   generate_benchmark_testcases(75000, 'testcases_75k.mat', ..., true)   % parallel (default if PCT available)
%   generate_benchmark_testcases(75000, 'testcases_75k.mat', ..., false)  % serial
%
% Distribution: 20% per nrev (0-4). Within each nrev bin:
%   15% LEO  (rp_alt & ra_alt in [350, 2000) km)
%   15% MEO  (rp_alt & ra_alt in [2000, 20200) km)
%   15% GEO  (rp_alt & ra_alt in [35286, 36286] km)
%   15% HEO  (rp_alt in [350, 1000) km, ra_alt in [20000, 42000) km)
%   40% hard (rp_alt in [200, 350) km, ra_alt in [rp_alt, 42000) km)
%
% Struct fields:
%   id, r1, v1, r2, v2, m, A, tm, dm, de, nrev, is_hard, regime

    if nargin < 1 || isempty(n_total),     n_total     = 75000;               end
    if nargin < 2 || isempty(output_file), output_file = 'testcases_75k.mat'; end
    if nargin < 3 || isempty(AbsTol),      AbsTol      = 1e-10;               end
    if nargin < 4 || isempty(RelTol),      RelTol      = 1e-12;               end
    if nargin < 5 || isempty(N),           N           = 16;                 end
    if nargin < 6 || isempty(delta),       delta       = 32;                 end
    if nargin < 7 || isempty(maxIter),     maxIter     = 200;                end
    if nargin < 8 || isempty(use_parallel)
        use_parallel = license('test', 'Distrib_Computing_Toolbox');
    end

    addpath(genpath('vallado/software/matlab'));

    c.mu      = 398600.5;
    c.Re      = 6378.137;
    c.J       = [0, 0.01082635854, -0.2532435346e-5, -0.1619331205e-5, ...
                 -0.2277161016e-6, 0.5396484906e-6];
    c.Cd      = 2.2;
    c.Cr      = 1.2;
    c.muM     = 4902.8;
    c.rhoSRP  = 4.56e-3;
    c.jdepoch = juliandate(2026, 1, 1, 12, 0, 0);

    ode_opts = struct('AbsTol',AbsTol,'RelTol',RelTol,'maxIter',maxIter,'N',N,'delta',delta);

    regimes  = {'LEO','MEO','GEO','HEO','hard'};
    fracs    = [0.15, 0.15, 0.15, 0.15, 0.40];
    nrev_set = 0:4;
    n_per_nrev = floor(n_total / 5);
    remainder  = n_total - n_per_nrev * 5;

    % Build flat list of 25 independent cells (5 nrev × 5 regime)
    n_cells      = 25;
    cell_nrev    = zeros(n_cells, 1);
    cell_regime  = cell(n_cells, 1);
    cell_ntarget = zeros(n_cells, 1);
    ci = 0;
    for ni = 1:5
        nrev_val = nrev_set(ni);
        n_this   = n_per_nrev + (ni <= remainder);
        for ri = 1:5
            ci = ci + 1;
            cell_nrev(ci)    = nrev_val;
            cell_regime{ci}  = regimes{ri};
            cell_ntarget(ci) = round(n_this * fracs(ri));
        end
    end

    tc_cells = cell(n_cells, 1);

    warning_state = warning('off', 'all');

    if use_parallel
        fprintf('Generating %d cases across %d cells (parallel)...\n\n', n_total, n_cells);
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(msg) fprintf('%s', msg));
        parfor ci = 1:n_cells
            warning('off', 'all');
            tc_cells{ci} = generate_cases(cell_ntarget(ci), cell_regime{ci}, ...
                                          cell_nrev(ci), c, ode_opts, dq, ci);
        end
    else
        fprintf('Generating %d cases across %d cells (serial)...\n\n', n_total, n_cells);
        for ci = 1:n_cells
            tc_cells{ci} = generate_cases(cell_ntarget(ci), cell_regime{ci}, ...
                                          cell_nrev(ci), c, ode_opts, [], ci);
        end
    end

    warning(warning_state);

    all_tc = vertcat(tc_cells{:});

    rng(42);
    all_tc = all_tc(randperm(length(all_tc)));
    for i = 1:length(all_tc)
        all_tc(i).id = i;
    end

    testcases = all_tc;
    fprintf('\nSaving to %s...\n', output_file);
    save(output_file, 'testcases', '-v7.3');
    print_statistics(testcases, c.Re);
    fprintf('\nTest case generation complete!\n');
end

% -------------------------------------------------------------------------

function tc_list = generate_cases(n_target, regime, nrev, c, ode_opts, dq, ci)
% Generate exactly n_target cases for the given regime and nrev.
% dq: DataQueue for real-time parallel progress, or [] to use fprintf directly.
% ci: cell index (1-25) used in progress labels.

    is_hard = strcmp(regime, 'hard');
    progress_interval = max(50, floor(n_target / 10));

    if isempty(dq)
        out = @(varargin) fprintf(varargin{:});
    else
        out = @(varargin) send(dq, sprintf(varargin{:}));
    end

    out('[%02d] %s nrev=%d  start  (target %d)\n', ci, regime, nrev, n_target);

    template = struct('id', 0, 'r1', zeros(1,3), 'v1', zeros(1,3), ...
                      'r2', zeros(1,3), 'v2', zeros(1,3), ...
                      'm', 0, 'A', 0, 'tm', 0, 'dm', '', 'de', '', ...
                      'nrev', 0, 'is_hard', false, 'regime', '');

    tc_list  = repmat(template, n_target, 1);
    accepted = 0;
    attempted = 0;

    while accepted < n_target
        attempted = attempted + 1;

        % Satellite properties
        m        = 0.1 + rand() * (500 - 0.1);
        A_over_m = 0.001 + rand() * (0.013 - 0.001);
        A        = A_over_m * m * 1e-6;

        % Orbital elements — regime-dependent
        switch regime
            case 'LEO'
                rp_alt = 350  + rand() * (2000 - 350);
                ra_alt = rp_alt + rand() * (2000 - rp_alt);
            case 'MEO'
                rp_alt = 2000 + rand() * (20200 - 2000);
                ra_alt = rp_alt + rand() * (20200 - rp_alt);
            case 'GEO'
                rp_alt = 35286 + rand() * 1000;
                ra_alt = 35286 + rand() * 1000;
                if ra_alt < rp_alt, [rp_alt, ra_alt] = deal(ra_alt, rp_alt); end
            case 'HEO'
                rp_alt = 350   + rand() * (1000 - 350);
                ra_alt = 20000 + rand() * (42000 - 20000);
            case 'hard'
                rp_alt = 200 + rand() * (350 - 200);
                ra_alt = rp_alt + rand() * (42000 - rp_alt);
        end

        rp  = c.Re + rp_alt;
        ra  = c.Re + ra_alt;
        a   = (rp + ra) / 2;
        ecc = (ra - rp) / (ra + rp);

        incl  = rand() * pi;
        Omega = rand() * 2 * pi;
        argp  = rand() * 2 * pi;
        nu0   = rand() * 2 * pi;

        T_orb = 2 * pi * sqrt(a^3 / c.mu);
        tm    = T_orb * (nrev + rand());
        ode_opts.Sec = T_orb / ode_opts.delta;

        [r0, v0] = posnvelos(a, ecc, incl, Omega, argp, nu0, c.mu);

        forcemodel = @(t,x) orbit_eq_J6_drag_SRP_moon(t, x, c.mu, c.Cd, A, m, ...
                                c.Re, c.J, c.jdepoch, c.rhoSRP, c.Cr, A, c.muM);

        try
            [~, xout, err] = odeMPCI(forcemodel, [0, tm], [r0, v0], ode_opts);
        catch
            continue;
        end
        if err == -1 || isempty(xout)
            continue;
        end

        rf = xout(end, 1:3);
        vf = xout(end, 4:6);

        % Earth-impact check (analytical perigee)
        [he, ~] = hitearth(100, r0, v0, rf, vf, nrev);
        if he
            continue;
        end

        h = cross(r0, v0);
        if dot(cross(r0, rf), h) > 0
            dm = 'S';
        else
            dm = 'L';
        end

        c_chord = norm(rf - r0);
        s       = (norm(r0) + norm(rf) + c_chord) / 2;
        if s <= 0 || (s - c_chord) < 0
            continue;
        end
        a_minE  = s / 2;
        beta    = 2 * asin(sqrt((s - c_chord) / s));
        sgn = 1;  if strcmp(dm, 'L'), sgn = -1; end
        dt_minE = sqrt(a_minE^3 / c.mu) * (2*nrev*pi + pi - sgn*(beta - sin(beta)));
        if tm < dt_minE,  de = 'H';  else,  de = 'L';  end

        accepted = accepted + 1;
        tc_list(accepted).r1      = r0;
        tc_list(accepted).v1      = v0;
        tc_list(accepted).r2      = rf;
        tc_list(accepted).v2      = vf;
        tc_list(accepted).m       = m;
        tc_list(accepted).A       = A;
        tc_list(accepted).tm      = tm;
        tc_list(accepted).dm      = dm;
        tc_list(accepted).de      = de;
        tc_list(accepted).nrev    = nrev;
        tc_list(accepted).is_hard = is_hard;
        tc_list(accepted).regime  = regime;

        if mod(accepted, progress_interval) == 0
            out('[%02d] %s nrev=%d  %d/%d (%.0f%%)  reject=%.0f%%\n', ...
                ci, regime, nrev, accepted, n_target, ...
                100*accepted/n_target, ...
                100*(attempted - accepted)/attempted);
        end
    end

    tc_list = tc_list(1:accepted);
    out('[%02d] %s nrev=%d  DONE  %d/%d attempted  reject=%.1f%%\n', ...
        ci, regime, nrev, accepted, attempted, 100*(attempted - accepted)/attempted);
end

% -------------------------------------------------------------------------

function print_statistics(testcases, Re)
    fprintf('\n========================================\n');
    fprintf('TEST CASE STATISTICS\n');
    fprintf('========================================\n\n');
    n = length(testcases);
    fprintf('Total test cases: %d\n\n', n);

    n_hard = sum([testcases.is_hard]);
    fprintf('Hard cases: %d (%.1f%%)  Non-hard: %d (%.1f%%)\n\n', ...
            n_hard, 100*n_hard/n, n-n_hard, 100*(n-n_hard)/n);

    fprintf('Revolution count distribution:\n');
    for nrev = 0:4
        k = sum([testcases.nrev] == nrev);
        fprintf('  nrev=%d: %d (%.1f%%)\n', nrev, k, 100*k/n);
    end
    fprintf('\n');

    fprintf('Regime distribution:\n');
    for rg = {'LEO','MEO','GEO','HEO','hard'}
        k = sum(strcmp({testcases.regime}, rg{1}));
        fprintf('  %-4s: %d (%.1f%%)\n', rg{1}, k, 100*k/n);
    end
    fprintf('\n');

    n_S  = sum(strcmp({testcases.dm}, 'S'));
    n_L  = sum(strcmp({testcases.dm}, 'L'));
    n_He = sum(strcmp({testcases.de}, 'H'));
    n_Le = sum(strcmp({testcases.de}, 'L'));
    fprintf('Direction mode:  S=%d (%.1f%%)  L=%d (%.1f%%)\n', ...
            n_S, 100*n_S/n, n_L, 100*n_L/n);
    fprintf('Energy mode:     H=%d (%.1f%%)  L=%d (%.1f%%)\n\n', ...
            n_He, 100*n_He/n, n_Le, 100*n_Le/n);

    alt1 = arrayfun(@(tc) norm(tc.r1) - Re, testcases);
    alt2 = arrayfun(@(tc) norm(tc.r2) - Re, testcases);
    fprintf('Departure altitude (km):  min=%.1f  mean=%.1f  max=%.1f\n', ...
            min(alt1), mean(alt1), max(alt1));
    fprintf('Arrival altitude (km):    min=%.1f  mean=%.1f  max=%.1f\n\n', ...
            min(alt2), mean(alt2), max(alt2));

    tms = [testcases.tm] / 3600;
    fprintf('Transfer time (h):  min=%.2f  mean=%.2f  max=%.2f\n\n', ...
            min(tms), mean(tms), max(tms));

    masses   = [testcases.m];
    areas_m2 = [testcases.A] * 1e6;
    fprintf('Mass (kg):   min=%.2f  mean=%.1f  max=%.1f\n', ...
            min(masses), mean(masses), max(masses));
    fprintf('Area (m²):   min=%.5f  mean=%.5f  max=%.5f\n', ...
            min(areas_m2), mean(areas_m2), max(areas_m2));
end
