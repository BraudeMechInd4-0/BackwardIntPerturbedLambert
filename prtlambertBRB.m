function [v1t, v2t, stats, history] = prtlambertBRB(forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTBRB  Perturbed Lambert solver — Broyden (select fwd/bwd) then Broyden.
%
% Preamble: propagate forward from r1 and backward from r2; choose the
%           direction with the smaller initial residual.
% Phase 1:  Broyden quasi-Newton in the chosen direction, limited to N_switch
%           iterations or until normdr < T_switch. J0 = eye(3).
% Phase 2:  Broyden forward, initialized with J from Phase 1 (plus one
%           secant refinement from the lambertb seed to the Phase 1 exit).
%
% History always tracks the FORWARD error (normdr2), even during backward phase.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertBRB( ...
%       forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%
% Outputs:
%   v1t     - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t     - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   stats   - struct: integration_number, integration_preamble, integration_phase1,
%             integration_phase2, runtime_preamble, runtime_phase1, runtime_phase2,
%             direction_chosen ('forward'|'backward'), hitearthvar, hitrad
%   history - struct (debug mode only)

mu = 398600.5;
r_GEO = 42164;
r_max_allowed = 1.5 * r_GEO;

if nargin < 9,  options = struct();  end
if ~isfield(options, 'Tolr'),    options.Tolr    = 1e-8;  end
if ~isfield(options, 'maxIter'), options.maxIter = 100;   end
if ~isfield(options, 'delta'),   options.delta   = 16;    end

Tolr    = options.Tolr;
maxIter = options.maxIter;
do_print      = isfield(options, 'debug') && options.debug;
build_history = nargout >= 4;

T_switch = Tolr * 1000;
N_switch = floor(maxIter * 0.05);
if isfield(options, 'T_switch'),  T_switch = options.T_switch;  end
if isfield(options, 'N_switch'),  N_switch = options.N_switch;  end

if do_print, fprintf('starting BRB\n'); end
if build_history
    history = struct('iterations', [], 'normdr2', [], 'phase', {{}}, ...
                     'direction_chosen', '', 'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_pre = tic;

%% ---- Preamble: lambertb + forward + backward propagations ----
[v1tnew, v2tnew] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1tnew)) || any(isnan(v1tnew)) || ~all(isreal(v2tnew)) || any(isnan(v2tnew))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_preamble = toc(t_pre);
    return
end

oe = kep_elements(r1, v1tnew, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options_fwd     = options;
options_fwd.Sec = T_orb / options.delta;
[~, xoutfwd]    = odeMPCI(forcemodel, [0 dtsec], [r1 v1tnew], options_fwd);

oe = kep_elements(r2, v2tnew, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options_bkw     = options;
options_bkw.Sec = T_orb / options.delta;
[~, xoutbkw]    = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2tnew], options_bkw);

integration_number = 2;
normdr1 = norm(xoutbkw(1,1:3) - r1);
normdr2 = norm(xoutfwd(end,1:3) - r2);

v1_A  = v1tnew;
dr2_A = xoutfwd(end,1:3) - r2;

rt_pre = toc(t_pre);
stats.integration_preamble = 2;
stats.runtime_preamble     = rt_pre;

if do_print, fprintf('init (BRB): normdr1=%g, normdr2=%g\n', normdr1, normdr2); end
if build_history
    history.iterations(end+1) = integration_number;
    history.normdr2(end+1)    = normdr2;
    history.phase{end+1}      = 'init';
end

%% ---- Select direction ----
use_backward = isfinite(normdr1) && (normdr1 < normdr2);
if use_backward
    stats.direction_chosen = 'backward';
else
    stats.direction_chosen = 'forward';
end
if do_print, fprintf('BRB: using %s\n', stats.direction_chosen); end
if build_history, history.direction_chosen = stats.direction_chosen; end

t_ph1 = tic;

%% ====================================================================
%% BACKWARD BRANCH — Broyden on v2, minimising dr1
%% ====================================================================
if use_backward
    v2_cur      = v2tnew;
    dr1_cur     = xoutbkw(1,1:3) - r1;
    xoutbkw_cur = xoutbkw;

    delta_size  = 1e-7;
    normv0_bkw  = norm(v2_cur);
    J_rev       = zeros(3);
    for jj = 1:3
        v2_pert     = v2_cur;
        v2_pert(jj) = v2_pert(jj) + delta_size * normv0_bkw;
        oe_p = kep_elements(r2, v2_pert, mu);  a_p = oe_p(1);
        if a_p <= 0 || ~isreal(a_p) || isnan(a_p) || a_p > 2*r_max_allowed
            T_p = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            T_p = sqrt(a_p^3/mu) * 2*pi;
        end
        opt_p = options_bkw;  opt_p.Sec = T_p / options.delta;
        [~, x_pert_bkw] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2_pert], opt_p);
        J_rev(:,jj) = (x_pert_bkw(1,1:3) - xoutbkw_cur(1,1:3))' / (delta_size * normv0_bkw);
    end
    integration_number = integration_number + 3;
    if ~all(isfinite(J_rev(:))),  J_rev = eye(3);  end

    n_broyden   = 0;
    normdr1_prev = norm(dr1_cur);

    while norm(dr1_cur) > T_switch && ...
            (n_broyden < N_switch || norm(dr1_cur) > normdr1_prev) && ...
            integration_number < maxIter

        v2_prev  = v2_cur;
        dr1_prev = dr1_cur;
        v2_cur   = v2_cur - (J_rev \ dr1_cur')';

        if ~all(isreal(v2_cur)) || any(~isfinite(v2_cur))
            if build_history, history.converged = false; end
            stats.integration_number = integration_number;
            stats.integration_phase1 = integration_number - 2;
            stats.runtime_phase1     = toc(t_ph1);
            v1t = NaN(1,3);  v2t = NaN(1,3);
            return
        end

        oe = kep_elements(r2, v2_cur, mu);
        a  = oe(1);
        if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
            T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            T_orb = sqrt(a^3/mu) * 2*pi;
        end
        options_bkw.Sec = T_orb / options.delta;
        [~, xoutbkw_cur] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2_cur], options_bkw);
        integration_number = integration_number + 1;
        if ~all(isfinite(xoutbkw_cur(1,:))), break, end

        normdr1_prev = norm(dr1_cur);
        dr1_cur = xoutbkw_cur(1,1:3) - r1;

        dx_r = (v2_cur  - v2_prev)';
        df_r = (dr1_cur - dr1_prev)';
        if norm(dx_r) > 0
            J_rev = J_rev + (df_r - J_rev*dx_r) * dx_r' / (dx_r'*dx_r);
        end
        n_broyden = n_broyden + 1;

        if build_history
            v1_check = xoutbkw_cur(1,4:6);
            oe_c = kep_elements(r1, v1_check, mu);
            a_c  = oe_c(1);
            if a_c <= 0 || ~isreal(a_c) || isnan(a_c) || a_c > 2*r_max_allowed
                T_c = 2*pi*sqrt((2*r_max_allowed)^3/mu);
            else
                T_c = sqrt(a_c^3/mu) * 2*pi;
            end
            opt_c = options;  opt_c.Sec = T_c / options.delta;
            [~, xfwd_c] = odeMPCI(forcemodel, [0 dtsec], [r1 v1_check], opt_c);
            integration_number = integration_number + 1;
            normdr2_check = norm(xfwd_c(end,1:3) - r2);
            history.iterations(end+1) = integration_number;
            history.normdr2(end+1)    = normdr2_check;
            history.phase{end+1}      = 'Broyden_bkw';
            if do_print, fprintf('integration %d (B bkw): normdr1=%g, normdr2=%g\n', ...
                    integration_number, norm(dr1_cur), normdr2_check); end
        elseif do_print
            fprintf('iteration %d (B bkw): normdr1=%g\n', n_broyden, norm(dr1_cur));
        end
    end

    %% Transition: backward J → forward J (same as TRB line 209 and MPSRB line 195)
    J = -J_rev';

    v1_B = xoutbkw_cur(1,4:6);
    oe = kep_elements(r1, v1_B, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options_fwd.Sec = T_orb / options.delta;
    [~, xout_B] = odeMPCI(forcemodel, [0 dtsec], [r1 v1_B], options_fwd);
    integration_number = integration_number + 1;
    dr2_B   = xout_B(end,1:3) - r2;
    xout    = xout_B;
    v1_cur  = v1_B;
    dr2_cur = dr2_B;

    if do_print, fprintf('switch to Broyden (bkw): normdr2=%g\n', norm(dr2_cur)); end

%% ====================================================================
%% FORWARD BRANCH — Broyden on v1, minimising dr2
%% ====================================================================
else
    v1_cur      = v1tnew;
    xout        = xoutfwd;
    dr2_cur     = xoutfwd(end,1:3) - r2;

    delta_size  = 1e-7;
    normv0_fwd  = norm(v1_cur);
    J           = zeros(3);
    for jj = 1:3
        v1_pert     = v1_cur;
        v1_pert(jj) = v1_pert(jj) + delta_size * normv0_fwd;
        oe_p = kep_elements(r1, v1_pert, mu);  a_p = oe_p(1);
        if a_p <= 0 || ~isreal(a_p) || isnan(a_p) || a_p > 2*r_max_allowed
            T_p = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            T_p = sqrt(a_p^3/mu) * 2*pi;
        end
        opt_p = options_fwd;  opt_p.Sec = T_p / options.delta;
        [~, x_pert_fwd] = odeMPCI(forcemodel, [0 dtsec], [r1 v1_pert], opt_p);
        J(:,jj) = (x_pert_fwd(end,1:3) - xout(end,1:3))' / (delta_size * normv0_fwd);
    end
    integration_number = integration_number + 3;
    if ~all(isfinite(J(:))),  J = eye(3);  end

    n_broyden   = 0;
    normdr2_prev = norm(dr2_cur);

    while norm(dr2_cur) > T_switch && ...
            (n_broyden < N_switch || norm(dr2_cur) > normdr2_prev) && ...
            integration_number < maxIter

        v1_prev  = v1_cur;
        dr2_prev = dr2_cur;
        v1_cur   = v1_cur - (J \ dr2_cur')';

        if ~all(isreal(v1_cur)) || any(~isfinite(v1_cur))
            if build_history, history.converged = false; end
            stats.integration_number = integration_number;
            stats.integration_phase1 = integration_number - 2;
            stats.runtime_phase1     = toc(t_ph1);
            v1t = NaN(1,3);  v2t = NaN(1,3);
            return
        end

        oe = kep_elements(r1, v1_cur, mu);
        a  = oe(1);
        if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
            T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            T_orb = sqrt(a^3/mu) * 2*pi;
        end
        options_fwd.Sec = T_orb / options.delta;
        [~, xout] = odeMPCI(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
        integration_number = integration_number + 1;
        if ~all(isfinite(xout(end,:))), break, end

        normdr2_prev = norm(dr2_cur);
        dr2_cur = xout(end,1:3) - r2;

        dx = (v1_cur  - v1_prev)';
        df = (dr2_cur - dr2_prev)';
        if norm(dx) > 0
            J = J + (df - J*dx) * dx' / (dx'*dx);
        end
        n_broyden = n_broyden + 1;

        if do_print, fprintf('integration %d (B fwd): normdr2=%g\n', integration_number, norm(dr2_cur)); end
        if build_history
            history.iterations(end+1) = integration_number;
            history.normdr2(end+1)    = norm(dr2_cur);
            history.phase{end+1}      = 'Broyden_fwd';
        end
    end

    v1_B    = v1_cur;
    dr2_B   = xout(end,1:3) - r2;
    dr2_cur = dr2_B;

    if do_print, fprintf('switch to Broyden (fwd): normdr2=%g\n', norm(dr2_cur)); end
end

%% Capture Phase 1 stats
i_phase1 = integration_number - 2;
rt_ph1   = toc(t_ph1);

%% Check if Phase 1 already converged
if norm(dr2_cur) <= Tolr
    if build_history, history.converged = true; history.final_error = norm(dr2_cur); end
    stats.integration_number = integration_number;
    stats.integration_phase1 = i_phase1;
    stats.runtime_phase1     = rt_ph1;
    v1t = v1_cur;
    v2t = xout(end,4:6);
    [hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
    stats.hitearthvar = hitearthvar;
    stats.hitrad      = hitrad;
    return
end

%% ---- Broyden init: refine J with secant from lambertb seed to Phase-1 exit ----
dx0 = (v1_B - v1_A)';
df0 = (dr2_B - dr2_A)';
if norm(dx0) > 0
    J = J + (df0 - J*dx0) * dx0' / (dx0'*dx0);
end

v1_cur = v1_B;

t_ph2 = tic;

%% ---- Phase 2: Broyden forward ----
while norm(dr2_cur) > Tolr && integration_number < maxIter
    v1_next  = v1_cur - (J \ dr2_cur')';
    v1_prev  = v1_cur;
    dr2_prev = dr2_cur;
    v1_cur   = v1_next;

    oe = kep_elements(r1, v1_cur, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options_fwd.Sec = T_orb / options.delta;

    [~, xout] = odeMPCI(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
    integration_number = integration_number + 1;
    dr2_cur = xout(end,1:3) - r2;

    dx = (v1_cur  - v1_prev)';
    df = (dr2_cur - dr2_prev)';
    if norm(dx) > 0
        J = J + (df - J*dx) * dx' / (dx'*dx);
    end

    if do_print, fprintf('integration %d (Broyden): normdr2=%g\n', integration_number, norm(dr2_cur)); end
    if build_history
        history.iterations(end+1) = integration_number;
        history.normdr2(end+1)    = norm(dr2_cur);
        history.phase{end+1}      = 'Broyden';
    end

    if ~all(isreal(v1_cur)) || any(~isfinite(v1_cur)) || ~isfinite(norm(dr2_cur))
        break
    end
end

rt_ph2 = toc(t_ph2);

%% Convergence check
if ~isfinite(norm(dr2_cur)) || norm(dr2_cur) > Tolr
    warning('prtlambertBRB: did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    stats.integration_number = integration_number;
    stats.integration_phase1 = i_phase1;
    stats.integration_phase2 = integration_number - 2 - i_phase1;
    stats.runtime_phase1     = rt_ph1;
    stats.runtime_phase2     = rt_ph2;
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

if build_history, history.converged = true; history.final_error = norm(dr2_cur); end

stats.integration_number = integration_number;
stats.integration_phase1 = i_phase1;
stats.integration_phase2 = integration_number - 2 - i_phase1;
stats.runtime_phase1     = rt_ph1;
stats.runtime_phase2     = rt_ph2;
v1t = v1_cur;
v2t = xout(end,4:6);
[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.hitearthvar = hitearthvar;
stats.hitrad      = hitrad;
end
