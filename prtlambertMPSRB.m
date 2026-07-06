function [v1t, v2t, stats, history] = prtlambertMPSRB(forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTMPSRB  Perturbed Lambert solver — MPS (select fwd/bwd) then Broyden.
%
% Preamble: propagate forward from r1 and backward from r2; choose the
%           direction with the smaller initial residual.
% Phase 1:  MPS in the chosen direction for N_switch_mps full iterations
%           or until normdr < T_switch.
% Phase 2:  Broyden's quasi-Newton method (always forward), initialized
%           from the secant between the lambertb seed and the MPS exit state.
%
% History always tracks the FORWARD error (normdr2), even during backward phase.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertMPSRB( ...
%       forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%
% Outputs:
%   v1t     - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t     - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   stats   - struct: integration_number, integration_preamble, integration_phase1,
%             integration_phase2, runtime_preamble, runtime_phase1, runtime_phase2,
%             direction_chosen ('forward'|'backward'), hitearthvar, hitrad
%   history - struct with convergence data (populated only when options.debug=true)

mu = 398600.5;
r_GEO = 42164;
r_max_allowed = 1.5 * r_GEO;

if nargin < 9,   options = struct();  end
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
N_switch_mps = max(1, floor(N_switch / 4));

if do_print, fprintf('starting MPSRB\n'); end
if build_history
    history = struct('iterations', [], 'integration_number', [], 'normdr2', [], ...
                     'phase', {{}}, 'direction_chosen', '', ...
                     'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_pre = tic;

%% ---- Preamble: lambertb + forward + backward propagation ----
[v1ref, v2ref] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1ref)) || any(isnan(v1ref)) || ~all(isreal(v2ref)) || any(isnan(v2ref))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_preamble = toc(t_pre);
    return
end

x1ref = [r1 v1ref];
oe    = kep_elements(r1, v1ref, mu);
a     = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    Torb_fwd = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    Torb_fwd = sqrt(a^3/mu) * 2*pi;
end
options_fwd     = options;
options_fwd.Sec = Torb_fwd / options.delta;
[~, xfwd_pre]   = odeMPCI(forcemodel, [0 dtsec], x1ref, options_fwd);

x2ref = [r2 v2ref];
oe    = kep_elements(r2, v2ref, mu);
a     = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    Torb_bkw = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    Torb_bkw = sqrt(a^3/mu) * 2*pi;
end
options_bkw     = options;
options_bkw.Sec = Torb_bkw / options.delta;
[~, xbkw_pre]   = odeMPCIrev(forcemodel, [0 dtsec], x2ref, options_bkw);

integration_number = 2;
normdr1 = norm(xbkw_pre(1,1:3) - r1);
normdr2 = norm(xfwd_pre(end,1:3) - r2);

v1_A  = v1ref;
dr2_A = xfwd_pre(end,1:3) - r2;

rt_pre = toc(t_pre);
stats.integration_preamble = 2;
stats.runtime_preamble     = rt_pre;

if do_print, fprintf('init (MPSRB): normdr1=%g, normdr2=%g\n', normdr1, normdr2); end
if build_history
    history.iterations(end+1)         = 0;
    history.integration_number(end+1) = integration_number;
    history.normdr2(end+1)            = normdr2;
    history.phase{end+1}              = 'init';
end

use_backward = isfinite(normdr1) && (normdr1 < normdr2);
if use_backward
    stats.direction_chosen = 'backward';
else
    stats.direction_chosen = 'forward';
end
if do_print, fprintf('MPSRB: using %s\n', stats.direction_chosen); end
if build_history, history.direction_chosen = stats.direction_chosen; end

delta_size = 1e-7;
xpert      = zeros(3,6);
n_mps      = 0;
J_phase1   = eye(3);

t_ph1 = tic;

%% ====================================================================
%% BACKWARD BRANCH
%% ====================================================================
if use_backward
    xref         = x2ref;
    xbkw_cur     = xbkw_pre;
    rdep_bkw           = r1 - xbkw_cur(1,1:3);
    rdep_bkw_prev_norm = norm(rdep_bkw);
    xref_v_prev        = xref(4:6);
    rdep_bkw_prev      = rdep_bkw;

    while norm(rdep_bkw) > T_switch && (n_mps < N_switch_mps || norm(rdep_bkw) > rdep_bkw_prev_norm) && integration_number + 3 < maxIter
        if n_mps == 0
            % Iteration 1: full FD Jacobian
            normv0    = norm(xref(4:6));
            deltardot = [zeros(3), eye(3)*delta_size*normv0];
            xprtrbd   = deltardot + xref;

            for j = 1:3
                [~, xtmp]  = odeMPCIrev(forcemodel, [0 dtsec], xprtrbd(j,:), options_bkw);
                xpert(j,:) = xtmp(1,:);
            end
            integration_number = integration_number + 3;
            if ~all(isfinite(xpert(:))), break, end

            DR       = (xpert(:,1:3) - xbkw_cur(1,1:3))';
            J_phase1 = DR / (delta_size * normv0);

            xref_v_prev   = xref(4:6);
            rdep_bkw_prev = rdep_bkw;
            alpha    = (DR \ rdep_bkw')';
            xref     = xref + alpha * deltardot;
        else
            % Iterations 2+: rank-1 Broyden update
            dx = (xref(4:6) - xref_v_prev)';
            df = (rdep_bkw_prev - rdep_bkw)';
            if norm(dx) > 0
                J_phase1 = J_phase1 + (df - J_phase1*dx) * dx' / (dx'*dx);
            end
            xref_v_prev   = xref(4:6);
            rdep_bkw_prev = rdep_bkw;
            dv   = (J_phase1 \ rdep_bkw')';
            xref(4:6) = xref(4:6) + dv;
        end
        oe = kep_elements(xref(1:3), xref(4:6), mu);
        a  = oe(1);
        if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
            Torb_bkw = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            Torb_bkw = sqrt(a^3/mu) * 2*pi;
        end
        options_bkw.Sec = Torb_bkw / options.delta;

        [~, xbkw_cur] = odeMPCIrev(forcemodel, [0 dtsec], xref, options_bkw);
        integration_number = integration_number + 1;
        if ~all(isfinite(xbkw_cur(1,:))), break, end
        rdep_bkw_prev_norm = norm(rdep_bkw);
        rdep_bkw = r1 - xbkw_cur(1,1:3);
        n_mps    = n_mps + 1;

        if build_history
            v1_check = xbkw_cur(1,4:6);
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
            history.iterations(end+1)         = n_mps;
            history.integration_number(end+1) = integration_number;
            history.normdr2(end+1)            = normdr2_check;
            history.phase{end+1}              = 'MPS_bkw';
            if do_print, fprintf('integration %d (MPS bkw): normdr1=%g, normdr2=%g\n', ...
                    integration_number, norm(rdep_bkw), normdr2_check); end
        elseif do_print
            fprintf('iteration %d (MPS bkw): normdr1=%g\n', n_mps, norm(rdep_bkw));
        end
    end

    J_phase1 = -J_phase1';
    v1_B = xbkw_cur(1,4:6);
    oe   = kep_elements(r1, v1_B, mu);
    a    = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options_fwd.Sec = T_orb / options.delta;
    [~, xout_B]    = odeMPCI(forcemodel, [0 dtsec], [r1 v1_B], options_fwd);
    integration_number = integration_number + 1;
    dr2_B   = xout_B(end,1:3) - r2;
    xout    = xout_B;

%% ====================================================================
%% FORWARD BRANCH
%% ====================================================================
else
    xref      = x1ref;
    xnew      = xfwd_pre;
    rdep           = r2 - xnew(end,1:3);
    rdep_prev_norm = norm(rdep);
    xref_v_prev    = xref(4:6);
    rdep_prev      = rdep;
    oe    = kep_elements(r1, v1ref, mu);
    a_ref = oe(1);
    if a_ref <= 0 || ~isreal(a_ref) || isnan(a_ref) || a_ref > 2*r_max_allowed
        Torbref = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        Torbref = sqrt(a_ref^3/mu) * 2*pi;
    end
    options.Sec = Torbref / options.delta;

    while norm(rdep) > T_switch && (n_mps < N_switch_mps || norm(rdep) > rdep_prev_norm) && integration_number + 3 < maxIter
        if n_mps == 0
            % Iteration 1: full FD Jacobian
            normv0    = norm(xref(4:6));
            deltardot = [zeros(3), eye(3)*delta_size*normv0];
            xprtrbd   = deltardot + xref;

            for j = 1:3
                [~, xtmp]  = odeMPCI(forcemodel, [0 dtsec], xprtrbd(j,:), options);
                xpert(j,:) = xtmp(end,:);
            end
            integration_number = integration_number + 3;
            if ~all(isfinite(xpert(:))), break, end

            DR       = (xpert(:,1:3) - xnew(end,1:3))';
            J_phase1 = DR / (delta_size * normv0);

            xref_v_prev = xref(4:6);
            rdep_prev   = rdep;
            alpha = (DR \ rdep')';
            xref  = xref + alpha * deltardot;
        else
            % Iterations 2+: rank-1 Broyden update
            dx = (xref(4:6) - xref_v_prev)';
            df = (rdep_prev - rdep)';
            if norm(dx) > 0
                J_phase1 = J_phase1 + (df - J_phase1*dx) * dx' / (dx'*dx);
            end
            xref_v_prev = xref(4:6);
            rdep_prev   = rdep;
            dv = (J_phase1 \ rdep')';
            xref(4:6)   = xref(4:6) + dv;
        end
        oe    = kep_elements(xref(1:3), xref(4:6), mu);
        a_ref = oe(1);
        if a_ref <= 0 || ~isreal(a_ref) || isnan(a_ref) || a_ref > 2*r_max_allowed
            Torbref = 2*pi*sqrt((2*r_max_allowed)^3/mu);
        else
            Torbref = sqrt(a_ref^3/mu) * 2*pi;
        end
        options.Sec = Torbref / options.delta;

        [~, xnew] = odeMPCI(forcemodel, [0 dtsec], xref, options);
        integration_number = integration_number + 1;
        if ~all(isfinite(xnew(end,:))), break, end
        rdep_prev_norm = norm(rdep);
        rdep  = r2 - xnew(end,1:3);
        n_mps = n_mps + 1;

        if do_print, fprintf('integration %d (MPS fwd): normdr2=%g\n', integration_number, norm(rdep)); end
        if build_history
            history.iterations(end+1)         = n_mps;
            history.integration_number(end+1) = integration_number;
            history.normdr2(end+1)            = norm(rdep);
            history.phase{end+1}              = 'MPS_fwd';
        end
    end

    v1_B  = xnew(1,4:6);
    dr2_B = -rdep;
    xout  = xnew;
    options_fwd = options;
end

%% Capture Phase 1 stats
i_phase1 = integration_number - 2;
rt_ph1   = toc(t_ph1);

%% Check if MPS already converged
if norm(dr2_B) <= Tolr
    if build_history
        history.converged   = true;
        history.final_error = norm(dr2_B);
    end
    stats.integration_number = integration_number;
    stats.integration_phase1 = i_phase1;
    stats.runtime_phase1     = rt_ph1;
    v1t = v1_B;
    v2t = xout(end,4:6);
    [hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
    stats.hitearthvar = hitearthvar;
    stats.hitrad      = hitrad;
    return
end

if do_print, fprintf('switch to Broyden: normdr2=%g\n', norm(dr2_B)); end

%% ---- Broyden initialization (J from Phase-1 FD Jacobian + secant refinement) ----
J   = J_phase1;
dx0 = (v1_B - v1_A)';
df0 = (dr2_B - dr2_A)';
if norm(dx0) > 0
    J = J + (df0 - J*dx0) * dx0' / (dx0'*dx0);
end

v1_cur  = v1_B;
dr2_cur = dr2_B;

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

    [~, xout]  = odeMPCI(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
    integration_number = integration_number + 1;
    dr2_cur = xout(end,1:3) - r2;

    dx = (v1_cur  - v1_prev)';
    df = (dr2_cur - dr2_prev)';
    if norm(dx) > 0
        J = J + (df - J*dx) * dx' / (dx'*dx);
    end

    if build_history
        history.iterations(end+1)         = n_mps + (integration_number - n_mps*4);
        history.integration_number(end+1) = integration_number;
        history.normdr2(end+1)            = norm(dr2_cur);
        history.phase{end+1}              = 'Broyden';
    end
    if do_print, fprintf('integration %d (Broyden): normdr2=%g\n', integration_number, norm(dr2_cur)); end

    if ~all(isreal(v1_cur)) || any(~isfinite(v1_cur)) || ~isfinite(norm(dr2_cur))
        break
    end
end

rt_ph2 = toc(t_ph2);

if ~isfinite(norm(dr2_cur)) || norm(dr2_cur) > Tolr
    warning('prtlambertMPSRB: did not converge');
    if build_history
        history.converged   = false;
        history.final_error = NaN;
    end
    stats.integration_number = integration_number;
    stats.integration_phase1 = i_phase1;
    stats.integration_phase2 = integration_number - 2 - i_phase1;
    stats.runtime_phase1     = rt_ph1;
    stats.runtime_phase2     = rt_ph2;
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

if build_history
    history.converged   = true;
    history.final_error = norm(dr2_cur);
end

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
