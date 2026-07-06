function [v1t, v2t, stats, history] = prtlambertMPSB(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTMPSB  Perturbed Lambert solver — MPS (forward) then Broyden.
%
% Phase 1: MPS (forward), runs for N_switch_mps = max(1,floor(N_switch/4))
%           full MPS iterations or until normdr2 < T_switch.
% Phase 2: Broyden's quasi-Newton method (forward), initialized from the
%           secant between the lambertb seed and the MPS exit state.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertMPSB( ...
%       forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
%
% Outputs:
%   v1t     - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t     - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   stats   - struct: integration_number, integration_preamble, integration_phase1,
%             integration_phase2, runtime_preamble, runtime_phase1, runtime_phase2,
%             direction_chosen, hitearthvar, hitrad
%   history - struct with convergence data (populated only when options.debug=true)

mu = 398600.5;
r_GEO = 42164;
r_max_allowed = 1.5 * r_GEO;

if nargin < 10,  options = struct();  end
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

if do_print, fprintf('starting MPSB\n'); end
if build_history
    history = struct('iterations', [], 'integration_number', [], 'normdr2', [], ...
                     'phase', {{}}, 'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

%% ---- Initial lambertb + forward propagation ----
[v1ref, ~] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1ref)) || any(isnan(v1ref)) || any(~isfinite(v1ref))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_phase1 = toc(t_ph1);
    return
end

xref  = [r1 v1ref];
oe    = kep_elements(r1, v1ref, mu);
a_ref = oe(1);
if a_ref <= 0 || ~isreal(a_ref) || isnan(a_ref) || a_ref > 2*r_max_allowed
    Torbref = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    Torbref = sqrt(a_ref^3/mu) * 2*pi;
end
options.Sec = Torbref / options.delta;

[~, xnew] = odesolver(forcemodel, [0 dtsec], xref, options);
integration_number = 1;
rdep = r2 - xnew(end,1:3);

v1_A  = v1ref;
dr2_A = -rdep;

if do_print, fprintf('integration %d (MPS init): normdr2=%g\n', integration_number, norm(rdep)); end
if build_history
    history.iterations(end+1)         = 0;
    history.integration_number(end+1) = integration_number;
    history.normdr2(end+1)            = norm(rdep);
    history.phase{end+1}              = 'MPS_init';
end

delta_size  = 1e-7;
xpert       = zeros(3,6);
n_mps          = 0;
J_mps          = eye(3);
rdep_prev_norm = norm(rdep);
xref_v_prev    = xref(4:6);
rdep_prev      = rdep;

%% ---- Phase 1: MPS forward ----
while norm(rdep) > T_switch && (n_mps < N_switch_mps || norm(rdep) > rdep_prev_norm) && integration_number + 3 < maxIter
    if n_mps == 0
        % Iteration 1: full FD Jacobian
        normv0    = norm(xref(4:6));
        deltardot = [zeros(3), eye(3)*delta_size*normv0];
        xprtrbd   = deltardot + xref;

        for j = 1:3
            [~, xtmp]  = odesolver(forcemodel, [0 dtsec], xprtrbd(j,:), options);
            xpert(j,:) = xtmp(end,:);
        end
        integration_number = integration_number + 3;
        if ~all(isfinite(xpert(:))), break, end

        DR    = (xpert(:,1:3) - xnew(end,1:3))';
        J_mps = DR / (delta_size * normv0);

        xref_v_prev = xref(4:6);
        rdep_prev   = rdep;
        alpha = (DR \ rdep')';
        xref  = xref + alpha * deltardot;
    else
        % Iterations 2+: rank-1 Broyden update (3 FD integrations saved per iteration)
        dx = (xref(4:6) - xref_v_prev)';
        df = (rdep_prev - rdep)';
        if norm(dx) > 0
            J_mps = J_mps + (df - J_mps*dx) * dx' / (dx'*dx);
        end
        xref_v_prev = xref(4:6);
        rdep_prev   = rdep;
        dv = (J_mps \ rdep')';
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

    [~, xnew] = odesolver(forcemodel, [0 dtsec], xref, options);
    integration_number = integration_number + 1;
    if ~all(isfinite(xnew(end,:))), break, end
    rdep_prev_norm = norm(rdep);
    rdep   = r2 - xnew(end,1:3);
    n_mps  = n_mps + 1;

    if do_print, fprintf('integration %d (MPS): normdr2=%g\n', integration_number, norm(rdep)); end
    if build_history
        history.iterations(end+1)         = n_mps;
        history.integration_number(end+1) = integration_number;
        history.normdr2(end+1)            = norm(rdep);
        history.phase{end+1}              = 'MPS';
    end
end

%% Capture Phase 1 stats
i_phase1 = integration_number;
rt_ph1   = toc(t_ph1);

v1_B  = xnew(1,4:6);
dr2_B = -rdep;

if norm(dr2_B) <= Tolr
    if build_history, history.converged = true; history.final_error = norm(dr2_B); end
    stats.integration_number = i_phase1;
    stats.integration_phase1 = i_phase1;
    stats.runtime_phase1     = rt_ph1;
    v1t = v1_B;
    v2t = xnew(end,4:6);
    [hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
    stats.hitearthvar = hitearthvar;
    stats.hitrad      = hitrad;
    return
end

if do_print, fprintf('switch to Broyden: normdr2=%g\n', norm(dr2_B)); end

%% ---- Broyden initialization (J from last MPS FD Jacobian + secant refinement) ----
J   = J_mps;
dx0 = (v1_B - v1_A)';
df0 = (dr2_B - dr2_A)';
if norm(dx0) > 0
    J = J + (df0 - J*dx0) * dx0' / (dx0'*dx0);
end

v1_cur  = v1_B;
dr2_cur = dr2_B;

t_ph2 = tic;

%% ---- Phase 2: Broyden forward ----
options_fwd = options;
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

    [~, xout]  = odesolver(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
    integration_number = integration_number + 1;
    dr2_cur = xout(end,1:3) - r2;

    dx = (v1_cur  - v1_prev)';
    df = (dr2_cur - dr2_prev)';
    if norm(dx) > 0
        J = J + (df - J*dx) * dx' / (dx'*dx);
    end

    if do_print, fprintf('integration %d (Broyden): normdr2=%g\n', integration_number, norm(dr2_cur)); end
    if build_history
        history.iterations(end+1)         = n_mps + (integration_number - n_mps*4);
        history.integration_number(end+1) = integration_number;
        history.normdr2(end+1)            = norm(dr2_cur);
        history.phase{end+1}              = 'Broyden';
    end

    if ~all(isreal(v1_cur)) || any(~isfinite(v1_cur)) || ~isfinite(norm(dr2_cur))
        break
    end
end

rt_ph2 = toc(t_ph2);

if ~isfinite(norm(dr2_cur)) || norm(dr2_cur) > Tolr
    warning('prtlambertMPSB: did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    stats.integration_number = integration_number;
    stats.integration_phase1 = i_phase1;
    stats.integration_phase2 = integration_number - i_phase1;
    stats.runtime_phase1     = rt_ph1;
    stats.runtime_phase2     = rt_ph2;
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

if build_history, history.converged = true; history.final_error = norm(dr2_cur); end

stats.integration_number = integration_number;
stats.integration_phase1 = i_phase1;
stats.integration_phase2 = integration_number - i_phase1;
stats.runtime_phase1     = rt_ph1;
stats.runtime_phase2     = rt_ph2;
v1t = v1_cur;
v2t = xout(end,4:6);
[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.hitearthvar = hitearthvar;
stats.hitrad      = hitrad;
end
