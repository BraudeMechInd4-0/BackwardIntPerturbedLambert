function [v1t, v2t, stats, history] = prtlambertTB(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
%PRTLAMBERTTB  Perturbed Lambert solver — Thompson (forward) then Broyden.
%
% Phase 1: Thompson shooting (forward), runs for N_switch iterations or
%           until normdr2 < T_switch (whichever comes first).
% Phase 2: Broyden's quasi-Newton method, initialized from the two-point
%           secant between the lambertb seed and the Thompson exit state.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertTB( ...
%       forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
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

if nargin < 10
    Learning_rate = 1;
end
if nargin < 11,  options = struct();  end
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

if do_print, fprintf('starting TB\n'); end
if build_history
    history = struct('iterations', [], 'normdr2', [], 'phase', {{}}, ...
                     'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

%% ---- Initial lambertb + forward propagation ----
[v1tnew, ~] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1tnew)) || any(isnan(v1tnew)) || any(~isfinite(v1tnew))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_phase1 = toc(t_ph1);
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

[~, xout] = odesolver(forcemodel, [0 dtsec], [r1 v1tnew], options_fwd);
integration_number = 1;
deltar2  = (xout(end,1:3) - r2) * Learning_rate;
normdr2  = norm(deltar2);

dr2_A  = xout(end,1:3) - r2;

if do_print, fprintf('integration %d (T init): normdr2=%g\n', integration_number, normdr2/Learning_rate); end
if build_history
    history.iterations(end+1) = integration_number;
    history.normdr2(end+1)    = normdr2 / Learning_rate;
    history.phase{end+1}      = 'Thompson';
end

prop_r2    = r2;
v1_cur     = v1tnew;
J          = eye(3);
v1_t_prev  = v1_cur;
dr2_t_prev = dr2_A;

%% ---- Phase 1: Thompson forward ----
n_thompson   = 0;
normdr2_prev = normdr2;
while normdr2 > T_switch && (n_thompson < N_switch || normdr2 > normdr2_prev) && integration_number < maxIter
    prop_r2  = prop_r2 - deltar2;
    [v1_cur, ~] = lambertb(r1, prop_r2, v1_cur, dm, de, nrev, dtsec);
    if ~all(isreal(v1_cur)) || any(isnan(v1_cur))
        if build_history, history.converged = false; end
        rt_ph1 = toc(t_ph1);
        stats.integration_number = integration_number;
        stats.integration_phase1 = integration_number;
        stats.runtime_phase1     = rt_ph1;
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

    [~, xout]   = odesolver(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
    integration_number = integration_number + 1;
    if ~all(isfinite(xout(end,:))), break, end
    normdr2_prev = normdr2;
    deltar2  = (xout(end,1:3) - r2) * Learning_rate;
    normdr2  = norm(deltar2);
    dr2_t_cur = xout(end,1:3) - r2;
    dx_t = (v1_cur - v1_t_prev)';
    df_t = (dr2_t_cur - dr2_t_prev)';
    if norm(dx_t) > 0
        J = J + (df_t - J*dx_t) * dx_t' / (dx_t'*dx_t);
    end
    v1_t_prev  = v1_cur;
    dr2_t_prev = dr2_t_cur;
    n_thompson = n_thompson + 1;

    if do_print, fprintf('integration %d (Thompson): normdr2=%g\n', integration_number, normdr2/Learning_rate); end
    if build_history
        history.iterations(end+1) = integration_number;
        history.normdr2(end+1)    = normdr2 / Learning_rate;
        history.phase{end+1}      = 'Thompson';
    end
end

%% Capture Phase 1 stats
i_phase1 = integration_number;
rt_ph1   = toc(t_ph1);

%% Check if Thompson already converged
if normdr2 <= Tolr
    if build_history, history.converged = true; history.final_error = normdr2 / Learning_rate; end
    stats.integration_number = i_phase1;
    stats.integration_phase1 = i_phase1;
    stats.runtime_phase1     = rt_ph1;
    v1t = v1_cur;
    v2t = xout(end, 4:6);
    [hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
    stats.hitearthvar = hitearthvar;
    stats.hitrad      = hitrad;
    return
end

%% ---- Broyden initialization (J already conditioned from Thompson) ----
v1_B    = v1_cur;
dr2_B   = xout(end,1:3) - r2;
dr2_cur = dr2_B;
v1_cur  = v1_B;

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
    warning('prtlambertTB: did not converge');
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
v2t = xout(end, 4:6);
[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.hitearthvar = hitearthvar;
stats.hitrad      = hitrad;
end
