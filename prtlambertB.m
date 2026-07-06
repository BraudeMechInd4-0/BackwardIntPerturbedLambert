function [v1t, v2t, stats, history] = prtlambertB(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTB  Perturbed Lambert solver — pure Broyden's method.
%
% Uses the unperturbed lambertb solution as the initial guess, then applies
% Broyden's quasi-Newton method (rank-1 Jacobian update) with J0 = eye(3).
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertB( ...
%       forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
%
% Outputs:
%   v1t     - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t     - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   stats   - struct: integration_number, integration_preamble, integration_phase1,
%             integration_phase2, runtime_preamble, runtime_phase1, runtime_phase2,
%             direction_chosen, hitearthvar, hitrad
%   history - struct (populated only when options.debug=true)

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

if do_print, fprintf('starting B\n'); end
if build_history
    history = struct('iterations', [], 'normdr2', [], 'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

%% Initial guess from unperturbed Lambert
[v1_cur, ~] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1_cur)) || any(isnan(v1_cur)) || any(~isfinite(v1_cur))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_phase1 = toc(t_ph1);
    return
end

%% First propagation
oe = kep_elements(r1, v1_cur, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options_fwd     = options;
options_fwd.Sec = T_orb / options.delta;

[~, xout] = odesolver(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
integration_number = 1;
dr2_cur = xout(end, 1:3) - r2;

if do_print, fprintf('integration %d (init): normdr2=%g\n', integration_number, norm(dr2_cur)); end
if build_history
    history.iterations(end+1) = integration_number;
    history.normdr2(end+1)    = norm(dr2_cur);
end

%% Broyden's method — J0 via finite differences (3 perturbations)
delta_size = 1e-7;
normv0     = norm(v1_cur);
J          = zeros(3);
for jj = 1:3
    v_pert    = v1_cur;
    v_pert(jj) = v_pert(jj) + delta_size * normv0;
    oe_p = kep_elements(r1, v_pert, mu);  a_p = oe_p(1);
    if a_p <= 0 || ~isreal(a_p) || isnan(a_p) || a_p > 2*r_max_allowed
        T_p = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_p = sqrt(a_p^3/mu) * 2*pi;
    end
    opt_p = options_fwd;  opt_p.Sec = T_p / options.delta;
    [~, x_pert] = odesolver(forcemodel, [0 dtsec], [r1 v_pert], opt_p);
    J(:,jj) = (x_pert(end,1:3) - xout(end,1:3))' / (delta_size * normv0);
end
integration_number = integration_number + 3;
if ~all(isfinite(J(:))),  J = eye(3);  end

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

    [~, xout] = odesolver(forcemodel, [0 dtsec], [r1 v1_cur], options_fwd);
    integration_number = integration_number + 1;
    dr2_cur = xout(end, 1:3) - r2;

    dx = (v1_cur  - v1_prev)';
    df = (dr2_cur - dr2_prev)';
    if norm(dx) > 0
        J = J + (df - J*dx) * dx' / (dx'*dx);
    end

    if do_print, fprintf('integration %d (Broyden): normdr2=%g\n', integration_number, norm(dr2_cur)); end
    if build_history
        history.iterations(end+1) = integration_number;
        history.normdr2(end+1)    = norm(dr2_cur);
    end

    if ~all(isreal(v1_cur)) || any(isnan(v1_cur)) || any(~isfinite(v1_cur)) || ...
       ~isfinite(norm(dr2_cur))
        break
    end
end

%% Check convergence
if ~isfinite(norm(dr2_cur)) || norm(dr2_cur) > Tolr
    warning('prtlambertB: Broyden did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    stats.integration_number = integration_number;
    stats.integration_phase1 = integration_number;
    stats.runtime_phase1     = toc(t_ph1);
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

if build_history, history.converged = true; history.final_error = norm(dr2_cur); end

stats.integration_number = integration_number;
stats.integration_phase1 = integration_number;
stats.runtime_phase1     = toc(t_ph1);
v1t = v1_cur;
v2t = xout(end, 4:6);
[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.hitearthvar = hitearthvar;
stats.hitrad      = hitrad;
end
