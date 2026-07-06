function [v1t, v2t, stats, history] = prtlambertBR(forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTBR  Perturbed Lambert solver — Broyden's method, always backward (Rear).
%
% Propagates backward from r2 using odeMPCIrev, applying Broyden's quasi-Newton
% method to adjust v2 until the backward trajectory arrives at r1.
% No direction selection — always integrates backward. No secondary phase.
% J0 initialized via 3 finite-difference perturbations of v2.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertBR( ...
%       forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%
% Outputs:
%   v1t     - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t     - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   stats   - struct: integration_number, integration_preamble, integration_phase1,
%             integration_phase2, runtime_preamble, runtime_phase1, runtime_phase2,
%             direction_chosen, hitearthvar, hitrad
%   history - struct (debug mode only): normdr2 is always the FORWARD error

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

if do_print, fprintf('starting BR\n'); end
if build_history
    history = struct('iterations', [], 'normdr2', [], 'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

%% Preamble: unperturbed Lambert solution
[~, v2tnew] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v2tnew)) || any(isnan(v2tnew)) || any(~isfinite(v2tnew))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_phase1 = toc(t_ph1);
    return
end

oe = kep_elements(r2, v2tnew, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options_bkw     = options;
options_bkw.Sec = T_orb / options.delta;

[~, xoutbkw] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2tnew], options_bkw);
integration_number = 1;
v2_cur  = v2tnew;
dr1_cur = xoutbkw(1,1:3) - r1;

if build_history
    oe_c = kep_elements(r1, xoutbkw(1,4:6), mu);
    a_c  = oe_c(1);
    if a_c <= 0 || ~isreal(a_c) || isnan(a_c) || a_c > 2*r_max_allowed
        T_c = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_c = sqrt(a_c^3/mu) * 2*pi;
    end
    opt_c = options;  opt_c.Sec = T_c / options.delta;
    [~, xfwd_c] = odeMPCI(forcemodel, [0 dtsec], [r1 xoutbkw(1,4:6)], opt_c);
    integration_number = integration_number + 1;
    history.iterations(end+1) = integration_number;
    history.normdr2(end+1)    = norm(xfwd_c(end,1:3) - r2);
    if do_print, fprintf('init (BR): normdr1=%g, normdr2=%g\n', norm(dr1_cur), norm(xfwd_c(end,1:3)-r2)); end
elseif do_print
    fprintf('init (BR): normdr1=%g\n', norm(dr1_cur));
end

%% J0 via finite differences — 3 perturbations of v2 through odeMPCIrev
delta_size = 1e-7;
normv2     = norm(v2_cur);
J          = zeros(3);
for jj = 1:3
    v_pert     = v2_cur;
    v_pert(jj) = v_pert(jj) + delta_size * normv2;
    oe_p = kep_elements(r2, v_pert, mu);  a_p = oe_p(1);
    if a_p <= 0 || ~isreal(a_p) || isnan(a_p) || a_p > 2*r_max_allowed
        T_p = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_p = sqrt(a_p^3/mu) * 2*pi;
    end
    opt_p = options_bkw;  opt_p.Sec = T_p / options.delta;
    [~, x_pert] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v_pert], opt_p);
    J(:,jj) = (x_pert(1,1:3) - xoutbkw(1,1:3))' / (delta_size * normv2);
end
integration_number = integration_number + 3;
if ~all(isfinite(J(:))),  J = eye(3);  end

%% Broyden loop: adjust v2 to minimise dr1 = xoutbkw(1,1:3) - r1
while norm(dr1_cur) > Tolr && integration_number < maxIter
    v2_prev  = v2_cur;
    dr1_prev = dr1_cur;
    v2_cur   = v2_cur - (J \ dr1_cur')';

    oe = kep_elements(r2, v2_cur, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options_bkw.Sec = T_orb / options.delta;

    [~, xoutbkw] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2_cur], options_bkw);
    integration_number = integration_number + 1;
    if ~all(isfinite(xoutbkw(1,:))), break, end
    dr1_cur = xoutbkw(1,1:3) - r1;

    dx = (v2_cur  - v2_prev)';
    df = (dr1_cur - dr1_prev)';
    if norm(dx) > 0
        J = J + (df - J*dx) * dx' / (dx'*dx);
    end

    if build_history
        v1_check = xoutbkw(1,4:6);
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
        history.iterations(end+1) = integration_number;
        history.normdr2(end+1)    = norm(xfwd_c(end,1:3) - r2);
        if do_print, fprintf('integration %d (BR): normdr1=%g, normdr2=%g\n', ...
                integration_number, norm(dr1_cur), norm(xfwd_c(end,1:3)-r2)); end
    elseif do_print
        fprintf('integration %d (BR): normdr1=%g\n', integration_number, norm(dr1_cur));
    end

    if ~all(isreal(v2_cur)) || any(~isfinite(v2_cur)) || ~isfinite(norm(dr1_cur))
        break
    end
end

%% Convergence check on backward residual
if norm(dr1_cur) > Tolr || ~isfinite(norm(dr1_cur))
    warning('prtlambertBR: Broyden did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    stats.integration_number = integration_number;
    stats.integration_phase1 = integration_number;
    stats.runtime_phase1     = toc(t_ph1);
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

%% Extract v1 from backward trajectory at t=0
v1t = xoutbkw(1, 4:6);

%% Final forward propagation to get v2t (and verify forward error for history)
oe = kep_elements(r1, v1t, mu);
a  = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    T_orb = sqrt(a^3/mu) * 2*pi;
end
options.Sec = T_orb / options.delta;
[~, xfwd] = odeMPCI(forcemodel, [0 dtsec], [r1 v1t], options);
integration_number = integration_number + 1;
v2t = xfwd(end, 4:6);

if build_history
    history.converged   = true;
    history.final_error = norm(xfwd(end,1:3) - r2);
end

stats.integration_number = integration_number;
stats.integration_phase1 = integration_number;
stats.runtime_phase1     = toc(t_ph1);
[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.hitearthvar = hitearthvar;
stats.hitrad      = hitrad;
end
