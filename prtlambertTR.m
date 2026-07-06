function [v1t, v2t, stats, history] = prtlambertTR(forcemodel, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
%PRTLAMBERTTR  Perturbed Lambert solver — Thompson shooting, always backward (Rear).
%
% Propagates backward from r2 using odeMPCIrev, applying Thompson's
% correction to the departure point until the backward trajectory lands at r1.
% No direction selection — always integrates backward. No Broyden phase.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertTR( ...
%       forcemodel, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
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

if nargin < 9
    Learning_rate = 1;
end
if nargin < 10,  options = struct();  end
if ~isfield(options, 'Tolr'),    options.Tolr    = 1e-8;  end
if ~isfield(options, 'maxIter'), options.maxIter = 100;   end
if ~isfield(options, 'delta'),   options.delta   = 16;    end

Tolr    = options.Tolr;
maxIter = options.maxIter;
do_print      = isfield(options, 'debug') && options.debug;
build_history = nargout >= 4;

if do_print, fprintf('starting TR\n'); end
if build_history
    history = struct('iterations', [], 'normdr2', [], 'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

%% Preamble: unperturbed Lambert solution
[v1tnew, v2tnew] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1tnew)) || any(isnan(v1tnew)) || ~all(isreal(v2tnew)) || any(isnan(v2tnew))
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
options.Sec = T_orb / options.delta;

[~, xoutbkw]  = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2tnew], options);
integration_number = 1;
prop_r1  = r1;
deltar1  = (xoutbkw(1,1:3) - r1) * Learning_rate;

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
    normdr2_init = norm(xfwd_c(end,1:3) - r2);
    history.iterations(end+1) = 0;
    history.normdr2(end+1)    = normdr2_init;
    if do_print, fprintf('init (TR): normdr1=%g, normdr2=%g\n', norm(deltar1), normdr2_init); end
elseif do_print
    fprintf('init (TR): normdr1=%g\n', norm(deltar1));
end

i = 1;

%% Backward Thompson loop
while norm(deltar1) > Tolr && i < maxIter
    prop_r1  = prop_r1 - deltar1;
    v1tnew   = xoutbkw(1,4:6);
    [v1tnew, v2tnew] = lambertb(prop_r1, r2, v1tnew, dm, de, nrev, dtsec);
    if ~all(isreal(v1tnew)) || any(isnan(v1tnew)) || ~all(isreal(v2tnew)) || any(isnan(v2tnew))
        if do_print, fprintf('TR: lambertb failed at iteration %d\n', i); end
        if build_history, history.converged = false; end
        stats.integration_number = integration_number;
        stats.integration_phase1 = integration_number;
        stats.runtime_phase1     = toc(t_ph1);
        v1t = NaN(1,3);  v2t = NaN(1,3);
        return
    end

    oe = kep_elements(r2, v2tnew, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options.Sec = T_orb / options.delta;

    [~, xoutbkw] = odeMPCIrev(forcemodel, [0 dtsec], [r2 v2tnew], options);
    integration_number = integration_number + 1;
    if ~all(isfinite(xoutbkw(1,:))), break, end
    deltar1 = (xoutbkw(1,1:3) - r1) * Learning_rate;

    if build_history
        v1_check = xoutbkw(1, 4:6);
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
        history.iterations(end+1) = i;
        history.normdr2(end+1)    = normdr2_check;
        if do_print, fprintf('iteration %d (TR): normdr1=%g, normdr2=%g\n', i, norm(deltar1), normdr2_check); end
    elseif do_print
        fprintf('iteration %d (TR): normdr1=%g\n', i, norm(deltar1));
    end

    i = i + 1;
end

%% Convergence check on backward residual
if norm(deltar1) > Tolr || ~isfinite(norm(deltar1))
    warning('prtlambertTR: did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    stats.integration_number = integration_number;
    stats.integration_phase1 = integration_number;
    stats.runtime_phase1     = toc(t_ph1);
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

%% Extract outputs: v1 is the velocity at r1 from the backward trajectory
v1t = xoutbkw(1, 4:6);

%% Final forward propagation to get v2t
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

if build_history, history.converged = true; history.final_error = norm(xfwd(end,1:3) - r2); end

stats.integration_number = integration_number;
stats.integration_phase1 = integration_number;
stats.runtime_phase1     = toc(t_ph1);
[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.hitearthvar = hitearthvar;
stats.hitrad      = hitrad;
end
