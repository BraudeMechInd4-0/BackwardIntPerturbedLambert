function [v1t, v2t, stats, history] = prtlambertT(forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
%PRTLAMBERTT  Perturbed Lambert solver — pure Thompson shooting (forward only).
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertT( ...
%       forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
%
% Inputs:
%   forcemodel   - function handle @(t,x) for the perturbing force model
%   odesolver    - function handle to the integrator, e.g. @odeMPCI
%   r1           - initial position (1x3, km)
%   r2           - target  position (1x3, km)
%   v1           - reference velocity for lambertb seed (1x3, km/s); [0 0 0] works
%   dm           - direction of motion ('L' or 'S')
%   de           - orbital energy    ('L' or 'H')
%   nrev         - number of revolutions
%   dtsec        - transfer time (s)
%   Learning_rate - step-size scaling in [0,1] (use 1 for standard Thompson)
%
% Outputs:
%   v1t     - transfer velocity at r1 (1x3, km/s); NaN if not converged
%   v2t     - transfer velocity at r2 (1x3, km/s); NaN if not converged
%   stats   - struct: integration_number, integration_preamble, integration_phase1,
%             integration_phase2, runtime_preamble, runtime_phase1, runtime_phase2,
%             direction_chosen, hitearthvar, hitrad
%   history - struct with convergence data (populated only when options.debug=true)
%             Fields: iterations, normdr2, converged, final_error

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

if do_print, fprintf('starting T\n'); end
if build_history
    history = struct('iterations', [], 'normdr2', [], 'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

deltar  = inf;
prop_r2 = r2;
i       = 1;

while norm(deltar) >= Tolr && i < maxIter
    [v1t, v2t] = lambertb(r1, prop_r2, v1, dm, de, nrev, dtsec);
    if ~all(isreal(v1t)) || ~all(isreal(v2t)) || all(isnan(v1t)) || all(isnan(v2t))
        if do_print, fprintf('stopped after %d integrations\n', i); end
        if build_history, history.converged = false; history.final_error = NaN; end
        v1t = NaN(1,3);
        v2t = NaN(1,3);
        stats.runtime_phase1 = toc(t_ph1);
        return
    end

    oe = kep_elements(r1, v1t, mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        T_orb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        T_orb = sqrt(a^3/mu) * 2*pi;
    end
    options.Sec = T_orb / options.delta;

    [~, x2] = odesolver(forcemodel, [0 dtsec], [r1 v1t], options);
    if ~all(isfinite(x2(end,:))), break, end
    r_new   = x2(end, 1:3);
    v2t     = x2(end, 4:6);
    deltar  = (r_new - r2) * Learning_rate;

    if do_print, fprintf('iteration %d: error = %g km\n', i, norm(r_new - r2)); end
    if build_history
        history.iterations(end+1) = i;
        history.normdr2(end+1)    = norm(r_new - r2);
    end

    prop_r2 = prop_r2 - deltar;
    i = i + 1;
end

if i >= maxIter || ~isfinite(norm(deltar)) || norm(deltar) > Tolr
    warning('prtlambertT: Thompson did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    v1t = NaN(1,3);
    v2t = NaN(1,3);
    stats.integration_number = i;
    stats.integration_phase1 = i;
    stats.runtime_phase1 = toc(t_ph1);
    return
end

if build_history
    history.converged   = true;
    history.final_error = norm(deltar) / Learning_rate;
end

[hitearthvar, hitrad] = hitearth(100, r1, v1t, r2, v2t, nrev);
stats.integration_number = i;
stats.integration_phase1 = i;
stats.runtime_phase1 = toc(t_ph1);
stats.hitearthvar = hitearthvar;
stats.hitrad = hitrad;
end
