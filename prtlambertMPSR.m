function [v1t, v2t, stats, history] = prtlambertMPSR(forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
%PRTLAMBERTMPSR  Perturbed Lambert solver — MPS, always backward (Rear).
%
% Applies the Method of Particular Solutions propagating backward from r2
% with odeMPCIrev, correcting the departure velocity at r2 until the backward
% trajectory lands at r1.  No direction selection — always backward. No Broyden.
%
% Usage:
%   [v1t,v2t,stats,history] = prtlambertMPSR( ...
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

if nargin < 9,   options = struct();  end
if ~isfield(options, 'Tolr'),    options.Tolr    = 1e-8;  end
if ~isfield(options, 'maxIter'), options.maxIter = 100;   end
if ~isfield(options, 'delta'),   options.delta   = 16;    end

Tolr    = options.Tolr;
maxIter = floor(options.maxIter / 4);   % 4 integrations per MPS iteration
do_print      = isfield(options, 'debug') && options.debug;
build_history = nargout >= 4;

if do_print, fprintf('starting MPSR\n'); end
if build_history
    history = struct('iterations', [], 'integration_number', [], 'normdr2', [], ...
                     'converged', false, 'final_error', NaN);
end

stats = struct('integration_number',0,'integration_preamble',0, ...
               'integration_phase1',0,'integration_phase2',0, ...
               'runtime_preamble',0,'runtime_phase1',0,'runtime_phase2',0, ...
               'direction_chosen',[],'hitearthvar',NaN,'hitrad',NaN);
t_ph1 = tic;

%% Preamble: unperturbed Lambert
[v1ref, v2ref] = lambertb(r1, r2, v1, dm, de, nrev, dtsec);
if ~all(isreal(v1ref)) || any(isnan(v1ref)) || ~all(isreal(v2ref)) || any(isnan(v2ref))
    v1t = NaN(1,3);  v2t = NaN(1,3);
    stats.runtime_phase1 = toc(t_ph1);
    return
end

x2ref = [r2 v2ref];
oe    = kep_elements(r2, v2ref, mu);
a     = oe(1);
if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
    Torb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
else
    Torb = sqrt(a^3/mu) * 2*pi;
end
options.Sec = Torb / options.delta;

[~, xbkw] = odeMPCIrev(forcemodel, [0 dtsec], x2ref, options);
integration_number = 1;
rdep_bkw = r1 - xbkw(1,1:3);

if build_history
    v1_check = xbkw(1,4:6);
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
    normdr2_init = norm(xfwd_c(end,1:3) - r2);
    history.iterations(end+1)         = 0;
    history.integration_number(end+1) = integration_number;
    history.normdr2(end+1)            = normdr2_init;
    if do_print, fprintf('init (MPSR): normdr1=%g, normdr2=%g\n', norm(rdep_bkw), normdr2_init); end
elseif do_print
    fprintf('init (MPSR): normdr1=%g\n', norm(rdep_bkw));
end

delta_size = 1e-7;
xpert      = zeros(3,6);
i = 1;

%% Backward MPS loop: correct v2ref until backward trajectory reaches r1
while norm(rdep_bkw) > Tolr && i < maxIter
    normv0    = norm(x2ref(4:6));
    deltardot = [zeros(3), eye(3)*delta_size*normv0];
    xprtrbd   = deltardot + x2ref;

    for j = 1:3
        [~, xtmp]  = odeMPCIrev(forcemodel, [0 dtsec], xprtrbd(j,:), options);
        xpert(j,:) = xtmp(1,:);
    end
    integration_number = integration_number + 3;
    if ~all(isfinite(xpert(:))), break, end

    DR = (xpert(:,1:3) - xbkw(1,1:3))';
    alpha = (DR \ rdep_bkw')';

    x2ref    = x2ref + alpha * deltardot;
    oe = kep_elements(x2ref(1:3), x2ref(4:6), mu);
    a  = oe(1);
    if a <= 0 || ~isreal(a) || isnan(a) || a > 2*r_max_allowed
        Torb = 2*pi*sqrt((2*r_max_allowed)^3/mu);
    else
        Torb = sqrt(a^3/mu) * 2*pi;
    end
    options.Sec = Torb / options.delta;

    [~, xbkw]  = odeMPCIrev(forcemodel, [0 dtsec], x2ref, options);
    integration_number = integration_number + 1;
    if ~all(isfinite(xbkw(1,:))), break, end
    rdep_bkw   = r1 - xbkw(1,1:3);

    if build_history
        v1_check = xbkw(1,4:6);
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
        history.iterations(end+1)         = i;
        history.integration_number(end+1) = integration_number;
        history.normdr2(end+1)            = normdr2_check;
        if do_print, fprintf('iteration %d (MPSR): normdr1=%g, normdr2=%g\n', i, norm(rdep_bkw), normdr2_check); end
    elseif do_print
        fprintf('iteration %d (MPSR): normdr1=%g\n', i, norm(rdep_bkw));
    end

    i = i + 1;
end

if norm(rdep_bkw) > Tolr || ~isfinite(norm(rdep_bkw))
    warning('prtlambertMPSR: did not converge');
    if build_history, history.converged = false; history.final_error = NaN; end
    stats.integration_number = integration_number;
    stats.integration_phase1 = integration_number;
    stats.runtime_phase1     = toc(t_ph1);
    v1t = NaN(1,3);  v2t = NaN(1,3);
    return
end

%% Extract v1 and compute v2 via final forward propagation
v1t = xbkw(1, 4:6);

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
