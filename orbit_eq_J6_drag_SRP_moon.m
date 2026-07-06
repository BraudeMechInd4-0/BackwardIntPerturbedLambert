function [dr] = orbit_eq_J6_drag_SRP_moon(t,r,muE,CD,Adrag,m,Re,J,jdepoch,rhoSRP,CR,ASRP,muM)%,pp_sun,pp_moon)
%% function [dr] = orbit_eq_J6_drag(t,r,mu,CD,A,m,Re,J,jdepoch,rhoSRP,CR,ASRP,muM) ?? ,pp_sun,pp_moon)
% definition of orbit eq with J2-J6 and drag


% d^2r/dt = -mu*r/norm(r);
[n,M] = size(r);
dr = zeros (size(r));
if numel(J) ~= 6
    error("need 6 harmonics")
end
AU = 149597870.7;

jdNow  = jdepoch + t(:)' / 86400;   % 1×N row vector
r_moon = moon_v(jdNow) * Re;        % N×3
r_sol  = sun_v(jdNow)  * AU;        % N×3

if n == 6

    normr = sqrt(sum(r(1:3,:).^2,1));
    normr2 = normr.^2;
    %normr3 = normr.^3;
    normr4 = normr.^4;
    normr5 = normr.^5;
    normr6 = normr.^6;
    normr7 = normr.^7;
    %normr8 = normr.^8;
    normr9 = normr.^9;

    r31 = r(3,:);
    r32 = r(3,:).^2;
    r33 = r(3,:).^3;
    r34 = r(3,:).^4;
    %r35 = r(3,:).^5;
    r36 = r(3,:).^6;
    %r37 = r(3,:).^7;
    %r38 = r(3,:).^8;
    %r39 = r(3,:).^9;
    
    r32tonormr2 = (r32)./(normr2);
    r34tonormr4 = (r34)./(normr4);
    r36tonormr6 = (r36)./(normr6);


    ad = drag_accel(r',CD,Adrag,m,Re)';

    aj2 = -(3*J(2)*muE*Re^2)./(2*normr5).*[...
        r(1,:).*(1-5*r32tonormr2),...
        r(2,:).*(1-5*r32tonormr2),...
        r(3,:).*(3-5*r32tonormr2)]';
    aj3 = -(5*J(3)*muE*Re^3)./(2*normr7).*[...
        r(1,:).*(3*r31-7*(r33)./(normr2)),...
        r(2,:).*(3*r31-7*(r33)./(normr2)),...
        6*r(3,:).^2-7*(r34)./(normr2)-3/5*normr2]';
    aj4 = (15*J(4)*muE*Re^4)./(8*normr7).*[...
        r(1,:).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(2,:).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(3,:).*(5-(70/3)*r32tonormr2+21*r34tonormr4)]';
    aj5 = (3*J(5)*muE*Re^5)./(8*normr9).*[...
        r(1,:).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r(2,:).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r32.*(105-315*r32tonormr2+231*r34tonormr4)-(15*J(5)*muE*Re^5)./(8*normr7)]';
    aj6 = -(J(6)*muE*Re^6)./(16*normr9).*[...
        r(1,:).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(2,:).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(3,:).*(245-2205*r32tonormr2+4851*r34tonormr4-3003*r36tonormr6)]';
    
    r_moon_T   = r_moon';                                          % 3×N
    r_sol_T    = r_sol';                                           % 3×N
    r_sol_norm = sqrt(sum(r_sol.^2, 2))';                          % 1×N
    aSRP       = -(rhoSRP*CR*ASRP/m) .* r_sol_T ./ r_sol_norm;    % 3×N

    rsatM      = r_moon_T - r(1:3,:);                              % 3×N
    normrsatM  = sqrt(sum(rsatM.^2,  1));                          % 1×N
    normr_moon = sqrt(sum(r_moon.^2, 2))';                         % 1×N
    amoon = muM*(rsatM./normrsatM.^3 - r_moon_T./normr_moon.^3);  % 3×N

    le = normr.^3;

    
    dr(1:3,:) = r(4:6,:);
    dr(4:6,:) = -muE*r(1:3,:)./le+ad+aj2+aj3+aj4+aj5+aj6+aSRP+amoon;
    
elseif M == 6
    normr = sqrt(sum(r(:,1:3).^2,2));
    normr2 = normr.^2;
    %normr3 = normr.^3;
    normr4 = normr.^4;
    normr5 = normr.^5;
    normr6 = normr.^6;
    normr7 = normr.^7;
    %normr8 = normr.^8;
    normr9 = normr.^9;

    r31 = r(:,3);
    r32 = r(:,3).^2;
    r33 = r(:,3).^3;
    r34 = r(:,3).^4;
    %r35 = r(:,3).^5;
    r36 = r(:,3).^6;
    %r37 = r(:,3).^7;
    %r38 = r(:,3).^8;
    %r39 = r(:,3).^9;
    
    r32tonormr2 = (r32)./(normr2);
    r34tonormr4 = (r34)./(normr4);
    r36tonormr6 = (r36)./(normr6);

    ad = drag_accel(r,CD,Adrag,m,Re);
    aj2 = -(3*J(2)*muE*Re^2)./(2*normr5).*[...
        r(:,1).*(1-5*r32tonormr2),...
        r(:,2).*(1-5*r32tonormr2),...
        r(:,3).*(3-5*r32tonormr2)];
    aj3 = -(5*J(3)*muE*Re^3)./(2*normr7).*[...
        r(:,1).*(3*r31-7*(r33)./(normr2)),...
        r(:,2).*(3*r31-7*(r33)./(normr2)),...
        6*r(:,3).^2-7*(r34)./(normr2)-3/5*normr2];
    aj4 = (15*J(4)*muE*Re^4)./(8*normr7).*[...
        r(:,1).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(:,2).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(:,3).*(5-(70/3)*r32tonormr2+21*r34tonormr4)];
    aj5 = (3*J(5)*muE*Re^5)./(8*normr9).*[...
        r(:,1).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r(:,2).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r32.*(105-315*r32tonormr2+231*r34tonormr4)-(15*J(5)*muE*Re^5)./(8*normr7)];
    aj6 = -(J(6)*muE*Re^6)./(16*normr9).*[...
        r(:,1).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(:,2).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(:,3).*(245-2205*r32tonormr2+4851*r34tonormr4-3003*r36tonormr6)];
    
    r_sol_norm = sqrt(sum(r_sol.^2, 2));                            % N×1
    aSRP       = -(rhoSRP*CR*ASRP/m) .* r_sol ./ r_sol_norm;      % N×3

    rsatM      = r_moon - r(:,1:3);                                % N×3
    normrsatM  = sqrt(sum(rsatM.^2,  2));                          % N×1
    normr_moon = sqrt(sum(r_moon.^2, 2));                          % N×1
    amoon = muM*(rsatM./normrsatM.^3 - r_moon./normr_moon.^3);    % N×3

    le =normr.^3;

    dr(:,1:3) = r(:,4:6);
    dr(:,4:6) = -muE*r(:,1:3)./le+ad(:,1:3)+aj2+aj3+aj4+aj5+aj6+aSRP+amoon;
else
    error("need 6 numbers in state")
end






