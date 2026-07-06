function [dr] = orbit_eq_J6_drag(t,r,mu1,CD,A,m,Re,J)
%% function [dr] = orbit_eq_J6_drag(t,r,mu,CD,A,m,Re,J)
% definition of orbit eq with J2-J6 and drag


% d^2r/dt = -mu*r/norm(r);
[n,M] = size(r);
dr = zeros (size(r));
if numel(J) ~= 6
    error("need 6 harmonics")
end

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


    ad = drag_accel(r',CD,A,m,Re)';

    aj2 = -(3*J(2)*mu1*Re^2)./(2*normr5).*[...
        r(1,:).*(1-5*r32tonormr2),...
        r(2,:).*(1-5*r32tonormr2),...
        r(3,:).*(3-5*r32tonormr2)]';
    aj3 = -(5*J(3)*mu1*Re^3)./(2*normr7).*[...
        r(1,:).*(3*r31-7*(r33)./(normr2)),...
        r(2,:).*(3*r31-7*(r33)./(normr2)),...
        6*r(3,:).^2-7*(r34)./(normr2)-3/5*normr2]';
    aj4 = (15*J(4)*mu1*Re^4)./(8*normr7).*[...
        r(1,:).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(2,:).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(3,:).*(5-(70/3)*r32tonormr2+21*r34tonormr4)]';
    aj5 = (3*J(5)*mu1*Re^5)./(8*normr9).*[...
        r(1,:).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r(2,:).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r32.*(105-315*r32tonormr2+231*r34tonormr4)-(15*J(5)*mu1*Re^5)./(8*normr7)]';
    aj6 = -(J(6)*mu1*Re^6)./(16*normr9).*[...
        r(1,:).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(2,:).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(3,:).*(245-2205*r32tonormr2+4851*r34tonormr4-3003*r36tonormr6)]';
    % ad = zeros(size(aj2));
    le = normr.^3;

    dr(1,:) = r(4,:);
    dr(2,:) = r(5,:);
    dr(3,:) = r(6,:);
    dr(4,:) = -mu1*r(1,:)./le + ad(1,:)+aj2(1,:)+aj3(1,:)+aj4(1,:)+aj5(1,:)+aj6(1,:);
    dr(5,:) = -mu1*r(2,:)./le + ad(2,:)+aj2(2,:)+aj3(2,:)+aj4(2,:)+aj5(2,:)+aj6(2,:);
    dr(6,:) = -mu1*r(3,:)./le + ad(3,:)+aj2(3,:)+aj3(3,:)+aj4(3,:)+aj5(3,:)+aj6(3,:);
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

    ad = drag_accel(r,CD,A,m,Re);
    aj2 = -(3*J(2)*mu1*Re^2)./(2*normr5).*[...
        r(:,1).*(1-5*r32tonormr2),...
        r(:,2).*(1-5*r32tonormr2),...
        r(:,3).*(3-5*r32tonormr2)];
    aj3 = -(5*J(3)*mu1*Re^3)./(2*normr7).*[...
        r(:,1).*(3*r31-7*(r33)./(normr2)),...
        r(:,2).*(3*r31-7*(r33)./(normr2)),...
        6*r(:,3).^2-7*(r34)./(normr2)-3/5*normr2];
    aj4 = (15*J(4)*mu1*Re^4)./(8*normr7).*[...
        r(:,1).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(:,2).*(1-14*r32tonormr2+21*r34tonormr4),...
        r(:,3).*(5-(70/3)*r32tonormr2+21*r34tonormr4)];
    aj5 = (3*J(5)*mu1*Re^5)./(8*normr9).*[...
        r(:,1).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r(:,2).*r31.*(35-210*r32tonormr2+231*r34tonormr4),...
        r32.*(105-315*r32tonormr2+231*r34tonormr4)-(15*J(5)*mu1*Re^5)./(8*normr7)];
    aj6 = -(J(6)*mu1*Re^6)./(16*normr9).*[...
        r(:,1).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(:,2).*(35-945*r32tonormr2+3465*r34tonormr4-3003*r36tonormr6),...
        r(:,3).*(245-2205*r32tonormr2+4851*r34tonormr4-3003*r36tonormr6)];
    % ad = zeros(size(aj2));
    le =normr.^3;

    dr(:,1) = r(:,4);
    dr(:,2) = r(:,5);
    dr(:,3) = r(:,6);
    dr(:,4) = -mu1*r(:,1)./le + ad(:,1)+aj2(:,1)+aj3(:,1)+aj4(:,1)+aj5(:,1)+aj6(:,1);
    dr(:,5) = -mu1*r(:,2)./le + ad(:,2)+aj2(:,2)+aj3(:,2)+aj4(:,2)+aj5(:,2)+aj6(:,2);
    dr(:,6) = -mu1*r(:,3)./le + ad(:,3)+aj2(:,3)+aj3(:,3)+aj4(:,3)+aj5(:,3)+aj6(:,3);
else
    error("need 6 numbers in state")
end






