function [dr] = orbit_eq(t,r,mu)
%% function [dr] = orbit_eq(t,r,mu)
% definition of orbit eq...


% d^2r/dt = -mu*r/norm(r);

[n,m] = size(r);
dr = zeros (size(r));

if n == 6
    
    le =sqrt(sum(r(1:3,:).^2,1)).^3;

    dr(1,:) = r(4,:);
    dr(2,:) = r(5,:);
    dr(3,:) = r(6,:);
    dr(4,:) = -mu*r(1,:)./le;
    dr(5,:) = -mu*r(2,:)./le;
    dr(6,:) = -mu*r(3,:)./le;
elseif m == 6

    le =sqrt(sum(r(:,1:3).^2,2)).^3;

    dr(:,1) = r(:,4);
    dr(:,2) = r(:,5);
    dr(:,3) = r(:,6);
    dr(:,4) = -mu*r(:,1)./le;
    dr(:,5) = -mu*r(:,2)./le;
    dr(:,6) = -mu*r(:,3)./le;
end






