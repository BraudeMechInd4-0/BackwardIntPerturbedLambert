function [a,e,i,O,w,f] = kep_elements(r,v,mu)
%%function [a,e,i,O,w,f] = kep_elements(r,v,mu)
%extracts the keplerian orbital elements of a known position and velocity
%in the ECI system
%THIS FUNCTION WORKS ON TWO SINGLE LINE 1x3 VECTORS!!!

if nargin < 3
    mu = 398600.5;
end

%% known vectors:
%ECI reference frame
iv = [1 0 0];
%jv = [0 1 0];
kv = [0 0 1];

norm_v = norm(v);
norm_r = norm(r);
%E = norm_v^2/2-mu/norm_r;
%h
h = cross(r,v);
norm_h = norm(h);
%n
n = cross(kv,h);
norm_n = norm(n);

%% a
a = (2/norm_r -norm_v^2/mu)^-1;

%% e
%vector
ev = cross(v,h)/mu-r/norm_r;
norm_ev = norm(ev);
e = norm(ev);

%% i
i = acos(dot(kv,h)/norm_h);

%% O
O = acos(dot(iv,n)/norm_n);
if n(2) < 0
    O = 2*pi - O;
end

%% w
w = acos(dot(n,ev)/(norm_n*norm_ev));
if ev(3) < 0
    w = 2*pi - w;
end

%% f
f = acos(dot(ev,r)/(norm_ev*norm_r));
if dot(r,v) < 0
    f = 2*pi - f;
end


if nargout == 1
    a = [a,e,i,O,w,f];
end


