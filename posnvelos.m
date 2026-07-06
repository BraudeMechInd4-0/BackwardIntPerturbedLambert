function [r,v] = posnvelos(a,e,i,O,w,f,mu)
%%[r,v] = posnvelos(a,e,i,O,w,f,mu)
%extracts the position and velocity of a set of known keplerian parameters
%in the ECI system (all angles are in radians)

%tet = wrapTo2Pi(f+w);
cosf = cos(f);
sinf = sin(f);
sinw = sin(w);
cosw = cos(w);
cosO = cos(O);
sinO = sin(O);
%costet = cos(tet);
%sintet = sin(tet);
cosi = cos(i);
sini = sin(i);
p = a*(1-e^2);
rPQW = [p*cosf/(1+e*cosf); p*sinf/(1+e*cosf);0];
vPQW = [-sqrt(mu/p)*sinf; sqrt(mu/p)*(e+cosf);0];

R = [cosO*cosw-sinO*sinw*cosi, -cosO*sinw-sinO*cosw*cosi, sinO*sini;...
     sinO*cosw+cosO*sinw*cosi, -sinO*sinw+cosO*cosw*cosi, -cosO*sini;...
     sinw*sini,                 cosw*sini,                cosi];

r = (R*rPQW)';
v = (R*vPQW)';