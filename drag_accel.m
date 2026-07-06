function [ad] = drag_accel(x,CD,A,m,Re)
%A function that computes the acceleration contribution due to drag
%Usage:
%   function [ad] = drag_accel(x,CD,A,m)
%       x - the curent state [rx ry rz vx vy vz]
%       CD - drag coeeficient
%       A - cross section of satellite
%       m - satellite mass
%       Re - earth's radius

v = x(:,4:6);
r = x(:,1:3);
rnorm = sqrt(sum(r.^2,2));
vnorm = sqrt(sum(v.^2,2));
ad = zeros(size(x,1),3);

h = rnorm - Re;

if any(h < 100)
    warning("Satellite decayed!")
end
h(h<=0) = 0;
[h0,r0,H] = expDens(h);
rho = r0.*exp(-(h-h0)./H)*1e9;%kg/km^3
ad = -0.5*CD*A/m*rho.*vnorm.*v; 


function [h0,r0,H] = expDens(h) %make sure you get the sizes of rho that are good for your equation (i.e. column or line)
if length(h) > 1
    h(h>1000) = 1000;
    htable = [0 25 30:10:150 180 200:50:500 600:100:1000 inf]';
    rtable = [1.225, 3.899e-2, 1.774e-2, 3.972e-2, 1.057e-3,3.206e-4, 8.770e-5, 1.905e-5,...
        3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9, 5.464e-10,...
        2.789e-10, 7.248e-11, 2.418e-11, 9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, ...
        3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15]'; %kg/m^3
    Htable = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799, 5.382, 5.877, 7.263,...
        9.473, 12.636, 16.149, 22.523, 29.740, 37.105, 45.546, 53.628, 53.298, 58.515,...
        60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00]';
    [~,I]=mink(h'>=htable,1,1); %a little ""hack"": create a matrix for which each column corresponds to a position in h.  
    I = max(I-1,1); % we look for the first 0, and then go back one
    h0 = htable(I);
    r0 = rtable(I);
    H = Htable(I);
else
    if h>1000
        h0 = 1000;
        r0 = 3.019e-15;
        H = 268;
    elseif h>900 && h<=1000
        h0 = 900;
        r0 = 5.245e-15;
        H = 181.05;
    elseif h>800 && h<=900
        h0 = 800;
        r0 = 1.170e-14;
        H = 124.64;
    elseif h>700 && h<=800
        h0 = 700;
        r0 = 3.614e-14;
        H  = 88.667;
    elseif h>600 && h<=700
        h0 = 600;
        r0 = 1.454e-13;
        H = 71.835;
    elseif h>500 && h<=600
        h0 = 500;
        r0 = 6.967e-13;
        H = 63.822;
    elseif h>450 && h<=500
        h0 = 450;
        r0 = 1.585e-12;
        H = 60.828;
    elseif h>400 && h<=450
        h0 = 400;
        r0 = 3.725e-12;
        H = 58.515;
    elseif h>350 && h<=400
        h0 = 350;
        r0 = 9.518e-12;
        H = 53.298;
    elseif h>300 && h<=350
        h0 = 300;
        r0 = 2.418e-11;
        H = 53.628;
    elseif h>250 && h<=300
        h0 = 250;
        r0 = 7.248e-11;
        H = 45.546;
    elseif h>200 && h<=250
        h0 = 200;
        r0 = 2.789e-10;
        H = 37.105;
    elseif h>180 && h<=200
        h0 = 180;
        r0 = 5.464e-10;
        H = 29.740;
    elseif h>150 && h<=180
        h0 = 150;
        r0 = 2.070e-9;
        H = 22.523;
    elseif h>140 && h<=150
        h0 = 140;
        r0 = 3.845e-9;
        H = 16.149;
    elseif h>130 && h<=140
        h0 = 130;
        r0 = 8.484E-9;
        H = 12.636;
    elseif h>120 && h<=130
        h0 = 120;
        r0 = 2.438e-8;
        H = 9.473;
    elseif h>110 && h<=120
        h0 = 110;
        r0 = 9.661e-8;
        H = 7.263;
    elseif h>100 && h<=110
        h0 = 100;
        r0 = 5.297e-7;
        H = 5.877;
    elseif h>90 && h<=100
        h0 = 90;
        r0 = 3.396e-6;
        H = 5.382;
    elseif h>80 && h<=90
        h0 = 80;
        r0 = 1.905e-4;
        H=5.779;
    elseif h>70 && h<=80
        h0 = 70;
        r0 = 8.770e-5;
        H = 6.549;
    elseif h>60 && h<=70
        h0 = 60;
        r0 = 3.206e-4;
        H=7.714;
    elseif h>50 && h<=60
        h0 = 50;
        r0 = 1.057e-3;
        H = 8.382;
    elseif h>40 && h<=50
        h0 = 40;
        r0 = 3.972e-3;
        H = 7.554;
     elseif h>30 && h<=40
        h0 = 30;
        r0 = 1.774e-2;
        H=6.682;
    elseif h>25 && h < 30
        h0 = 25;
        r0 = 3.899e-2;
        H=6.349;
    elseif h >= 0 && h <=25
        h0 = 0;
        r0 = 1.225;
        H=7.249;
    else
        error('bad input');
    end
end