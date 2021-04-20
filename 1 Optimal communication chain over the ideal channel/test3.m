clear all; clc; close all;

f = 3.6*10^9;
c = 3*10^8;
lambda = c/f;
er = 4.5

d = sqrt(40^2 + 10^2)+sqrt(30^2 + 5^2)
o1 = atan(5/30)
o2 = atan(35/70)
% o1 = pi/2;
o2 = pi/2;
d_0 = sqrt(70^2 + 15^2)
Dr = d - d_0


t = d/c
g1 = (cos(o1) - sqrt(er*(1-(1/er)*sin(o1)^2)))/(cos(o1) + sqrt(er*(1-(1/er)*sin(o1)^2)))
g2 = (cos(o2) - sqrt(er*(1-(1/er)*sin(o2)^2)))/(cos(o2) + sqrt(er*(1-(1/er)*sin(o2)^2)))

nu = sqrt(4*Dr/lambda)

F2dB = -6.9 -20*log10(sqrt((nu-0.1)^2 +1) +nu -0.1) 

F2 = 10^(F2dB/20)
Fa = mod((-1/4 -(nu^2)/2),2)

mod_e = (sqrt(F2)*g1*g2*sqrt(60)/d)
arg_e = mod((2*d)/(lambda),2)

arg_E = arg_e + Fa
