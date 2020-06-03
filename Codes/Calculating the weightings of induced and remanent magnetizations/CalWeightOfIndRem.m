%% Calculating the weightings of induced and remanent magnetizations 
%% Edited by Shuang Liu and Jamaledin Baniamerian
%%
clc
clear

Ii=46;
Di=-4;

Ir=50.1337;
Dr=-53.9564;

It=51;
Dt=-30;

Q=1.17;

dpi=pi/180;

i=[cos(Ii*dpi)*cos(Di*dpi),cos(Ii*dpi)*sin(Di*dpi),sin(Ii*dpi)];

r=[cos(Ir*dpi)*cos(Dr*dpi),cos(Ir*dpi)*sin(Dr*dpi),sin(Ir*dpi)];

t=[cos(It*dpi)*cos(Dt*dpi),cos(It*dpi)*sin(Dt*dpi),sin(It*dpi)];

cosa=i*t';

cosb=r*t';

w_i=cosa/(cosa+Q*cosb)

w_r=Q*cosb/(cosa+Q*cosb)