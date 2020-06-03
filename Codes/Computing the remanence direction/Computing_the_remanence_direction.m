%% Computing the remanence direction
%% Edited by Shuang Liu and Jamaledin Baniamerian
 
 
clc
clear


dpi=pi/180;
H0=0.5;

%%%%%%%%%%Input%%%%%%
I0=46;
D0=-4;

It=51;
Dt=-30;

Q=1.17;
%%%%%%%%%Input%%%%%%

Ir=0;Dr=0;
UnitVecInd=[cos(I0*dpi)*cos(D0*dpi);cos(I0*dpi)*sin(D0*dpi);sin(I0*dpi)];
UnitVecTot=[cos(It*dpi)*cos(Dt*dpi);cos(It*dpi)*sin(Dt*dpi);sin(It*dpi)];
UnitVecR=[cos(Ir*dpi)*cos(Dr*dpi);cos(Ir*dpi)*sin(Dr*dpi);sin(Ir*dpi)];
ht=UnitVecInd'*UnitVecTot;

if Q>=1 
    tp=ht+sqrt(ht^2-1+Q^2);
    UnitVecRes=(tp*UnitVecTot-UnitVecInd )/Q;
else 
    tp=ht+sqrt(ht^2-1+Q^2);   %%%%%%% two solutions tp=ht + - sqrt(ht^2-1+Q^2);
    UnitVecRes=(tp*UnitVecTot-UnitVecInd )/Q;
end

%UnitVecRes
length= UnitVecRes(1,1)*UnitVecRes(1,1)+UnitVecRes(2,1)*UnitVecRes(2,1)+UnitVecRes(3,1)*UnitVecRes(3,1);
Ir=asin(UnitVecRes(3,1))/dpi;
Dr=-acos(UnitVecRes(1,1)/sqrt(UnitVecRes(1,1)*UnitVecRes(1,1)+UnitVecRes(2,1)*UnitVecRes(2,1)) )/dpi;
RemID=[Ir; Dr]

