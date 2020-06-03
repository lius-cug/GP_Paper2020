%% Estimating the direction of total magnetization using VG-TG method
%% Edited by Shuang Liu and Jamaledin Baniamerian
%% Reference: Dannemiller, N. and Y. Li, 2006, A new method for determination of magnetization direction:  Geophysics, 71, L69-L73, doi: 10.1190/1.2356116.
%%

clc
clear

I0=46;
D0=-4;

Inc=46;%%%%%%%%%%%%%%%%
Dec=-4;%%%%%%%%%%%%%%%%

[tpX,tpY,TR,zmin,zmax,xin,yin] =read_ASCII_grid('Za.grd');
deltax = xin(2)-xin(1);deltay = yin(2)-yin(1);
[nx,ny] = size(TR);
nmax=max([nx ny]);
npts=2^nextpow2(nmax);

subplot(2,2,1);
contourf(TR);figure(gcf);
axis equal;
axis image;
colorbar();
title ('Magnetic anomaly');

% out = RTP(TR,Inc,Dec,deltax,deltay);%A///////////////////////////////////
RToP=RTP_grid(TR,I0,D0,Inc,Dec,50,0,deltax,'sp0');
subplot(2,2,2);
contourf(RToP);figure(gcf);
axis equal;
axis image;
colorbar();
title ('RTP');


maxr=-10^38;
for i=1:91
    for j=1:181
        incl=i*1-1+0.1
        decl=j*1-91+0.1

        tpRTP=ZcompRTP(TR,I0,D0,incl,decl,50,0,xin,yin,'sp0');      
        [f1z,f2z,Gx1,Gx2,Gy1,Gy2]=VDerGridd_ST(tpRTP,xin,yin,50,'sp0');
        AS=sqrt(Gx1.*Gx1+Gy1.*Gy1+f1z.*f1z);
        
        
        tpR=corrcoef(f1z, AS);
        if  tpR(1,2)>maxr
            inc_opt=incl;
            dec_opt=decl;
            maxr=tpR(1,2);
        end
        R(i,j)=power(tpR(1,2),11);      
    end
       
end
R(1,1)=0;
inc_opt
dec_opt


subplot(2,2,3);
contourf(R);figure(gcf);
set(gca,'xtick',[0:10:90]);
set(gca,'ytick',[0:10:90]);
xlabel('dec');
ylabel('inc');
axis equal;
axis image;
colorbar();
title ('orrelation');
grid on;

str=num2str([Inc,Dec,inc_opt,dec_opt],'Results: inc_true=%.1f, dec_true=%.1f, inc_cal=%.1f, dec_cal=%.1f\n');
set(gcf,'Name',str);

dectp=-90:90;
inctp=0:90;
write_ASCII_grid('Correlation.grd',R,dectp,inctp) ;




