%% RTP transform
%%%% Inputs
% 
% f= field
% If,Df: Inc and Dec of the field
% Im,Dm: Inc and Dec of Magnetization
% cell is grid data cell size
% expand= grid extension in percent. different type of extension could be used
% A  azimuth of grid respect to geographic north. usually A=0
%      mode=
%     'zpd' zero extension.
%     'sp0' smooth extension of order 0.
%     'spd' (or 'sp1') smooth extension of order 1.
%     'sym' (or 'symh') symmetric extension (half-point).
%     'symw' symmetric extension (whole-point).
%     'asym' (or 'asymh') antisymmetric extension (half-point).
%     'asymw' antisymmetric extension (whole-point).
%     'ppd' periodized extension (1).
%     'per' periodized extension (2).
%     'inp1','inp2','inp3','inp4' (inpaint_Nans)
%     expand is given in per cent
%% Ref: Blakely 1995.
%%
function [frtp]=RTPgrid(f,If,Df,Im,Dm,expand,A,Cell,mode)

Im=Im*pi/180; Dm=Dm*pi/180;  If=If*pi/180; Df=Df*pi/180; A=A*pi/180;
kk=sqrt(-1);

%%
[fexp,T1]=grid_extension(f,expand,mode);
 %%
[mnew,nnew]=size(fexp);

%,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
F=fftshift(fft2(fexp));
[~,~,KX,KY,K]=kvalue2(mnew,nnew,Cell,Cell);

[ind1,ind2]=find(K==0); 
% if isnan(ind1)
%     hhhhh=1;
% else
%    
K(ind1,ind2)=mean(K(ind1,ind2-1)+K(ind1,ind2+1)+K(ind1-1,ind2)+K(ind1+1,ind2));
KX(ind1,ind2)=0.5*mean(abs(KX(ind1,ind2-1))+abs(KX(ind1,ind2+1))+abs(KX(ind1-1,ind2))+abs(KX(ind1+1,ind2)));
KY(ind1,ind2)=0.5*mean(abs(KY(ind1,ind2-1))+abs(KY(ind1,ind2+1))+abs(KY(ind1-1,ind2))+abs(KY(ind1+1,ind2)));
% end
%,,,,,,,,,,,,,,,,,,,,,,, RTP Operator

a1=sin(If)*sin(Im)-cos(If)*cos(Df-A)*cos(Im)*cos(Dm-A);
a2=sin(If)*sin(Im)-cos(If)*sin(Df-A)*cos(Im)*sin(Dm-A);
a3=-cos(If)*cos(Df-A)*cos(Im)*sin(Dm-A)-cos(If)*sin(Df-A)*cos(Im)*cos(Dm-A);
b1=cos(Im)*cos(Dm-A)*sin(If)+sin(Im)*cos(If)*cos(Df-A);
b2=cos(Im)*sin(Dm-A)*sin(If)+sin(Im)*cos(If)*sin(Df-A);

rtpop=(K.^2)./(a1*(KY.^2)+a2*(KX.^2)+a3*(KX.*KY)+kk*K.*(b1*KY+b2*KX));

Fourier_rtp=F.*rtpop;
frtp=real(ifft2(fftshift(Fourier_rtp)));

frtp=frtp(T1(1):T1(2),T1(3):T1(4));
end

function [fexp,t]=grid_extension(f,expand,mode)
method=mode;

[m,n]=size(f);


num=max(m,n);
newsize=floor(num+(num*expand/100));

if mod(newsize,2)==0
    twopower=newsize;
elseif mod(newsize,2)==1
    twopower=newsize+1;
end
[fexp,t]=bordowav(f,method,twopower);
end



function [ky,kx,KX,KY,K]=kvalue2(m,n,dx,dy)
dkx=2*pi/(dx*n);
dky=2*pi/(dy*m);
kx=dkx*(-n/2:n/2-1);
ky=dky*(-m/2:m/2-1);
[KX,KY]=meshgrid(kx,ky);
K=sqrt(KX.^2+KY.^2);
end


function [grid4,t]=bordowav(a,mode,twopower)
% function [grid4,t]=bordowav(a,mode,twopower)
% output:
%   grid4 = matrice bordata
%   t = indici riga e colonna per debordare
%   a=matrice da bordare
%   twopower=potenza di 2 (o altra dimensione)
%mode='zpd' zero extension.
%     'sp0' smooth extension of order 0.
%     'spd' (or 'sp1') smooth extension of order 1.
%     'sym' (or 'symh') symmetric extension (half-point).
%     'symw' symmetric extension (whole-point).
%     'asym' (or 'asymh') antisymmetric extension (half-point).
%     'asymw' antisymmetric extension (whole-point).
%     'ppd' periodized extension (1).
%     'per' periodized extension (2).
%     'inp1','inp2','inp3','inp4' (inpaint_Nans)
if mod(twopower-size(a,1),2)~= 0;
    s1=(twopower-size(a,1)+1)/2;
    t1=s1+1;t2=t1+size(a,1)-1;
else s1=(twopower-size(a,1))/2;
    t1=s1+1;t2=t1+size(a,1)-1;
end
if mod(twopower-size(a,2),2)~= 0;
    s2=(twopower-size(a,2)+1)/2;
    t3=s2+1;t4=t3+size(a,2)-1;
else s2=(twopower-size(a,2))/2;
    t3=s2+1;t4=t3+size(a,2)-1;
end
if strcmp(mode, 'inp1')
  c=NaN(twopower,twopower);
  c(t1:t2,t3:t4)=a(:,:);
  grid4=inpaint_nans(c,1);  
elseif strcmp(mode, 'inp2')
  c=NaN(twopower,twopower);
  c(t1:t2,t3:t4)=a(:,:);
  grid4=inpaint_nans(c,2);  
elseif strcmp(mode, 'inp3')
  c=NaN(twopower,twopower);
  c(t1:t2,t3:t4)=a(:,:);
  grid4=inpaint_nans(c,3);  
elseif strcmp(mode, 'inp4')
  c=NaN(twopower,twopower);
  c(t1:t2,t3:t4)=a(:,:);
  grid4=inpaint_nans(c,4);  
else   
  grid4=wextend(2,mode,a,[s1,s2]);
  grid4=grid4(1:twopower,1:twopower);
end
t=[t1,t2,t3,t4];
end
