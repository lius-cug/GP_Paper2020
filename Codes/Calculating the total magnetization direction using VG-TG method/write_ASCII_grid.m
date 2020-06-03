% **************************************************************
% Potensoft - Gravity and Magnetic Data Mapping, Processing, 
% Modeling and Filtering Software
% Tested using Matlab 7.6 Release 2008b under Windows XP
% **************************************************************
%           AUTHORS
% M. Ozgu Arisoy & Unal Dikmen - 2009
% arisoy@eng.ankara.edu.tr
% dikmen@eng.ankara.edu.tr
%
% Ankara University, Faculty of Engineering, 
% Department of Geophysical Engineering, 06100, Ankara-TURKEY
% http://eng.ankara.edu.tr/~arisoy/potensoft.htm
% **************************************************************
% Function : write_ASCII_grid
% Creates a Surfer ASCII grid file
% **************************************************************

function write_ASCII_grid(grid_file,data,xin,yin) 

[nrr,ncc] = size(data);
if length(xin) == nrr
    xinb = xin;
    yinb = yin;     
    yin = xinb;
    xin = yinb;
end

fid = fopen(grid_file,'w');
fprintf(fid,'DSAA\n');
siz=size(data');
nx = siz(1);
ny = siz(2);
fprintf(fid,'%i ',siz);
fprintf(fid,'\n%i %i',min(xin),max(xin));
fprintf(fid,'\n%i %i',min(yin),max(yin));
fprintf(fid,'\n%i %i',min(data(:)),max(data(:)));
fprintf(fid,'\n');
for jj=1:ny                                 
    for ii=1:nx
        fprintf(fid,'%g %c',data(jj,ii),' ');
    end
    fprintf(fid,'\n');
end
fclose(fid);