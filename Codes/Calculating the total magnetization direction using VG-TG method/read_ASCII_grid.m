% ***********************************************************
% Potensoft - Gravity and Magnetic Data Mapping, Processing, 
% Modeling and Filtering Software
% Tested using Matlab 7.6 Release 2008b under Windows XP
% ***********************************************************
%           AUTHORS
% M. Ozgu Arisoy & Unal Dikmen - 2009
% arisoy@eng.ankara.edu.tr
% dikmen@eng.ankara.edu.tr
%
% Ankara University, Faculty of Engineering, 
% Department of Geophysical Engineering, 06100, Ankara-TURKEY
% http://eng.ankara.edu.tr/~arisoy/potensoft.htm
% ***********************************************************
% Function : read_ASCII_grid
% Reads a ASCII Surfer grid file
% ***********************************************************

function [X,Y,TR,zmin,zmax,x,y] = read_ASCII_grid(grid_file)

fid=fopen(grid_file,'r');
if fid == -1
     messtext = [' Error: Grid file must be in the project path.'];
    errordlg(messtext,'Error!')
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    return
end
s=fgetl(fid);
char_s = ischar(s);
if char_s == 1
    [nr_s,nc_s] = size(s);
    if nr_s == 1 && nc_s == 4
        if ~strcmp(s,'DSAA')  
        X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
        error_mess;
        return
        end
    else
        X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
        error_mess;
        return
    end
else
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    error_mess;
    return
end
gridsize=fscanf(fid,'%f',2);
[nr_grid,nc_grid] = size(gridsize);
if nr_grid ~= 2 || nc_grid ~= 1
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    error_mess;
    return
end
mx=fscanf(fid,'%f',2);mx(3)=mx(2)-mx(1);
[nr_mx,nc_mx] = size(mx);
if nr_mx ~= 3 || nc_mx ~= 1
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    error_mess;
    return
end
my=fscanf(fid,'%f',2);my(3)=my(2)-my(1);
[nr_my,nc_my] = size(my);
if nr_my ~= 3 || nc_my ~= 1
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    error_mess;
    return
end
hmmm=fscanf(fid,'%f',2);
[nr_hmmm,nc_hmmm] = size(hmmm);
if nr_hmmm ~= 2 || nc_hmmm ~= 1
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    error_mess;
    return
end
zmin=hmmm(1);
zmax=hmmm(2);
if zmax < zmin
    X=[];    Y=[];     TR=[];    zmin=0;    zmax=0;     x=[];     y=[];
    error_mess;
    return
end
x=mx(1)+mx(3)*(0:gridsize(1)-1)/(gridsize(1)-1);
y=my(1)+my(3)*(0:gridsize(2)-1)/(gridsize(2)-1);
A=fscanf(fid,'%f',inf);
[nr_A,nc_A] = size(A);
size_A = gridsize(1)*gridsize(2);
if nr_A == size_A
    A=reshape(A,gridsize')';
    fclose(fid);
    TR = A;
    [X,Y] = meshgrid(x,y);
else
    X=[];    Y=[];     TR=[]; zmin=0; zmax=0;
    error_mess;
    return
end

function error_mess

messtext = ' Error: It is not a valid grid file                ';                 
    errordlg(messtext,'Error!')

