% 
% --> function [A,B,Fe]=Nemoh(w, dir, depth)
%
% Purpose: Matlab wrapper for calculation of hydrodynamic coefficients using Nemoh
% 
% Inputs :
% - w     : Vector length(w) of wave frequencies (rad/s)
% - dir   : Wave direction (degrees)
% - depth : water depth (m), 0 for deep water.
%
% Outputs :
% - A  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of added mass coefficients
% - B  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of radiation damping coefficients
% - Fe : Matrix (6xnBodies)xlength(w) of excitation forces (complex
% values)
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.
%
function [A,B,Fe]=Nemoh(w, dir, depth)
% Preparation du calcul
rep='.';

fid=fopen([rep,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
fclose(fid);
fid=fopen([rep,filesep,'Nemoh.cal'],'r');
n=1;
clear textline;
textline={};
while (~feof(fid))
    textline(n)={fgetl(fid)};
    if (n == 4) 
        textline(n)={sprintf('%f                 ! DEPTH			! M		! Water depth',depth)};
    end
    if ((mod(n,18) == 9) && ((n-9)/18 <= nBodies))
        temp=cell2mat(textline(n));
        temp2=[];
        ntemp=length(temp);
        k=1;
        for i=1:ntemp
            if (temp(i) == '\')
                temp2=[temp2,temp(k:i),'\'];
                k=i+1;
            end;            
        end
        temp2=[temp2,temp(k:ntemp)];
        textline(n)={temp2};
        cell2mat(textline(n));
    end
    if (n == 9+18*nBodies) 
        textline(n)={sprintf('%g %f %f       ! Number of wave frequencies, Min, and Max (rad/s)',length(w),w(1),w(length(w)))};
    end
     if (n == 10+18*nBodies) 
        textline(n)={sprintf('%g %f %f		! Number of wave directions, Min and Max (degrees)',1,dir,dir)};
    end
    n=n+1;
end
fclose(fid);
fid = fopen([rep,filesep,'Nemoh.cal'], 'w'); 
for i=1:n-1
    fprintf(fid, [cell2mat(textline(i)),'\n']);
end
fclose(fid);
fid=fopen([rep,filesep,'input.txt'],'wt');
fprintf(fid,' \n 0 \n');
status=fclose(fid);
% Calcul des coefficients hydrodynamiques
system(['capytaine ' rep filesep 'Nemoh.cal']);
%% Lecture des resultats CA CM Fe
clear Periode A B Famp Fphi Fe;
fid=fopen([rep,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
for i=1:2+18*nBodies
    ligne=fgetl(fid);
end
nw=fscanf(fid,'%g',1);
fclose(fid);
fid=fopen([rep,filesep,'results',filesep,'ExcitationForce.tec'],'r');
ligne=fgetl(fid);
for c=1:6*nBodies
    ligne=fgetl(fid);
end;
ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+12*nBodies);
    w(i)=ligne(1);
    for j=1:6*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end;
end;
status=fclose(fid);
fid=fopen([rep,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
ligne=fgetl(fid);
for i=1:6*nBodies
    ligne=fgetl(fid);
end;
for i=1:nBodies*6
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+12*nBodies);
        for j=1:6*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end;
        ligne=fgetl(fid);
    end;
end;
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
i=sqrt(-1);
Fe=Famp.*(cos(Fphi)+i*sin(Fphi));
end
