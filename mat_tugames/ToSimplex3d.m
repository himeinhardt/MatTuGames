function [sim,sX,oX]=ToSimplex3d(X)
% TOSIMPLEX3D projects data from 4d to 3d as well as from 3d to 2d.
%
% Source: http://stackoverflow.com/questions/3506982/projecting-points-from-4d-space-into-3d-space-in-mathematica
%
% Usage: [sim,sX,oX]=ToSimplex3d(X)
%
% Define variables:
% output:
% sim      -- Data projected on the 3d- or 2d-simplex.
% sX       -- Raw data projected on the 3d- or 2d-simplex in 
%             the 4d- or 3d-space. Holding the last column constant.
% oX       -- Simplex data reprojected from 3d to 4d or from 3d to 2d. 
%             Must reproduce the data set X.
%
% input:   
% X        -- A data set (vector/matrix) of 4d or 3d.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/12/2015        0.7             hme
%


[sr sc]=size(X);
if sc==4
   A=[-(1/2), -(1/(2*sqrt(3))), -(1/(2*sqrt(6))),1/sqrt(4); ...
     1/2, -(1/(2*sqrt(3))), -(1/(2*sqrt(6))),1/sqrt(4); ...
     0, -(1/(2*sqrt(3))) + sqrt(3)/2, -(1/(2*sqrt(6))),1/sqrt(4); ...
     0, 0, sqrt(2/3) - 1/(2*sqrt(6)), 1/sqrt(4)];
   sX=X*A;
   B=A\eye(4);
   oX=sX*B;
   sim=sX;
   sim(:,4)=[];
elseif sc==3
   k=0:2;
   ip=pi/2 + 2*pi/3 + 2*k*pi/3;
   A1=sqrt(2/3)*cos(ip);
   A2=sqrt(2/3)*sin(ip);
   A3=sqrt(2/3)*sqrt(1/2)*ones(1,3);
   A=[A1;A2;A3]';
   sX=X*A;
   B=A\eye(3);
   oX=sX*B;
   sim=sX;
   sim(:,3)=[];
else
   sim=X;
   sX=[];
   oX=[];
   error('Data set has not the correct dimension! Required is a 4d or 3d data set.');
end


