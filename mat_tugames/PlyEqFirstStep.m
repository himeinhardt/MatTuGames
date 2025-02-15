function PEF=PlyEqFirstStep(th,w_vec,MmW)
% PLYEQFIRSTSTEP returns the last player that is equivalent to the first step.
%
% Source: Sudhoelter (1996), Star-shapedness of the kernel for homogeneous games.    
%
%
% Usage: PEF=PlyEqFirstStep(th,w_vec);
%
% Define field variables:
%  output:
%  lp       -- Returns last player equivalent to the first step.
%  fst      -- Returns first step.
%
%  input:
%  th       -- Threshold to pass a bill (positive number).
%  w_vec    -- Vector of sorted weights in descend order.
%  MmW      -- minimal winning coalitions in matrix form (optional).
%              Use function min_homogrep to retrieve correct format.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/10/2021        1.9             hme
%   05/24/2022        1.9.1           hme
%

w_svec=sort(w_vec,'descend');
if nargin < 3
   PPS=GetPlayersCharacter(th,w_svec);
else 
   PPS=GetPlayersCharacter(th,w_svec,MmW);
end	
idx=PPS.steps(1);
rs = max(PPS.steps(w_svec(idx)==w_svec(PPS.steps)));
if idx < rs
   PEF.lp=rs;
   PEF.fst=idx;
else
   %PEF.lp=[];
   PEF.lp=rs;
   PEF.fst=idx;
end	
