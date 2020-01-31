function [eqtQ,eqt,sbs]=p_equal_treatmentQ(v,x,tol)
% P_EQUAL_TREATMENTQ checks if the vector x satisfies the equal
% treatment property (ETP) using Matlab's PCT.
%
% Usage: [eqtQ,eqt]=p_equal_treatmentQ(v,x,tol)
%
% Define variables:
%  output:
%  eqtQ       -- Returns true (1) if the solution x satisfies ETP,
%                otherwise false (0).   
%  eqt        -- Returns an array of ones (true )and/or zeros
%                (false) for each investigated pair of substitutes.
%  sbs        -- A matrix of substitutes. 
%                Shows in each row the pair that are substitutes.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n) (optional)
%  tol      -- Tolerance value. Its default value is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/06/2015        0.6             hme
% 


if nargin < 3
    tol=10^6*eps;
end    

sbs_v=p_substitutes(v);
if isempty(sbs_v)==0
   ls=size(sbs_v,1);
   eqt=false(1,ls);
   for k=1:ls
       y=x(sbs_v(k,:));
       eqt(k)=abs(y(1)-y(2))<tol;
   end    
   eqtQ=all(eqt);
else
   eqtQ=false;
   eqt=[];
end    
