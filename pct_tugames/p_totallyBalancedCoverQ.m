function [tbcQ,v1]=p_totallyBalancedCoverQ(v,tol)
% P_TOTALLYBALANCEDCOVERQ checks whether the characteristic function v is equal to a 
% totally balanced cover using MATLAB's PCT. 
%
%
% Usage: [tbcQ,v1]=p_totallyBalancedCoverQ(v,tol)
% Define structure variables:
%  output:
%  tbcQ      -- Returns true (1) whenever the characteristic function v is equal to a totally balanced cover, otherwise false.
%  v1        -- Returns a totally balanced cover. This is a game.  
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/09/2022        1.9.1           hme
%
    
if nargin<2
 tol=10^8*eps;
end

N=length(v);
[~, n]=log2(N);
v1=zeros(1,N);
parfor S=1:N
 vS=SubGame(v,S);
 if length(vS)==1
    v1(S)=vS;
 elseif all(vS==0)==1
    v1(S)=vS(end); 
 else  
   [bcQ,bcv]=balancedCoverQ(vS,tol);
   v1(S)=bcv(end);
 end  
end
tbcQ=all(abs(v1-v)<=tol);
