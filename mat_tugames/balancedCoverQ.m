function [bcQ,v1,lS]=balancedCoverQ(v,tol)
% BALANCEDCOVERQ checks whether the characteristic function v is equal to a balanced cover.
%
%
% Usage: [bcQ,v1]=balancedCoverQ(v,tol)
% Define structure variables:
%  output:
%  bcQ      -- Returns true (1) whenever the characteristic function v is equal to a balanced cover, otherwise false.
%  v1       -- Returns a balanced cover. This is a game.  
%  lS       -- Weights (dual solution).
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
%   08/24/2022        1.9.1           hme
%
    
if nargin<2
 tol=10^8*eps;
end

N=length(v);
[~, n]=log2(N);
%try
%  [crQ,~,lS]=CddCoreQ(v,tol); % getting not the correct dual-solution!!
%catch
  [crQ,~,lS]=coreQ(v,tol);
%end	
sS=find(lS);
v1=v;
v1(N)=lS(sS)*v(sS)';
%v1(N)=lS*v'
bcQ=all(abs(v1-v)<=tol);
