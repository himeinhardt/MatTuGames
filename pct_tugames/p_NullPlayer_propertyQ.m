function NPPQ=p_NullPlayer_propertyQ(v,x,tol)
% NULLPLAYER_PROPERTYQ verifies if x satisfies the null player property 
% using Matlab's PCT.
%
%  Usage: COV=p_COV_propertyQ(v,x,m,t,str)
%
% Define variables:
%  output: Fields
%  NppQ     -- Returns true (1), if x satisfies the null player property, otherwise
%              false (0).
%  snp      -- Set of null players.
%  pzp      -- Players with zero payoff.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
%              

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/18/2015        0.7             hme
%

if nargin<3
 tol=10^6*eps;
end

N=length(v);
[~, n]=log2(N);

S=1:N;
npQ=false(1,n);

parfor k=1:n
    a=bitget(S,k)==1;
    Swk=S(a==0);
    Sk=S(a);
    vS=min(v(Swk));
    npQ(k)=all(abs(v(Sk)-vS)<10^6*eps);
end    

NpQ=any(npQ);
J=1:n;
snp=J(npQ);

if NpQ==1
   zpQ=abs(x(snp))<10^6*eps;
   pzp=J(zpQ);
   NppQ=any(zpQ);
else
   NppQ=false;
   pzp=[];
end    

NPPQ=struct('NppQ',NppQ,'snp',snp,'pzp',pzp);

