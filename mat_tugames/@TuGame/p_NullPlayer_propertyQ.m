function NPPQ=p_NullPlayer_propertyQ(clv,x,tol)
% NULLPLAYER_PROPERTYQ verifies if x satisfies the null player property 
% using Matlab's PCT.
%
%  Usage: NPPQ=clv.p_COV_propertyQ(x,tol)
%
% Define variables:
%  output: Fields
%  NppQ     -- Returns true (1), if x satisfies the null player property, otherwise
%              false (0).
%  snp      -- Set of null players.
%  pzp      -- Players with zero payoff.
%
%  input:
%  clv      -- TuGame class object.
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
%   12/27/2020        1.9             hme
%


N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;

if nargin<3
 tol=10^6*eps;
end
pl=1:n;
S=1:N;
npQ=false(1,n);
NppQ=false;
parfor k=1:n
    a=bitget(S,k)==1;
    Swk=S(a==0);
    Sk=S(a);
    vS=[0,v(Swk)];
    npQ(k)=all(abs(v(Sk)-vS)<10^6*eps);
end    
NpQ=any(npQ);
snp=pl(npQ);
nnp=pl(npQ==0);
nlp=pl(x(pl)~=0);
if NpQ==1
   zpQ=abs(x(snp))<10^6*eps;
   pzp=J(zpQ);
   NppQ=any(zpQ);
elseif any(nlp) %% contrapositive
   slc=ismember(nnp,nlp);
   NppQ=all(nnp(slc)==nlp);
   pzp=pl(x(pl)==0);
else
   NppQ=false;
   pzp=[];
end    

NPPQ=struct('NppQ',NppQ,'snp',snp,'pzp',pzp);

