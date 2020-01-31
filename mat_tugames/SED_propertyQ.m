function SED=SED_propertyQ(v,x,tol)
%SED_PROPERTYQ checks wheter the solution x satisfies small excess difference property.
%
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica"
%
% Usage: SED=SED_propertyQ(v,x,tol) 
%
% Define variables:
%
%  Output structure variables:
%  propQ    -- Returns true (1) whenever the solution has small excess difference,
%              otherwise false (0).
%  sexd     -- Retruns the value of the small excess difference.
%  x        -- Replicates the input variable x.
%  scl      -- Returns coalition in integer representation having the smallest excess.
%  lbexcl   -- Returns coalition in integer representation having the largest bi-excess.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n). Must be efficient.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%              (optional) 
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   03/01/2018        1.0             hme
%


if nargin<3
 tol=10^6*eps;
elseif nargin==3
 if ischar(tol) 
    error('Tolerance value is a string character! Must be a number!');
 end
end
N=length(v);
[~, n]=log2(N);
if N==1
   SED.propQ=true;
   SED.lexd=-inf; %% -inf - (inf)
   SED.x=x;
  return
end


dv=dual_game(v);
ex_v=excess(v,x);
dfv=v-dv;
dfv(end)=[];
ex_v(end)=[];
[mex_v,scl]=min(ex_v);
[mdfv,lbexcl]=max(dfv);
ldf=mdfv-mex_v;
ldQ=ldf<=tol;
SED.propQ=ldQ;
SED.sexd=ldf;
SED.x=x;
SED.scl=scl;
SED.lbexcl=lbexcl;
