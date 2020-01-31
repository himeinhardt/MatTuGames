function LED=LED_propertyQ(clv,x,tol)
%LED_PROPERTYQ checks wheter the solution x satisfies large excess difference property.
%
% Source:  H. I. Meinhardt. The Modiclus Reconsidered. Technical report, Karlsruhe Institute of Technology (KIT), Karlsruhe, Germany,
%          2018. URL http://dx.doi.org/10.13140/RG.2.2.32651.75043.
%
%          Meinhardt (2018), "Analysis of Cooperative Games with Matlab and Mathematica".
%
%          Sudhoelter (1997), The modified nucleolus: Properties and axiomatizations. International Journal of Game Theory, 26
%          (2):147–182, Jun 1997. ISSN 1432-1270. doi: 10.1007/BF01295846. URL https://doi.org/10.1007/BF01295846.
%
%
% Usage: LED=clv.LED_propertyQ(x,tol) 
%
% Define variables:
%
%  Output structure variables:
%  propQ    -- Returns true (1) whenever the solution has large excess difference,
%              otherwise false (0).
%  lexd     -- Retruns the value of the large excess difference.
%  x        -- Replicates the input variable x.
%  lcl      -- Returns coalition in integer representation having the largest excess.
%  sbexcl   -- Returns coalition in integer representation having the smallest bi-excess.
%
%  Input:
%  clv      -- TuGame class object.
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
tol=-tol;
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

if N==1
   LED.propQ=true;
   LED.lexd=inf; %% inf - (-inf)
   LED.x=x;
  return
end

dv=clv.dual_game();
ex_v=clv.excess(x);
dfv=v-dv;
dfv(end)=[];
ex_v(end)=[];
[mex_v,lcl]=max(ex_v);
[mdfv,sbexcl]=min(dfv);
ldf=mdfv-mex_v;
ldQ=ldf>=tol;
LED.propQ=ldQ;
LED.lexd=ldf;
LED.x=x;
LED.lcl=lcl;
LED.sbexcl=sbexcl;
