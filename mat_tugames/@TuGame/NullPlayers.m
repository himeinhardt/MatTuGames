function NP=NullPlayers(clv,tol)
% NULLPLAYERS returns the list of null players of game v. 
%
%  Usage: NP=clv.NullPlayers(tol)
%
% Define variables:
%  output: Fields
%  NP       -- Set of null players.
%
%  input:
%  clv      -- TuGame class object.
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
%   12/27/2020        1.9             hme
%

if nargin<2 
 tol=10^6*eps;    
end
ize;
n=clv.tuplayers;
v=clv.tuvalues;
pl=1:n;
S=1:N;
npQ=false(1,n);
NppQ=false;
for k=1:n
    a=bitget(S,k)==1;
    Swk=S(a==0);
    Sk=S(a);
    vS=[0,v(Swk)];
    npQ(k)=all(abs(v(Sk)-vS)<tol);
end    
J=1:n;
NP.lnpl=npQ;
NP.npl=J(npQ);

