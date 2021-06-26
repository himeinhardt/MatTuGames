function DP=p_DummyPlayers(clv,tol)
% P_DUMMYPLAYERS returns the list of dummy players of game v using Matlab's PCT. 
%
%  Usage: DP=clv.p_DummyPlayers(tol)
%
% Define variables:
%  output: Fields
%  DP       -- Set of dummy players.
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
N=clv.tusize;
n=clv.tuplayers;
v=clv.tuvalues;
pl=1:n;
cli=2.^(pl-1);
S=1:N;
npQ=false(1,n);
NppQ=false;
parfor k=1:n
    a=bitget(S,k)==1;
    Swk=S(a==0);
    Sk=S(a);
    vS=[0,v(Swk)];
    npQ(k)=all(abs(v(Sk)-vS-v(cli(k)))<tol);
end    
J=1:n;
DP=J(npQ);


