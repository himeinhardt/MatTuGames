function DP=DummyPlayers(v,tol)
% DUMMYPLAYERS returns the list of dummy players of game v. 
%
%  Usage: DP=DummyPlayers(v,tol)
%
% Define variables:
%  output: Fields
%  DP       -- Set of dummy players.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
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

N=length(v);
[~, n]=log2(N);
pl=1:n;
cli=2.^(pl-1);
S=1:N;
npQ=false(1,n);
NppQ=false;
for k=1:n
    a=bitget(S,k)==1;
    Swk=S(a==0);
    Sk=S(a);
    vS=[0,v(Swk)];    
    npQ(k)=all(abs(v(Sk)-vS-v(cli(k)))<tol);
end    
J=1:n;
DP.ldpl=npQ;
DP.dpl=J(npQ);

