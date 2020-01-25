function ctv1=critical_value1(clv)
%CRITICAL_VALUE1 computes the biggest gain that any group of players can
%ensure while forming a coalition.
%Source: Maschler, Peleg and Shapley (1979, p. 306).
% 
% Usage: ctv1=critical_value1(clv)
%
%
% Define variables:
%  output:
%  ctv1     -- Critical epsilon value w.r.t. the strong epsion-core.
%
%  input:
%  clv        -- TuGame class object.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/29/2012        0.3             hme
%

    
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
S=1:N-1;
k=1:n;
vi=v(bitset(0,k));
Vm=vi(1); for ii=2:n, Vm=[Vm vi(ii) Vm+vi(ii)]; end
Vm(N)=[];

ctv1=max(v(S)-Vm);
