function ctv_star=critical_value_star(clv)
%CRITICAL_VALUE_STAR computes the critical value at which the strong
%epsilon core just contains the intersection of the imputation and reasonable set.
%Source: Maschler, Peleg and Shapley (1979, p. 320).
% 
% Usage: ctv_star=critical_value_star(clv)
%
%
% Define variables:
%  output:
%  ctv_star -- Critical epsilon value w.r.t. the strong epsion-core.
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

r=reasonable_outcome(clv);
Rm=r(1); for k=2:n Rm=[Rm r(k) Rm+r(k)]; end
Rm(N)=[];
sr=sum(r);
srv=sr-Rm;

ctv_star=max(min(v(S)-Vm,v(S)-v(N)+srv));
