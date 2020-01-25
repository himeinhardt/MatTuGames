function ctv2=critical_value2(clv)
%CRITICAL_VALUE2 computes a critical value w.r.t. the strong epsilon-core.
%Source: Maschler, Peleg and Shapley (1979, p. 313).
% 
% Usage: ctv2=critical_value2(clv)
%
%
% Define variables:
%  output:
%  ctv2     -- Critical epsilon value w.r.t. the strong epsion-core.
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
    
r=reasonable_outcome(clv);
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
S=1:N-1;

Rm=r(1); for k=2:n Rm=[Rm r(k) Rm+r(k)]; end
Rm(N)=[];
sr=sum(r);
srv=sr-Rm;
ctv2=max(v(S)-v(N)+srv);
