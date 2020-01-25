function ctv2=critical_value2(v)
%CRITICAL_VALUE2 computes a critical value w.r.t. the strong epsilon-core.
%Source: Maschler, Peleg and Shapley (1979, p. 313).
% 
% Usage: ctv2=critical_value2(v)
%
%
% Define variables:
%  output:
%  ctv2     -- Critical epsilon value w.r.t. the strong epsion-core.
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/14/2012        0.2 beta        hme
%   09/18/2012        0.2             hme
%   10/27/2012        0.3             hme
%
    
r=reasonable_outcome(v);
N=length(v);
[~, n]=log2(N);
S=1:N-1;

Rm=r(1); for k=2:n Rm=[Rm r(k) Rm+r(k)]; end
Rm(N)=[];
sr=sum(r);
srv=sr-Rm;
ctv2=max(v(S)-v(N)+srv);
