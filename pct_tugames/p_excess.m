function ex=p_excess(v,x)
% P_EXCESS computes the excess vector of game v w.r.t. x.
% using Matlab's PCT.
%
% Poor Performance, use excess() instead.
%
%
% Usage: ex=p_excess(v,x)
% Define variables:
% output:
% ex         -- excess vector of game v w.r.t. x.
%
% input: 
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- payoff vector of size(1,n).
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/19/2010        0.1 beta        hme
%   06/10/2012        0.2 beta        hme
%   05/05/2014        0.5             hme
%                


% Computing the excess vector w.r.t. x.
N=length(v); n=length(x);
S=1:N;
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end
if poolsize < 6
  it=0:-1:1-n;
  PlyMat=x*rem(floor(S(:)*pow2(it)),2)';
  ex=v-PlyMat;
else
  PlyMat=zeros(n,N);
  parfor i = 1:n, PlyMat(i,:) = bitget(S,i); end
  clear S;
  spmd
    pym=codistributed(PlyMat);
    cv=codistributed(v);
    PlyMat=x*pym;
    ex1=cv-PlyMat;
  end
  ex=gather(ex1);
end
