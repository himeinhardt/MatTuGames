function ex=p_excess(v,x)
% P_EXCESS computes the excess vector of game v w.r.t. x.
% using Matlab's PCT.
%
%
% Usage: ex=excess(v,x)
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
%                

% Computing the excess vector w.r.t. x.
N=length(v); n=length(x);
ex=zeros(1,N); % the excesses of x.
S=1:N;
if matlabpool('size') < 6
  it=0:-1:1-n;
  PlyMat=rem(floor(S(:)*pow2(it)),2);
else
  PlyMat=zeros(N,n);
  parfor i = 1:n, PlyMat(:,i) = bitget(S,i); end
end
PayMat=PlyMat*x';
ex=(v-PayMat');
