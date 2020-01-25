function ex=p_excess(clv,x)
% P_EXCESS computes the excess vector of game v w.r.t. x.
% using Matlab's PCT.
%
%
% Usage: ex=excess(clv,x)
%
% Define variables:
% output:
% ex         -- excess vector of game v w.r.t. x.
%
% input: 
%  clv        -- TuGame class object.
%  x        -- payoff vector of size(1,n).
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

% Computing the excess vector w.r.t. x.
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
