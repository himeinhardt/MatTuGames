function [w_sh sh_uS]=weightedShapley(v,w_vec);
% weightedSHAPLEY computes the weighted Shapley-value of a TU-game v.
%
% Usage: [w_sh sh_uS]=weightedShapley(v,w_vec)
% Define variables:
%  output:
%  w_sh     -- The weighted Shapley-value of a TU-game v.
%  sh_us    -- The Shapley value matrix of unanimity_games.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%  w_vec    -- A vector of positive weights.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/19/2010        0.1 beta        hme
%   05/18/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%   06/25/2013        0.4             hme
%                
if min(w_vec)<=0
  error('The vector of weights cannot contain zero or negative components!')
end
N=length(v);
[~, n]=log2(N);
[u_coord sutm]=unanimity_games(v);
k=1:n;
sh_uS=zeros(N,n);

for S=1:N;
    uw=zeros(1,n);
    clS=bitget(S,k)==1;
    sum_wS=clS*w_vec';
    plS=k(clS);
    uw(plS)=w_vec(plS);
    sh_uS(S,:)=uw/sum_wS;
end
w_sh=u_coord*sh_uS;

