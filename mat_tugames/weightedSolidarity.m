function [w_sd, sd_uS]=weightedSolidarity(v,w_vec)
% WEIGHTEDSOLIDARITY computes the weighted Solidarity value of a TU-game v.
% Needs some time to complete for n>10.
% 
% Usage: [w_sd sd_uS]=weightedSolidarity(v,w_vec)
% Define variables:
%  output:
%  w_sd     -- The weighted Solidarity value of a TU-game v.
%  sd_us    -- The Solidarity value matrix of unanimity_games.
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
%   07/19/2013        0.4             hme
%                
if min(w_vec)<=0
  error('The vector of weights cannot contain zero or negative components!')
end
N=length(v);
[~, n]=log2(N);
[u_coord, sutm]=basis_coordinates(v);
k=1:n;
sd_uS=zeros(N,n);

for S=1:N;
    uw=zeros(1,n);
    clS=bitget(S,k)==1;
    sum_wS=clS*w_vec';
    plS=k(clS);
    uw(plS)=w_vec(plS);
    sd_uS(S,:)=(uw/sum_wS)*sutm(N,S);
end
w_sd=u_coord*sd_uS;

