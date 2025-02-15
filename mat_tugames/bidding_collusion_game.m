function v=bidding_collusion_game(bd)
% BIDDING_COLLUSION_GAME determines a bidding collusion game from a bid vector.
%
% Usage: v=bidding_collusion_game(bd)
%
% Define variables:
%  output:
%  v         -- The bidding collusion game of lenght 2^n-1.
%
%  input:
%  bd        -- A bidding vector of length n.
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/08/2025        1.9.2             hme
%    
    
%[bd,idx] = sort(bd,'ascend');
n=length(bd);
N=2^n-1;
v=zeros(1,N);
bm=max(bd);
v(N)=bm;
k=1:n;
for sS=1:N-1
    v(sS)=bm-max(bd(bitget(sS,k)==0));
end
