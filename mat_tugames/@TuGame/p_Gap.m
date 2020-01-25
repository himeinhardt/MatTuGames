function [g, bv, lv]=p_Gap(clv)
% P_GAP computes the gap function from game v using Matlab's PCT.
%
% Usage: [g bv lv]=p_Gap(clv)
%
% Define variables:
%  output:
%  g        -- The gap function of game v. A vector of length 2^n-1.
%  bv       -- The upper vector/payoff of game v.
%  lv       -- The lower/concession vector of game v.
%
%  input:
%  clv      -- TuGame class object.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   10/30/2012        0.3              hme
%                


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
Si=clv.tuSi;
% upper vector
bv=v(N)-v(Si);


% Computing the gap function w.r.t. v.
S=1:N;
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=Bm-v;

lv=zeros(1,n);

% concession vector
parfor i=1:n
   a_i=bitget(S,i)==1;
   lv(i)=min(g(a_i));
end

