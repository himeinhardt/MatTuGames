function [g bv lv]=p_Anti_Gap(clv)
% P_NTI_GAP computes the anti-gap function from game v using Matlab's PCT.
%
% Usage: [g bv lv]=clv.p_Anti_Gap()
% Define variables:
%  output:
%  g        -- The anti-gap function of game v. A vector of length 2^n-1.
%  bv       -- The anti-upper vector/payoff of game v.
%  lv       -- The anti-lower/concession vector of game v.
%
%  input:
%  clv        -- TuGame class object.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/07/2021        1.9              hme
%                



v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
Si=clv.tuSi;
% anti-upper vector
bv=v(N)-v(Si);

% Computing the anti-gap function w.r.t. v.
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=v-Bm;

lv=zeros(1,n);

% anti-concession vector
S=1:N;
parfor i=1:n
   a_i=bitget(S,i)==1;
   lv(i)=min(g(a_i));
end

lv=gather(lv);
