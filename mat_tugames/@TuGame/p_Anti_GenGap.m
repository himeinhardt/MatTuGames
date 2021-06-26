function [g bv lv]=p_Anti_GenGap(clv)
% P_ANTI_GENGAP computes the anti-generalized gap function from game v using Matlab's PCT.
%
% Usage: [g bv lv]=clv.p_Anti_GenGap()
% Define variables:
%  output:
%  g        -- The anti-generalized gap function of game v. A vector of length 2^n-1.
%  bv       -- The generalized upper vector/payoff of game v.
%  lv       -- The anti-generalized lower/concession vector of game v.
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
%   04/07/2021        1.9             hme
%                



v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% anti-gneralized upper vector
S=1:N;
bv=zeros(1,n);
k=1:n;
a=cell(n,1);
%Mgc=AllMarginalContributions(v);
parfor ii=1:n
   a{ii}=bitget(S,ii)==1;
   si=S(a{ii});
   fe=si(1);
   nsi=bitset(si,ii,0);
   si(1)=[];
   nsi(1)=[];
   bv(ii)=min(min(v(si)-v(nsi)),v(fe));
end


% Computing the anti-generalized gap function w.r.t. v.
g=zeros(1,N); % the gap vector w.r.t. v.
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=v-Bm;

lv=zeros(1,n);
% anti-generalized concession vector
parfor ii=1:n
   lv(ii)=min(g(a{ii}));
end


