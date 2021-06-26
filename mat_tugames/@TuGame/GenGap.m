function [g bv lv]=GenGap(clv)
% GENGAP computes the generalized gap function from game v.
%
% Usage: [g bv lv]=clv.GenGap()
% Define variables:
%  output:
%  g        -- The generalized gap function of game v. A vector of length 2^n-1.
%  bv       -- The generalized upper vector/payoff of game v.
%  lv       -- The generalized lower/concession vector of game v.
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
%   08/27/2020        1.9             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

% upper vector
S=1:N;
bv=zeros(1,n);
k=1:n;
a=cell(n,1);
for ii=1:n
   a{ii}=bitget(S,ii)==1;
   si=S(a{ii});
   fe=si(1);
   nsi=bitset(si,ii,0);
   si(1)=[];
   nsi(1)=[];
   bv(ii)=max(max(v(si)-v(nsi)),v(fe));
end


% Computing the gap function w.r.t. v.
g=zeros(1,N); % the gap vector w.r.t. v.
Bm=bv(1); for ii=2:n, Bm=[Bm bv(ii) Bm+bv(ii)]; end
g=Bm-v;

lv=zeros(1,n);
% generalized concession vector
for ii=1:n
   lv(ii)=min(g(a{ii}));
end


