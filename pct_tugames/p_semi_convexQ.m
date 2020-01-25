function scvQ=p_semi_convexQ(v)
% P_SEMI_CONVEXQ checks whether the game v is semi-convex using Matlab's PCT.
%
%
% Usage: scvQ=p_semi_convexQ(v)
% Define variables:
%  output:
%  ck       -- Returns 1 (true) whenever the game v is semi-convex, 
%              otherwise 0 (false).
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/19/2011        0.1 alpha        hme
%   10/27/2012        0.3              hme
%   05/16/2014        0.5              hme
%                


N=length(v);
[~, n]=log2(N);
lrQ=cell(n,1);

[gv bv lv]=p_Gap(v);

k=1:n;
sC=bitset(0,k);
S=1:N;


grZ=false(n,1);
grQ=false(n,1);
scvQ=false(n,1);

parfor i=1:n
  Si=bitget(S,i)==1;
  grZ(i)=gv(sC(i))>=0;
  lrQ{i}=gv(Si)>=gv(sC(i));
  grQ(i)=all(lrQ{i});
  scvQ(i)=grQ(i).*grZ(i);
end


scvQ=all(scvQ);
