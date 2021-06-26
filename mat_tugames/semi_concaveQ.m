function scvQ=semi_concaveQ(v)
% SEMI_CONCAVEQ checks whether the game v is semi-concave.
%
%
% Usage: scvQ=semi_concaveQ(v)
% Define variables:
%  output:
%  scvQ       -- Returns 1 (true) whenever the game v is semi-concave, 
%                otherwise 0 (false).
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
%   03/30/2021        1.9             hme
%                


N=length(v);
[~, n]=log2(N);
lrQ=cell(n,1);

gv=Anti_Gap(v);

k=1:n;
sC=bitset(0,k);
grQ=false(1,n);
scvQ=false(1,n);
S=1:N;
for i=1:n
  Si=bitget(S,i)==1;
  grZ(i)=gv(sC(i))>=0;
  lrQ{i}=gv(Si)>=gv(sC(i));
  grQ(i)=all(lrQ{i});
  scvQ(i)=grQ(i).*grZ(i);
end

scvQ=all(scvQ);

