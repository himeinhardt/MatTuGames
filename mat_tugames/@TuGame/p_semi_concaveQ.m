function scvQ=p_semi_concaveQ(clv)
% P_SEMI_CONCAVEQ checks whether the game v is semi-concave using Matlab's PCT.
%
%
% Usage: scvQ=clv.semi_concaveQ()
% Define variables:
%  output:
%  scvQ       -- Returns 1 (true) whenever the game v is semi-concave, 
%                otherwise 0 (false).
%
%  input:
%  clv        -- TuGame class object.
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


N=clv.tusize;
n=clv.tuplayers;
lrQ=cell(n,1);

gv=clv.p_Gap();

k=1:n;
sC=bitset(0,k);
grQ=false(1,n);
scvQ=false(1,n);
S=1:N;
parfor i=1:n
  Si=bitget(S,i)==1;
  grZ(i)=gv(sC(i))<=0;
  lrQ{i}=gv(Si)<=gv(sC(i));
  grQ(i)=all(lrQ{i});
  scvQ(i)=grQ(i).*grZ(i);
end

scvQ=all(scvQ);

