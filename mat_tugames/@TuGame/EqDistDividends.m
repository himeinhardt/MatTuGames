function shv=EqDistDividends(clv)
% EqDistDividends computes the equally distributed dividends of coalitions to which players belong, i.e., Shapley-value, of a TU-game v. 
%
% Usage: shv=EqDistDividends(v)
% Define variables:
%  output:
%  shv      -- The Shapley-value of a TU-game v.
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   04/22/2024        1.9.1           hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

hd=clv.harsanyi_dividends(); 
it=0:-1:1-n;
sS=1:N;
mat=rem(floor(sS(:)*pow2(it)),2);
% cardinality of coalitions.
cd=sum(mat,2)';
% per capita dividends.
chd=hd./cd; 
for ii=1:n
    cbii=bitget(sS,ii)==1;
    shv(ii)=sum(chd(cbii));
end    
