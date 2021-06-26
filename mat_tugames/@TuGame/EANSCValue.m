function ean=EANSCValue(clv)
% EANSCValue computes the Equal Allocation of Non-Separable Contribution/Cost Value 
%
% Source: Moulin, H. (1985). The separability axiom and equal-sharing methods, 
% Journal of Economic Theory 36(1): 120-148. 
%
% Usage: ean=clv.EANSCValue()
%
% Define variables:
%  output:
%  ean       -- Equal Allocation of Non-Separable Contribution/Cost Value 
%               
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
%   05/25/2020        1.9             hme
%
v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if N==1
  ean=v;return;
 else
end
k=1:n;
Ni=N-2.^(k-1);
% upper payoff/Separable Contribution to the Grand Coalition
bv=v(N)-v(Ni);
% Non-Separable Contribution/Cost (non-negative/non-positive)
nsc=v(N)-sum(bv);
% EANSC Value
ean=bv+nsc/n;
