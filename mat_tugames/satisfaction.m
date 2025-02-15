function SATIS=satisfaction(v,x)
% SATISFACTION computes the satisfaction of a partition and of its
% anti partition, if those exist.
%
% Usage: SATIS=satisfaction(v,x)
%
% Define variables:
%  output:
%  FC       -- This structure element returns the satisfaction of a partition.
%  AFC      -- This structure element returns satisfaction of an anti partition.    
%  ptn      -- This structure element returns the set of most effective coalitions 
%              with largest excess (best coalitions).
%  ptnQ     -- Structure element ptnQ returns 1 (true) if ptn forms a partition of the player set, 
%              otherwise 0 (false).
%  aptn     -- The anti partition if pnt is a partition, otherwise none.
%  CIS      -- The center of imputation set.
%  EC       -- The excess of the whole player partition.
%  OB       -- The objector set.
%  setD     -- Coalitions other than the grand coalitons with largest nontrivial excess.
%  lex      -- Largest nontrivial excess.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- payoff vector of size(1,n).
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/04/2013        0.5             hme
%   01/25/2022        1.9.1           hme
%

N=length(v);
[f1, n]=log2(N);    
    
PartQ=FindPartition(v,x);
if PartQ.ptnQ==1
    lgp=length(PartQ.ptn);
    algp=length(PartQ.aptn);
    vS=v(PartQ.ptn);
    vaS=v(PartQ.aptn);
    FC=(v(N)-sum(vS))/lgp;
    AFC=((algp-1)*v(N)-sum(vaS))/algp;
else
  msg='No partition found, no satisfaction can be computed!';
  warning(msg);
  FC='none';
  AFC='none';
end    
k=1:n;
sC=2.^(k-1);
EC=(v(N)-sum(v(sC)))/n;
CIS=v(sC) + EC;
exc=excess(v,CIS);
exc(end)=[];
obj=find(exc>-EC);
SATIS=struct('FC',FC,'AFC',AFC,'ptn',PartQ.ptn,'ptnQ',PartQ.ptnQ,'aptn',PartQ.aptn,'CIS',CIS,'EC',-EC,'OB',obj,'setD',PartQ.setD,'lex',PartQ.lex);
