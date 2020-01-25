function ad_vl=p_ADvalue(clv,cs)
% P_AD_VALUE computes the Aumann-Dreze value w.r.t. 
% a coalition structure using Matlab's PCT.
%
% Usage: ad_vl=p_ADvalue(clv,cs)
%
% Define variables:
%  output:
%  ad_vl    -- The Aumann-Dreze value w.r.t to a coalition structure.
%
%  input:
%  clv      -- TuGame class object.
%  cs       -- A coalition structure like [3 4 8]
%              for {[1,2],[3],[4]}. Must be a partition.    
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/26/2013        0.4             hme
%   05/12/2014        0.5             hme
%

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
if isa(clv,'TuVal')
   ptn = clv.tu_ptn;
elseif isa(clv,'p_TuVal')
   ptn = clv.tu_ptn;
else
   ptn='';
end

if nargin < 2
   if isempty(ptn)
       error('A game and coalition structure must be given!');
   else
      cs = ptn;
   end
end

J=1:n;
int=0:-1:1-n;
csm=rem(floor(cs(:)*pow2(int)),2)==1;
lcs=length(cs);
ad_vl=zeros(1,n);
iad=cell(1,lcs);
idx=cell(1,lcs);
parfor k=1:lcs
    Tk=csm(k,:);
    idx{k}=J(Tk);
    sbv=SubGame(v,cs(k));
    iad{k}=ShapleyValue(sbv);
end
for j=1:lcs 
  ad_vl(idx{j})=iad{j};
end
    
