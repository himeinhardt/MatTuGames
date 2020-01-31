function ad_vl=ADvalue(v,cs)
% AD_VALUE computes the Aumann-Dreze value w.r.t. 
% a coalition structure.
%
% Usage: ad_vl=ADvalue(v,cs)
%
% Define variables:
%  output:
%  ad_vl    -- The Aumann-Dreze value w.r.t to a coalition structure.
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
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
%   07/17/2013        0.4             hme
%
if nargin < 2
   error('A game and coalition structure must be given!'); 
elseif nargin==2
    N=length(v);
    [~, n]=log2(N);
    if (2^n-1)~=N
      error('Game has not the correct size!');
    end
    if iscell(cs)
       cs=clToMatlab(cs);
    end
end    
    
J=1:n;
int=0:-1:1-n;
csm=rem(floor(cs(:)*pow2(int)),2)==1;
lcs=length(cs);
ad_vl=zeros(1,n);
for k=1:lcs
    Tk=csm(k,:);
    idx=J(Tk);
    sbv=SubGame(v,cs(k));
    ad_vl(idx)=ShapleyValue(sbv);
end    
