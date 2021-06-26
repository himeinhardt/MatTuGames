function [sC,d2b]=SubCoalitions(Ar,n,csz)
% SUBCOALITIONS computes the power set (subsets) from an array/vector Ar. 
%
% Example: 
%     Set n=8;
%     Returns all subsets of set/coalition Ar=[3 5 6 8]:
%     sC=SubCoalitions(Ar,n)
%
% sC =
%
%     4    16    20    32    36    48    52   128   132   144   148   160   164   176   180.
%
% Usage: sC=SubCoalitions(Ar,n)
% Define variables:
%  output:
%  sC       -- Subsets of the vector Ar.
%
%  input:
%  Ar       -- A vector of players not larger than n, which represents a 
%              subcoalition. 
%  n        -- A positive number, which indicates the number of 
%              players in a Tu-game

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/22/2014        0.6             hme
%

if nargin < 3
   csz='';
end

Ar=sort(Ar);
k=1:n;
lA=length(Ar);
ic=2.^(Ar-1);
nS=2^lA-1;
S=1:nS;
it=0:-1:1-lA;
mS2=rem(floor(S(:)*pow2(it)),2);
ric=repmat(ic,[nS,1]);
sS=sum(mS2 .* ric,2);
if isempty(csz)==0
  sl=sum(mS2,2)==csz; 
  sC=sS(sl)';
else
  sC=sS;
end

if nargout == 2
  d2b=dec2bin(sC,n);
end


