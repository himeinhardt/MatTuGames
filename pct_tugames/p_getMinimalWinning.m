function [mW,wC]=p_getMinimalWinning(v)
% P_GETMINIMALWINNING computes from a simple game
% the minimal winning coalitions using Matlab's PCT.
%
% Usage: [mW,wC]=p_getMinimalWinning(v);
%
% Define variables:
%  output:
%  mW       -- The list/vector of minimal winning coalitions.
%  wC       -- The list of winning coalitions. 
%  input:
%  v        -- A simple TU game array of length 2^n-1.


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   09/29/2020        1.9             hme
%                    

N=length(v);
[~, n]=log2(N);
S=1:N;
v=logical(v);
wC=S(v);
lwC=length(wC);

parfor k=1:lwC
  cl=wC(k);
  sS=SubSets(cl,n);
  [vl, idx]=max(v(sS));
  smS=sS(idx);
  if smS==wC(k)
      mWS(k)=cl;
  end
end

mW=mWS(logical(mWS>0));
