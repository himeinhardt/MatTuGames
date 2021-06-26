function [mW,wC]=getMinimalWinning(clv)
% GETMINIMALWINNING computes from a simple game
% the minimal winning coalitions.
%
% Usage: [mW,wC]=clv.getMinimalWinning();
%
% Define variables:
%  output:
%  mW       -- The list/vector of minimal winning coalitions.
%  wC       -- The list of winning coalitions. 
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
%   09/29/2020        1.9             hme
%                    

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;
if strcmp(gt,'sv')
else
   error('Wrong game type!. Game must be a simple game!')
end

S=1:N;
v=logical(v);
wC=S(v);
lwC=length(wC);
sS=cell(lwC,1);

for k=1:lwC
  cl=wC(k);
  sS=SubSets(cl,n);
  [vl, idx]=max(v(sS));
  smS=sS(idx);
  if smS==wC(k)
      mWS(k)=cl;
  end
end

mW=mWS(logical(mWS>0));
