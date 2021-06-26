function [uv,mv]=p_UpperPayoff(clv)
% UTOPIAPAYOFF computes the utopia and minimum claim vector of game v using MATLAB's PCT.
%
% Usage: [uv mv]=clv.p_UpperPayoff()
% Define variables:
%  output:
%  uv       -- The upper vector/payoff of game v.
%  mv       -- The minimum claim (disagreement) vector of game v.
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
%   08/27/2020        1.9             hme
%                


v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
% upper vector
k=1:n;
Si=bitset(N,k,0);
uv=v(N)-v(Si);

UV=uv(1); for k=2:n,UV=[UV uv(k) UV+uv(k)]; end


mv=zeros(1,n);
S=1:N;
parfor k=1:n
  Swk=S(bitget(S,k)==1);
  Snk=Swk-2^(k-1);
  if Snk(1)==0
     Snk(1)=[];
     av=[0,UV(Snk)];
  else
     av=UV(Snk);
  end
  mv(k)=max(v(Swk)-av);
end


