function [uv,mv,utpQ]=p_UtopiaPayoff(clv)
% P_UTOPIAPAYOFF computes the utopia and minimum claim vector of game v using MATLAB's PCT.
%
% Usage: [uv mv,utpQ]=clv.p_UtopiaPayoff()
% Define variables:
%  output:
%  uv       -- The upper vector/payoff of game v.
%  mv       -- The minimum claim (concession) vector of game v.
%  utpQ     -- Returns a true (1), if the upper payoff is even the
%              utopia payoff, otherwise false (0).
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
%   08.27.2020        1.9             hme
%                



v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

% upper vector
k=1:n;
Si=bitset(N,k,0);
uv=v(N)-v(Si);

gp=clv.p_Gap();
giQ=gp(Si)>=0;
gnQ=gp(N)>0;

utpQ=giQ & gnQ;

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


