function tauv=TauValue(clv)
% TAU_VALUE computes the tau-value of a TU-game v.
%
% Usage: tauv=TauValue(v)
% Define variables:
%  output:
%  tauv     -- The Tau-Value of a TU-Game v.
%
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
%   10/30/2012        0.3             hme
%                


N=clv.tusize;
n=clv.tuplayers;


tauv=zeros(1,n);

[g, bv, lv]=Gap(clv);


if g(N)==0
   tauv=bv;
elseif g(N)>0
   sumlv=lv*ones(n,1);
   if sumlv==0
     warning('Tau:ZeroConc','The sum of the concession vector is zero. Assuming simple game!');
     [~, lvp]=veto_players(clv);
     J=1:n;
     slcvp=J(lvp);
     lgtv=length(slcvp);
     if lgtv>0
        tauv=lvp/lgtv;
     else
     error('Simple game has no veto player! No Tau-Value computed. Sorry!');
     end  
   else
     tauv=bv-(g(N)*lv)/sumlv;
   end
else
  error('Game is not quasi-balanced! No Tau-Value computed. Sorry!');
end
