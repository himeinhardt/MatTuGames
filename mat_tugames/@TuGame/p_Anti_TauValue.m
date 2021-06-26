function tauv=p_Anti_TauValue(clv)
% P_ANTI_TAUVALUE computes the anti-tau-value of a TU-game v.
% Using Matlab's PCT.
%
% Usage: tauv=clv.p_Anti_TauValue()
% Define variables:
%  output:
%  tauv     -- The Anti-Tau-Value of a TU-Game v.
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
%   04/07/2021        1.9              hme
%                


v=clv.tuvalues;
n=clv.tuplayers;
N=clv.tusize;

tauv=zeros(1,n);

[g,bv,lv]=p_Anti_Gap(clv);


if g(N)==0
   tauv=bv;
elseif g(N)>0
   sumlv=lv*ones(n,1);
   if sumlv==0
     warning('Tau:ZeroConc','The sum of the anti-concession vector is zero. Assuming simple game!');
     [~, lvp]=veto_players(v);
     J=1:n;
     slcvp=J(lvp);
     lgtv=length(slcvp);
     if lgtv>0
        tauv=lvp/lgtv;
     else
     error('Simple game has no veto player! No Tau-Value computed. Sorry!');
     end  
   else
     tauv=bv+(g(N)*lv)/sumlv;
   end
else
  error('Game is not anti-quasi-balanced! No Tau-Value computed. Sorry!');
end
