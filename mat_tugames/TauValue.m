function tauv=TauValue(v)
% TAU_VALUE computes the tau-value of a TU-game v.
%
% Usage: tauv=TauValue(v)
% Define variables:
%  output:
%  tauv     -- The Tau-Value of a TU-Game v.
%
%  input:
%  v        -- A TU-Game of length 2^n-1.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/10/2010        0.1 beta        hme
%   05/18/2012        0.2 beta        hme
%   10/27/2012        0.3             hme
%                


if nargin<1
    error('At least the game must be given!');  
else
    N=length(v);
    [~, n]=log2(N);  
    if (2^n-1)~=N, error('Game has not the correct size!'); end
end


tauv=zeros(1,n);

[g bv lv]=Gap(v);


if g(N)==0
   tauv=bv;
elseif g(N)>0
   sumlv=lv*ones(n,1);
   if sumlv==0
     warning('Tau:ZeroConc','The sum of the concession vector is zero. Assuming simple game!');
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
     tauv=bv-(g(N)*lv)/sumlv;
   end
else
  error('Game is not quasi-balanced! No Tau-Value computed. Sorry!');
end
