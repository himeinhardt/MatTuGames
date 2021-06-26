function atauv=Anti_TauValue(v)
% ANTI_TAU_VALUE computes the anti-tau-value of a TU-game v.
%
% Usage: atauv=Anti_TauValue(v)
% Define variables:
%  output:
%  atauv     -- The Anti-Tau-Value of a TU-Game v.
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
%   03/30/2021        1.9             hme
%                


if nargin<1
    error('At least the game must be given!');  
else
    N=length(v);
    [~, n]=log2(N);  
    if (2^n-1)~=N, error('Game has not the correct size!'); end
end


tauv=zeros(1,n);

[g bv lv]=Anti_Gap(v);

if g(N)==0
   atauv=bv;
elseif g(N)>0
   sumlv=lv*ones(n,1);
   if sumlv==0
     warning('Tau:ZeroConc','The sum of the anti-concession vector is zero. Assuming simple game!');
     [~, lvp]=veto_players(v);
     J=1:n;
     slcvp=J(lvp);
     lgtv=length(slcvp);
     if lgtv>0
        atauv=lvp/lgtv;
     else
     error('Simple game has no veto player! No Anti-Tau-Value computed. Sorry!');
     end  
   else
     atauv=bv+(g(N)*lv)/sumlv;
   end
else
  error('Game is not anti-quasi-balanced! No Anti-Tau-Value computed. Sorry!');
end
