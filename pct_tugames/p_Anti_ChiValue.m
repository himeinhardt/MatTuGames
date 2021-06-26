function chiv=p_Anti_ChiValue(v)
% P_ANTI_ChiVALUE computes the anti-chi-value of a TU-game v using MATLAB's PCT. This is a generalized 
% anti Tau value.
%
% Resource: Notes on a new compromise value: The Chi-Value, by G. Bergantinos & J. Masso (1996).
%
% Usage: chiv=p_Anti_ChiValue(v)
% Define variables:
%  output:
%  chiv      -- The anti Chi-Value of a TU-Game v.
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
%   04/07/2021        1.9             hme
%                


if nargin<1
    error('At least the game must be given!');  
else
    N=length(v);
    [~, n]=log2(N);  
    if (2^n-1)~=N, error('Game has not the correct size!'); end
end


xiv=zeros(1,n);

[g,bv,lv]=p_Anti_GenGap(v);


if g(N)==0
   chiv=bv;
elseif g(N)>0
   sumlv=lv*ones(n,1);
   if sumlv==0
     warning('Xi:ZeroConc','The sum of the anti-concession vector is zero. Assuming veto-rich game!');
     [~, lvp]=p_veto_rich_players(v);
     J=1:n;
     slcvp=J(lvp);
     lgtv=length(slcvp);
     if lgtv>0
        chiv=lvp/lgtv;
     else
     error('Game has no veto player! No Anti-Chi-Value computed. Sorry!');
     end  
   else
     chiv=bv+(g(N)*lv)/sumlv;
   end
else
  error('Game is not anti-essential! No Anti-Chi-Value computed. Sorry!');
end
