function [vp vp_vec]=p_veto_rich_players(v);
% VETO_RICH_PLAYERS returns a list of veto players for the TU-game v using MATLAB's PCT.
%
% Source: J. Arin and V. Feltkamp. The nucleolus and kernel of veto-rich transferable 
%         utility games." Int J of Game Theory, 26:61-73, 1997.
%
% Define variables:
%  output:
%  vp       -- The list of veto rich players for the TU-game v. It returns
%              a list of zeros whenever the game has no veto player.
%  vp_vec   -- The list of veto players in logical format.
%
%  input:
%   v       -- A TU-game of length 2^n-1.
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

N=length(v);
[~, n]=log2(N);

vp=zeros(1,n);
S=1:N-1;

parfor k=1:n
  Si=S(bitget(S,k)==0);
  sp=v(Si)==0;
  if all(sp==1)
    vp(k)=k;	  
    else
  end
end 
vp_vec=vp>0;
