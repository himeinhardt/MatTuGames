function [uv,mv]=AntiUtopiaPayoff(v)
% ANTIUTOPIAPAYOFF computes the anti-utopia and minimum claim vector of game v.
%
% Usage: [uv mv]=AntiUtopiaPayoff(v)
% Define variables:
%  output:
%  uv       -- The anti-upper vector/payoff of game v.
%  mv       -- The agreement vector of game v.
%
%  input:
%  v        -- A TU-game of length 2^n-1.
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/18/2015        0.7             hme
%                

% Use with care, the result of mv might not be correct!!!

N=length(v);
[~, n]=log2(N);
% anti-upper vector
k=1:n;
Si=2.^(k-1);
uv=v(Si);

UV=uv(1); for k=2:n,UV=[UV uv(k) UV+uv(k)]; end

mv=zeros(1,n);
S=1:N;
for k=1:n
  Swk=S(bitget(S,k)==1);
  Snk=Swk-2^(k-1);
  if Snk(1)==0
     Snk(1)=[];
    av=[0,UV(Snk)];
  else
     av=UV(Snk);
  end
  mv(k)=min(v(Swk)-av);
end


