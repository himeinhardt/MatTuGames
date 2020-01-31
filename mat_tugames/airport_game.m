function [v,gC]=airport_game(C)
% AIRPORT_GAME computes from an airport problem the associated
% airport game.
%
% Usage: [v,gC]=airport_game(C)
%
%
% Define variables:
%  output:
%  v        -- The saving game v of the cost game gC of length 2^n-1.
%  gC       -- The associated cost game gC of length 2^n-1.
%
%  input:
%  c        -- A cost vector of length n.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/26/2015        0.6             hme
%
n=length(C);
N=2^n-1;
S=1:N;

for k=1:n, A(:,k) = bitget(S,k)==1;end

for ss=1:N
    gC1(ss)=max(C(A(ss,:)));
end
v=savings_game(gC1);
gC=-gC1;
    
