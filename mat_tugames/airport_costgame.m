function [gC,v]=airport_costgame(C,nj)
% AIRPORT_COSTGAME computes from an airport problem the associated
% airport cost game.
%
% Usage: [gC,v]=airport_costgame(C,nj)
%
%
% Define variables:
%  output:
%  gC       -- The associated cost game gC of length 2^n-1.
%  v        -- The saving game v of the cost game gC of length 2^n-1.
%
%  input:
%  c        -- A cost vector of length n.
%  nj       -- Number of airplane movements/number of players.
%

%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/04/2023        1.9.1             hme
%
n=sum(nj); % Determining the number of players.
if n<=48
% Expanding cost vector from types to players.
  exc=[];	
  for k=1:length(nj) 
	  exc=[exc,ones(1,nj(k)).*C(k)];  
  end
  C=exc;
  n1=length(C);
  if n~=n1
     error("Game is not consistent!");
     gC=[];
     v=[];
     return
  end	  
  N=2^n-1;
  S=1:N;
  A=false(N,n);
  for k=1:n, A(:,k) = bitget(S,k)==1;end
  gC=zeros(1,N);
  for ss=1:N
      gC(ss)=max(C(A(ss,:)));
  end
  v=savings_game(gC);
%  gC=-gC1;
else
  warning("The game is too large! The number of players should not exceed 48!");
  v=[];
  gC=[];
end	
