function [v,x]=grMaxFlowGame(E,s,t,n)
% GRMAXFLOWGAME computes from a flow problem (E,s,t,n) a TU flow game.
% Cannot be used with ownership vector. This function requires the 
% grTheory - Graph Theory Toolbox that can be found here:
% 
% http://www.mathworks.com/matlabcentral/fileexchange/4266-grtheory-graph-theory-toolbox?s_tid=srchtitle
%
% Usage: v=grMaxFlowGame(E,s,t,n)
%
% Define variables:
% output:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- A matrix of optimal solution vectors for all
%              coalition S.
% input: 
%  E        -- An edge matrix of size (lx3) for the details see grMaxFlows(). 
%              The third column must be used for the capacities. The first 
%              two columns are for the edges. 
%  s        -- A positive integer to indicate the source.
%  t        -- A positive integer to indicate the sink.
%  n        -- A positive integer to indicate the number of players.
%
%


%
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   02/16/2016        0.8             hme
%

N=2^n-1;
k=1:n;
m=size(E);
x=zeros(N,m(1));
v=zeros(1,N);
% if the source in E is given by zero.
me=min(E);
if me(1)==0
   E0=E(:,1:2);
   E0=E0+s;
   c=E(:,3);
   E=[E0,c];
end

for S=1:N
    bc=k(bitget(S,k)==0);
    idx=k(ismember(k,bc)==1);
    E1=E;
    E1(idx,3)=0;
    [y,v(S)]=grMaxFlows(E1,s,t);
    x(S,:)=y';
end
