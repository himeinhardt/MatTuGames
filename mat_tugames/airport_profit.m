function v=airport_profit(C,b)
% AIRPORT_PROFIT computes from a cost and benefit vector (C,b) the corresponding 
% airport profit game. 
%
% Usage: v=airport_profit(C,b)
%
%
% Define variables:
%  output:
%  v        -- The profit game v of length 2^n-1.
%
%  input:
%  c        -- A cost vector of length n.
%  b        -- A benefit vector of length n.
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

if nargin < 2
   error('An airport profit problem must be specified!'); 
end    
    
n=length(b);
nc=length(C);
if n~=nc
   error('The cost and benefit vector must have equal size!'); 
end    
N=2^n-1;
k=1:n;
v=zeros(1,N);
S=1:N;
gb=additive_game(b);
for k=1:n, A(:,k) = bitget(S,k)==1;end

for ss=1:N
    gC(ss)=max(C(A(ss,:)));
    sC=SubSets(ss,n);
    v(ss)=max([gb(sC)-gC(sC),0]); % Zero is the empty set.
end


    
