function delvT=diffOperator(v,T)
% DIFFOPERATOR computes the difference operator of the game w.r.t. coalition T.
%
% Source: Cores of Convex Games, By Lloyd S. Shapley, 1971.
%
% Usage: sdQ=diffOperator(v,T)
% Define variables:
%  output:
%  dlvT     -- Returns the list of diffOperators w.r.t. coalition T of length 2^n-1.
%
%  output:
%  v        -- A TU-Game of length 2^n-1.
%  T        -- A coalition represented by its unique integer representation.  
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   12/27/2021        1.9.1           hme
%        
    
N=length(v);
[~, n]=log2(N);
sS=1:N;
ST=bitor(sS,T);
SwT=sS-bitand(sS,T);
v1=zeros(1,N);
v1(SwT>0)=v(SwT(SwT>0));
delvT=v(ST)-v1;
