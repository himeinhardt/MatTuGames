function sdQ=secDiffOperatorQ(v,tol)
% SECDIFFOPERATOR checks whether the second order differences are all positive which indicates convexity.
%
% Source: Cores of Convex Games, By Lloyd S. Shapley, 1971.
%
% Usage: sdQ=secDiffOperatorQ(v)
% Define variables:
%  output:
%  sdQ      -- Returns one (true) whenever all second order differences are positive.
%
%  output:
%  v        -- A TU-Game of length 2^n-1.
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
%   12/27/2021        1.9.1           hme
%    

if nargin<2
   tol=10^8*eps;
end	


N=length(v);
[~, n]=log2(N);
SDQ=false(N,N);
for tT=1:N; 
    dvT=diffOperator(v,tT);
    for sS=1:N
        SDQ(tT,sS)=all(diffOperator(dvT,sS)+tol>=0); 
    end        
end
%SDQ
sdQ=all(SDQ(:));


