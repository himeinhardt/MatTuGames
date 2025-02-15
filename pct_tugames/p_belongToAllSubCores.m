function bacQ=p_belongToAllSubCores(v,x,tol)
% P_BELONGTOALLSUBCORES checks whether all projections of an imputation x are a member of 
% an associated core of a subgame using MATLAB's PCT.
%
%
% Usage: bacQ=p_belongToAllSubCores(v,x,tol)
%
% Define structure variables:
%  output:
%  Q         -- If Q=1 all projection of x are a member of an associated core 
%               of a subgame, otherwise false.
%  sbcrQ     -- An array of ones (true) and/or false indicating if the 
%               the projection xS belongs to the core of the associated subgame vS. 
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1.
%  x        -- payoff vector of size(1,n).
%  tol      -- Tolerance value. Its default value is set to 10^8*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/04/2021        1.9.1           hme
%

if nargin<3
 tol=10^8*eps; % Change this value if the solution is not correct.
end

N=length(v);
[~, n]=log2(N);
bcQ=false(1,N);
k=1:n;

parfor S=1:N
    vS=SubGame(v,S); 
    xS=x(logical(bitget(S,k))); % projection of x on S.
    bcQ(S)=belongToCoreQ(vS,xS,'rat',tol);
end

bacQ.Q=all(bcQ);
bacQ.sbcrQ=bcQ;
