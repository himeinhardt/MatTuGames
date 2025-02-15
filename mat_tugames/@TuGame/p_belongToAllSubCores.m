function bacQ=p_belongToAllSubCores(clv,x,tol)
% P_BELONGTOALLSUBCORES checks whether all projections of an imputation x are a member of 
% an associated core of a subgame using MATLAB's PCT.
%
%
% Usage: bacQ=clv.p_belongToAllSubCores(x,tol)
%
% Define structure variables:
%  output:
%  Q         -- If Q=1 all projection of x are a member of an associated core 
%               of a subgame, otherwise false.
%  sbcrQ     -- An array of ones (true) and/or false indicating if the 
%               the projection xS belongs to the core of the associated subgame vS. 
%
%  input:
%  clv      -- TuGame class object.
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

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

if nargin<3
 tol=10^8*eps; % Change this value if the solution is not correct.
end

bcQ=false(1,N);
k=1:n;

parfor S=1:N
    vS=SubGame(v,S); 
    xS=x(logical(bitget(S,k))); % projection of x on S.
    bcQ(S)=belongToCoreQ(vS,xS,'rat',tol);
end

bacQ.Q=all(bcQ);
bacQ.sbcrQ=bcQ;
