function [cQ,wpn_v,nwex,wexQ]=WSysNuclQ(v,pS,tol)
% WSYSPRENUCLQ checks the correctness of the nucleolus of a weight system pS.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
%
% Source: Proposition 2.2 (Kleppe, J.; Reijnierse, J.H.; Sudhölter, P., 2013)
%
% Usage: [cQ,wpn_v,nwex,wexQ]=WSysPreNuclQ(v,pS,tol) 
%
% Define variables:
%  output:
%  cQ       -- Returns true (1) whenever the weighted nucleolus w.r.t.
%              the Tu-game v and a weight system pS was correctly evaluated, 
%              otherwise false (0).
%  wpn_v    -- Returns the weighted nucleolus w.r.t. v and pS
%  nwex     -- Returns the nucleolus of the excess game wex. 
%  wexQ     -- Returns true (1) if the pre-nucleolus of the excess game
%              is a null vector, otherwise false (0).
%  
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  pS       -- A vector of weights of length 2^n-1.
%  tol      -- Specifies a tolerance value, default is 10^7*eps.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/24/2014        0.5               hme
%   01/27/2024        1.9.1             hme
%
%
narginchk(2,3); % check for legal number of input arguments.


if nargin<3
  N=length(v);
  [~, n]=log2(N);
  if (2^n-1)~=N
    error('Game has not the correct size!');
  end
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
  tol=10^7*eps;
else
  N=length(v);
  [~, n]=log2(N);
  if (2^n-1)~=N
    error('Game has not the correct size!');
  end
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
end

try
  wpn_v=cplex_weightedNucl(v,pS);
  wex=weightedExcess(v,wpn_v,pS);
  nwex=cplex_nucl(wex);
  zv=zeros(1,n);
  zex=excess(wex,zv);
  cQ=all(all(abs(zv-nwex)<tol));
  wexQ=all(all(abs(zex-wex)<tol));
catch
  wpn_v=weightedNucl(v,pS);
  wex=weightedExcess(v,wpn_v,pS);
  nwex=Nucl(wex);
  zv=zeros(1,n);
  zex=excess(wex,zv);
  cQ=all(all(abs(zv-nwex)<tol));
  wexQ=all(all(abs(zex-wex)<tol));
end
