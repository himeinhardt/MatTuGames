function [cQ,wpn_v,nwex,wexQ]=WSysNuclQ(clv,pS,tol)
% WSYSPRENUCLQ checks the correctness of the nucleolus of a weight system pS.
% Requires Matlab's Optimization toolbox (default), otherwise CPLEX.
%
% Source: Proposition 2.2 (Kleppe, J.; Reijnierse, J.H.; Sudhölter, P., 2013)
%
% Usage: [cQ,wpn_v,nwex,wexQ]=clv.WSysPreNuclQ(pS,tol) 
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
%  clv      -- TuGame class object.
%  pS       -- A vector of weights of length 2^n-1.
%  tol      -- Specifies a tolerance value, default is 10^7*eps.

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/24/2014        0.5             hme
%
%
if nargin<3
  N=clv.tusize;
  n=clv.tuplayers;
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
  tol=10^7*eps;
else
  N=clv.tusize;
  n=clv.tuplayers;
  np=length(pS);
  if np~=N
    error('Weight vector has not the correct size!');
  end
end

try
  wpn_v=clv.cplex_weightedNucl(pS);
  wex=clv.weightedExcess(wpn_v,pS);
  nwex=cplex_nucl(wex);
  zv=zeros(1,n);
  zex=excess(wex,zv);
  cQ=all(all(abs(zv-nwex)<tol));
  wexQ=all(all(abs(zex-wex)<tol));
catch
  wpn_v=clv.weightedNucl(pS);
  wex=clv.weightedExcess(wpn_v,pS);
  nwex=Nucl(wex);
  zv=zeros(1,n);
  zex=excess(wex,zv);
  cQ=all(all(abs(zv-nwex)<tol));
  wexQ=all(all(abs(zex-wex)<tol));
end
